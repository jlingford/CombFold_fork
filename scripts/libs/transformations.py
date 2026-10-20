"""Pairwise rigid transformations between subunits, extracted from the input predictions with AF2trans."""

import logging
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import gemmi
import numpy as np
import polars as pl

from libs.structures import ChainRecord, ca_coords_and_bfactors, write_chain_pdb
from libs.subunits import Subunit

log = logging.getLogger(__name__)

# Same interface definition as the original CombFold: CA-CA distance below 8 A
INTERFACE_CA_DIST = 8.0


@dataclass
class PairTask:
    subunit1: Subunit
    subunit2: Subunit
    record1: ChainRecord
    record2: ChainRecord
    sample1: Path
    sample2: Path


@dataclass
class Transformation:
    subunit1: str
    subunit2: str
    source_file: Path
    source_chains: tuple[str, str]
    score: float
    rmsd1: float
    rmsd2: float
    trans: str  # "rx ry rz tx ty tz" as printed by AF2trans


def rep_pdb_name(subunit: Subunit) -> str:
    return f"{subunit.name}_{subunit.chain_id}.pdb"


def write_representatives(subunits: list[Subunit], structures: dict[Path, gemmi.Structure], rep_dir: Path) -> None:
    for subunit in subunits:
        rep = subunit.representative
        log.info("Representative of %s: %s chain %s (mean pLDDT %.1f)",
                 subunit.name, rep.file.name, rep.chain, rep.mean_plddt)
        write_chain_pdb(structures[rep.file][0][rep.chain], subunit.chain_id, rep_dir / rep_pdb_name(subunit))


def interface_score(chain1: gemmi.Chain, chain2: gemmi.Chain) -> float | None:
    """Mean pLDDT of interface CAs from both chains (the original CombFold score); None if no interface."""
    coords1, b1 = ca_coords_and_bfactors(chain1)
    coords2, b2 = ca_coords_and_bfactors(chain2)
    in_interface1 = np.zeros(len(coords1), dtype=bool)
    in_interface2 = np.zeros(len(coords2), dtype=bool)
    # Chunked to keep the distance matrix small for multi-thousand residue chains
    for start in range(0, len(coords1), 512):
        d2 = ((coords1[start:start + 512, None, :] - coords2[None, :, :]) ** 2).sum(axis=-1)
        close = d2 < INTERFACE_CA_DIST ** 2
        in_interface1[start:start + 512] |= close.any(axis=1)
        in_interface2 |= close.any(axis=0)
    if not in_interface1.any():
        return None
    return float(np.concatenate([b1[in_interface1], b2[in_interface2]]).mean())


def run_af2trans(af2trans_bin: Path, rep_dir: Path, task: PairTask) -> tuple[float, float, str]:
    """AF2trans superimposes each representative onto its copy in the prediction and returns the
    transformation of subunit2's representative relative to subunit1's representative."""
    out = subprocess.run([str(af2trans_bin),
                          str(rep_dir / rep_pdb_name(task.subunit1)), str(rep_dir / rep_pdb_name(task.subunit2)),
                          str(task.sample1), str(task.sample2)],
                         check=True, capture_output=True, text=True).stdout.strip()
    fields = out.split(" | ")
    if len(fields) != 4:
        raise RuntimeError(f"Unexpected AF2trans output: {out!r}")
    _, desc1, desc2, trans = fields
    return float(desc1.split("_")[0]), float(desc2.split("_")[0]), trans


def compute_transformation(task: PairTask, structures: dict[Path, gemmi.Structure], af2trans_bin: Path,
                           rep_dir: Path) -> Transformation | None:
    model = structures[task.record1.file][0]
    score = interface_score(model[task.record1.chain], model[task.record2.chain])
    if score is None:
        log.info("%s: no interface between %s and %s, skipping",
                 task.record1.file.name, task.subunit1.name, task.subunit2.name)
        return None
    rmsd1, rmsd2, trans = run_af2trans(af2trans_bin, rep_dir, task)
    return Transformation(task.subunit1.name, task.subunit2.name, task.record1.file,
                          (task.record1.chain, task.record2.chain), score, rmsd1, rmsd2, trans)


def build_pair_tasks(subunits: list[Subunit], structures: dict[Path, gemmi.Structure],
                     sample_dir: Path) -> list[PairTask]:
    """One task per pair of subunits co-occurring in an input file; also writes the sample chains."""
    file_index = {path: i for i, path in enumerate(structures)}
    by_file: dict[Path, list[tuple[Subunit, ChainRecord, Path]]] = defaultdict(list)
    for subunit in subunits:
        for record in subunit.occurrences:
            # Index prefix: AF3 job folders often reuse the same file name (e.g. model.cif)
            sample = sample_dir / f"{file_index[record.file]:04d}_{record.file.name.split('.')[0]}_{record.chain}.pdb"
            write_chain_pdb(structures[record.file][0][record.chain], subunit.chain_id, sample)
            by_file[record.file].append((subunit, record, sample))
    return [PairTask(s1, s2, r1, r2, p1, p2)
            for members in by_file.values()
            for (s1, r1, p1), (s2, r2, p2) in combinations(members, 2)]


def extract_transformations(subunits: list[Subunit], structures: dict[Path, gemmi.Structure], cache_dir: Path,
                            af2trans_bin: Path, threads: int) -> list[Transformation]:
    rep_dir = cache_dir / "assembly"
    sample_dir = cache_dir / "samples"
    sample_dir.mkdir(parents=True, exist_ok=True)
    tasks = build_pair_tasks(subunits, structures, sample_dir)
    log.info("Extracting %d pairwise transformations with %d threads", len(tasks), threads)
    # Threads suffice: the heavy lifting happens in the AF2trans subprocess
    with ThreadPoolExecutor(max_workers=threads) as pool:
        results = pool.map(lambda t: compute_transformation(t, structures, af2trans_bin, rep_dir), tasks)
        return [r for r in results if r is not None]


def write_transformation_files(transformations: list[Transformation], subunits: list[Subunit],
                               trans_dir: Path) -> None:
    """One file per subunit pair, best interface score first, in the format read by the assembler
    (BBContainer.cc): 'index | score | description | rx ry rz tx ty tz'."""
    trans_dir.mkdir(parents=True, exist_ok=True)
    chain_ids = {s.name: s.chain_id for s in subunits}
    by_pair: dict[tuple[str, str], list[Transformation]] = defaultdict(list)
    for t in transformations:
        by_pair[(t.subunit1, t.subunit2)].append(t)
    for (name1, name2), pair_transformations in by_pair.items():
        pair_transformations.sort(key=lambda t: t.score, reverse=True)
        log.info("%d transformations between %s and %s", len(pair_transformations), name1, name2)
        lines = []
        for i, t in enumerate(pair_transformations, start=1):
            # The assembler splits on ':', '|' and tabs, so keep them out of the description
            desc = f"{t.rmsd1}_{t.rmsd2}_{t.source_chains[0]}{t.source_chains[1]}_{t.source_file.name}"
            desc = desc.translate(str.maketrans({":": "-", "|": "-", "\t": "-", " ": "-"}))
            lines.append(f"{i} | {t.score} | {desc} | {t.trans}\n")
        (trans_dir / f"{name1}_{chain_ids[name1]}_plus_{name2}_{chain_ids[name2]}").write_text("".join(lines))


def write_transformations_tsv(transformations: list[Transformation], path: Path) -> None:
    rows = [{"subunit1": t.subunit1,
             "subunit2": t.subunit2,
             "source_file": str(t.source_file),
             "source_chain1": t.source_chains[0],
             "source_chain2": t.source_chains[1],
             "interface_plddt": round(t.score, 2),
             "rep_rmsd1": t.rmsd1,
             "rep_rmsd2": t.rmsd2}
            for t in transformations]
    schema = {"subunit1": pl.String, "subunit2": pl.String, "source_file": pl.String, "source_chain1": pl.String,
              "source_chain2": pl.String, "interface_plddt": pl.Float64, "rep_rmsd1": pl.Float64,
              "rep_rmsd2": pl.Float64}
    pl.DataFrame(rows, schema=schema).sort(["subunit1", "subunit2", "interface_plddt"],
                                           descending=[False, False, True]).write_csv(path, separator="\t")
