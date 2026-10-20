"""Top-level steps shared by the CLI scripts."""

import logging
import shutil
from pathlib import Path

import gemmi
import polars as pl

from libs import assembler, reconstruct, subunits as su, transformations as tr
from libs.structures import ChainRecord, protein_chains, read_structure

log = logging.getLogger(__name__)

CACHE = "cache"
DETECT_OUTPUTS = ["subunits.json", f"{CACHE}/combfold_subunits.tsv"]
# Listed individually so that assemble.py keeps the detection outputs it is run on
ASSEMBLY_OUTPUTS = ["combfold_scores.tsv", "combfold_model_*.cif", f"{CACHE}/combfold_transformations.tsv",
                    f"{CACHE}/assembly", f"{CACHE}/samples", f"{CACHE}/transformations"]


def prepare_out_dir(out_dir: Path, outputs: list[str], overwrite: bool) -> None:
    """Refuse to clobber earlier results unless asked; only ever delete files this pipeline writes."""
    existing = [path for pattern in outputs for path in sorted(out_dir.glob(pattern))]
    if existing and not overwrite:
        raise FileExistsError(f"Outputs already exist in {out_dir} (use --overwrite): {[str(p.relative_to(out_dir)) for p in existing]}")
    for path in existing:
        shutil.rmtree(path) if path.is_dir() else path.unlink()
    (out_dir / CACHE).mkdir(parents=True, exist_ok=True)


def load_structures(paths: list[Path]) -> tuple[dict[Path, gemmi.Structure], list[ChainRecord]]:
    paths = sorted({p.resolve() for p in paths})
    if len(paths) < 2:
        raise ValueError("At least two input structures are needed")
    structures = {p: read_structure(p) for p in paths}
    records = [r for p, st in structures.items() for r in protein_chains(p, st)]
    log.info("Read %d protein chains from %d structures", len(records), len(paths))
    return structures, records


def detect(records: list[ChainRecord], out_dir: Path, names_fasta: Path | None) -> list[su.Subunit]:
    names = su.read_names_fasta(names_fasta) if names_fasta else None
    subunits = su.define_subunits(records, names)
    su.check_connected(subunits, su.cooccurrence_edges(subunits), "shared input structures")
    su.write_subunits_json(subunits, out_dir / "subunits.json")
    su.write_subunits_tsv(subunits, out_dir / CACHE / "combfold_subunits.tsv")
    for s in subunits:
        log.info("  %s (chain %s, %d aa) in %d structure(s)", s.name, s.chain_id, len(s.sequence), len(s.occurrences))
    return subunits


def assemble(subunits: list[su.Subunit], structures: dict[Path, gemmi.Structure], out_dir: Path,
             params: assembler.AssemblerParams, bin_dir: Path, crosslinks: Path | None, top_n: int,
             threads: int) -> pl.DataFrame:
    af2trans_bin = assembler.find_binary(bin_dir, "AF2trans.out")
    assembler.find_binary(bin_dir, "CombinatorialAssembler.out")
    cache_dir = out_dir / CACHE
    rep_dir = cache_dir / "assembly"
    trans_dir = cache_dir / "transformations"
    rep_dir.mkdir(parents=True)
    su.write_subunits_tsv(subunits, cache_dir / "combfold_subunits.tsv")

    tr.write_representatives(subunits, structures, rep_dir)
    transformations = tr.extract_transformations(subunits, structures, cache_dir, af2trans_bin, threads)
    tr.write_transformations_tsv(transformations, cache_dir / "combfold_transformations.tsv")
    # Pairs without an interface give no transformation, which can disconnect the subunit graph
    su.check_connected(subunits, [(t.subunit1, t.subunit2) for t in transformations], "predicted interfaces")
    tr.write_transformation_files(transformations, subunits, trans_dir)

    assembler.write_chain_list(subunits, rep_dir)
    res_path = assembler.run_assembler(bin_dir, rep_dir, trans_dir, params, crosslinks)
    return write_models(res_path, subunits, rep_dir, out_dir, top_n)


def write_models(res_path: Path, subunits: list[su.Subunit], rep_dir: Path, out_dir: Path,
                 top_n: int) -> pl.DataFrame:
    chain_files = [rep_dir / tr.rep_pdb_name(s) for s in subunits]
    scores = assembler.parse_res(res_path).head(top_n)
    model_files = []
    for row in scores.iter_rows(named=True):
        model_path = out_dir / f"combfold_model_{row['rank']}.cif"
        st = reconstruct.build_model(assembler.parse_transforms(row["transforms"]), chain_files, model_path.stem)
        reconstruct.write_cif(st, model_path)
        model_files.append(str(model_path))
    scores = (scores.with_columns(model=pl.Series(model_files),
                                  n_subunits_total=pl.lit(len(subunits)))
              .rename({"size": "n_subunits_placed"})
              .select("rank", "model", "weightedTransScore", "n_subunits_placed", "n_subunits_total",
                      "transScore", "backBonePen", "restraintsRatio"))
    scores.write_csv(out_dir / "combfold_scores.tsv", separator="\t")
    log.info("Wrote %d models to %s; best weightedTransScore %.2f",
             len(scores), out_dir, scores["weightedTransScore"].max())
    return scores
