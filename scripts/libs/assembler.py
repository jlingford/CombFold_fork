"""Running the C++ CombinatorialAssembler and parsing its .res output."""

import argparse
import logging
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from libs.subunits import Subunit
from libs.transformations import rep_pdb_name

log = logging.getLogger(__name__)

DEFAULT_BIN_DIR = Path(__file__).resolve().parents[2] / "CombinatorialAssembler"
# multPen_, singlePen_, diffPen and maxPen_ are also printed in .res lines but never assigned in SuperBB.cc
# (uninitialised values), so they are not reported
RES_SCORE_KEYS = ["size", "weightedTransScore", "transScore", "backBonePen", "restraintsRatio"]


@dataclass
class AssemblerParams:
    trans_num: int  # max transformations read per subunit pair
    best_k: int  # assemblies kept per subset of subunits at each step
    max_backbone_collision: float  # max fraction of a chain's backbone allowed to clash (-b)
    min_plddt: float  # atoms below this pLDDT are ignored in clash checks (-t)


def add_assembler_args(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("assembler parameters (defaults as in the original CombFold)")
    group.add_argument("--trans-num", type=int, default=900, help="max transformations read per subunit pair")
    group.add_argument("--best-k", type=int, default=100, help="assemblies kept per subunit subset")
    group.add_argument("--max-backbone-collision", type=float, default=0.05,
                       help="max fraction of a chain's backbone allowed to clash")
    group.add_argument("--min-plddt", type=float, default=80,
                       help="atoms with lower pLDDT are ignored when checking clashes")
    group.add_argument("--bin-dir", type=Path, default=DEFAULT_BIN_DIR,
                       help="directory with CombinatorialAssembler.out and AF2trans.out")


def params_from_args(args: argparse.Namespace) -> AssemblerParams:
    return AssemblerParams(args.trans_num, args.best_k, args.max_backbone_collision, args.min_plddt)


def find_binary(bin_dir: Path, name: str) -> Path:
    path = bin_dir / name
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found; build it with `make -C {bin_dir}`")
    return path


def write_chain_list(subunits: list[Subunit], rep_dir: Path) -> None:
    """Line order defines the subunit indices used in the .res output."""
    (rep_dir / "chain.list").write_text("".join(f"{rep_pdb_name(s)}\n" for s in subunits))


def run_assembler(bin_dir: Path, rep_dir: Path, trans_dir: Path, params: AssemblerParams,
                  crosslinks: Path | None) -> Path:
    """Run the assembler in rep_dir and return the path of the result file to use."""
    xlinks = rep_dir / "xlink_consts.txt"
    if crosslinks is not None:
        shutil.copy(crosslinks, xlinks)
    else:
        xlinks.write_text("")
    cmd = [str(find_binary(bin_dir, "CombinatorialAssembler.out")), "chain.list", f"{trans_dir.resolve()}/",
           str(params.trans_num), str(params.best_k), xlinks.name,
           "-b", str(params.max_backbone_collision), "-t", str(params.min_plddt)]
    log.info("Running assembler: %s", " ".join(cmd))
    with open(rep_dir / "output.log", "w") as log_file:
        result = subprocess.run(cmd, cwd=rep_dir, stdout=log_file, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        tail = (rep_dir / "output.log").read_text().splitlines()[-20:]
        raise RuntimeError(f"Assembler failed (exit {result.returncode}), see {rep_dir / 'output.log'}:\n"
                           + "\n".join(tail))
    return find_result_file(rep_dir)


def find_result_file(rep_dir: Path) -> Path:
    clustered = rep_dir / "output_clustered.res"
    if clustered.is_file() and clustered.stat().st_size > 0:
        return clustered
    # Without a complete assembly the assembler only writes the largest partial ones (cb_<n>_output.res)
    partial = sorted(rep_dir.glob("cb_*_output.res"), key=lambda p: int(p.name.split("_")[1]))
    if not partial:
        raise RuntimeError(f"Assembler produced no results, see {rep_dir / 'output.log'}")
    log.warning("No complete assembly found; using partial assemblies from %s", partial[-1].name)
    return partial[-1]


def parse_res(path: Path) -> pl.DataFrame:
    """One row per assembly. A line looks like:
    size_ 5 transScore_ 335.5 ... weightedTransScore 84.2 [1(0 0 0 0 0 0),2(rx ry rz tx ty tz),...]  ..."""
    rows = []
    for line in path.read_text().splitlines():
        if "[" not in line:
            continue
        head, rest = line.split("[", 1)
        tokens = head.split()
        values = {k.rstrip("_"): float(v) for k, v in zip(tokens[::2], tokens[1::2])}
        row = {key: values[key] for key in RES_SCORE_KEYS}
        row["size"] = int(row["size"])
        row["transforms"] = rest.split("]", 1)[0]
        rows.append(row)
    return pl.DataFrame(rows).with_row_index("rank")


def parse_transforms(transforms: str) -> list[tuple[int, list[float]]]:
    """'0(rx ry rz tx ty tz),1(...)' -> [(chain.list index, [6 numbers]), ...]"""
    return [(int(idx), [float(x) for x in nums.split()])
            for idx, nums in re.findall(r"(\d+)\(([^)]*)\)", transforms)]
