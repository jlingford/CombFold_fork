#!/usr/bin/env python3
"""Assemble a complex from overlapping sub-complex predictions, given a subunits.json.

Writes combfold_model_<rank>.cif and combfold_scores.tsv; other tables and intermediate files go to cache/.
"""

import argparse
import logging
from pathlib import Path

from libs import assembler, pipeline
from libs.subunits import assign_occurrences, read_subunits_json


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("structures", type=Path, nargs="+", help="predicted structures (.cif/.pdb, optionally .gz)")
    parser.add_argument("-s", "--subunits", type=Path, required=True, help="subunits.json")
    parser.add_argument("-o", "--out-dir", type=Path, required=True)
    parser.add_argument("--crosslinks", type=Path, help="crosslinks file using the subunits.json chain names")
    parser.add_argument("--top-n", type=int, default=5, help="number of models to write")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--overwrite", action="store_true", help="replace existing outputs")
    assembler.add_assembler_args(parser)
    return parser.parse_args()


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = parse_args()
    pipeline.prepare_out_dir(args.out_dir, pipeline.ASSEMBLY_OUTPUTS, args.overwrite)
    structures, records = pipeline.load_structures(args.structures)
    subunits = assign_occurrences(read_subunits_json(args.subunits), records)
    pipeline.assemble(subunits, structures, args.out_dir, assembler.params_from_args(args), args.bin_dir,
                      args.crosslinks, args.top_n, args.threads)


if __name__ == "__main__":
    main()
