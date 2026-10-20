#!/usr/bin/env python3
"""Detect the subunits (unique protein chains) shared across overlapping sub-complex predictions.

Writes subunits.json (CombFold schema, hand-editable, input of assemble.py) and cache/combfold_subunits.tsv.
"""

import argparse
import logging
from pathlib import Path

from libs import pipeline


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("structures", type=Path, nargs="+", help="predicted structures (.cif/.pdb, optionally .gz)")
    parser.add_argument("-o", "--out-dir", type=Path, required=True)
    parser.add_argument("--names", type=Path, help="FASTA naming subunits by sequence (header = name)")
    parser.add_argument("--overwrite", action="store_true", help="replace existing outputs")
    return parser.parse_args()


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = parse_args()
    pipeline.prepare_out_dir(args.out_dir, pipeline.DETECT_OUTPUTS, args.overwrite)
    _, records = pipeline.load_structures(args.structures)
    pipeline.detect(records, args.out_dir, args.names)


if __name__ == "__main__":
    main()
