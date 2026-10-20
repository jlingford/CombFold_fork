# CombFold (assembly-only fork)

A trimmed fork of [CombFold](https://github.com/dina-lab3D/CombFold)
([paper](https://www.nature.com/articles/s41592-024-02174-0)), keeping only the
**combinatorial assembly** step. Pass it several predicted structures of
overlapping sub-complexes (e.g. AlphaFold3 `.cif` files); it detects the shared
subunits, extracts pairwise rigid transformations, runs the C++ assembler and
writes the assembled complex as mmCIF together with a score table.

Current scope: **heteromers only**. Every protein occurs once in the complex,
and inputs overlap by sharing whole chains with *identical* sequences
(e.g. `job1 = A+B`, `job2 = B+C`, `job3 = C+D`). Homo-oligomers and chains split
into fragments across predictions are not supported yet.

## Installation

```bash
micromamba create -f envs/combfold.yaml
micromamba activate combfold
make -C CombinatorialAssembler
```

The Makefile takes Boost from the active env (`$CONDA_PREFIX`). To use another
Boost install, run e.g. `make BOOST_INCLUDE=/opt/homebrew/include BOOST_LIB=/opt/homebrew/lib`.

## Usage

### One step

```bash
python scripts/combfold_auto.py job1_model.cif job2_model.cif job3_model.cif.gz -o out/ --threads 8
```

You can pass several models per job (e.g. all AF3 seeds and samples). Every
pair of chains in every input gives one candidate transformation, ranked by
interface pLDDT, and the assembler chooses among them.

### Two steps (to inspect or edit the subunits)

```bash
python scripts/detect_subunits.py inputs/*.cif -o out/ [--names names.fasta]
# optional: edit out/cache/subunits.json (names, chain IDs)
python scripts/assemble.py inputs/*.cif -s out/cache/subunits.json -o out/
```

### Useful options

- `--names names.fasta`: name the subunits by sequence. The header's first word
  becomes the name; the default names are `S1`, `S2`, …
- `--crosslinks xl.txt`: distance restraints, one per line, in the form
  `res1 chain1 res2 chain2 [minDist] maxDist [weight]`. Chain IDs are the
  detected ones (see `cache/subunits.json`), and residues are numbered from 1 within
  each chain.
- `--top-n`: how many models to write.
- Assembler parameters: `--trans-num`, `--best-k`, `--max-backbone-collision`,
  `--min-plddt`.

### HPC (SLURM)

```bash
sbatch slurm/combfold_auto.sbatch out/ jobs/*/*_model.cif
COMBFOLD_ARGS="--names names.fasta" sbatch slurm/combfold_auto.sbatch out/ jobs/*/*_model.cif
```

## Outputs

| Path | Content |
|---|---|
| `combfold_model_<rank>.cif` | Assembled models, best first (rank 0 = best). B-factors hold the representatives' pLDDT |
| `combfold_scores.tsv` | One row per model: `weightedTransScore` (the CombFold confidence, roughly the mean interface pLDDT of the transformations used), `transScore` (their sum), `n_subunits_placed`, `backBonePen`, `restraintsRatio` (fraction of crosslinks satisfied) |
| `cache/subunits.json` | Detected subunits, in the original CombFold schema |
| `cache/combfold_subunits.tsv` | Each occurrence of each subunit: source file and chain, mean pLDDT, and which occurrence is the representative |
| `cache/combfold_transformations.tsv` | Every pairwise transformation: source file and chains, interface pLDDT, representative superposition RMSDs |
| `cache/` (rest) | Intermediate files: representative and sample PDBs, transformation files, `chain.list`, raw `.res` output, `output.log` |

If no complete assembly is found, the largest partial assemblies are written
instead, with a warning; check `n_subunits_placed`.

Input requirements:
- Input structures can be `.cif` or `.pdb`, optionally gzipped.
- Non-protein chains (ligands, nucleic acids) are ignored.
- The inputs must form a connected graph through their shared subunits, and each
  pair of co-predicted subunits needs a predicted interface (CA–CA < 8 Å). If
  either condition fails, the run stops early and lists the disconnected groups.

## Layout

```
scripts/combfold_auto.py, detect_subunits.py, assemble.py   CLI entry points
scripts/libs/                 structures, subunits, transformations, assembler, reconstruct, pipeline
CombinatorialAssembler/       original C++ assembler + AF2trans
example/example_xlinks/       heteromer test set (5 proteins, all pairwise AFM dimers, crosslinks)
example/pdbs/                 original homomer example (not supported by the new scripts yet)
```
