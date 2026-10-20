"""Structure I/O helpers built on gemmi (mmCIF/PDB, optionally gzipped)."""

import logging
from dataclasses import dataclass
from pathlib import Path

import gemmi
import numpy as np

log = logging.getLogger(__name__)

PROTEIN_TYPES = (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD)


@dataclass
class ChainRecord:
    """One protein chain in one input structure."""
    file: Path
    chain: str
    sequence: str
    mean_plddt: float


def read_structure(path: Path) -> gemmi.Structure:
    """Read the first model of a .cif/.pdb(.gz) file, without ligands, waters or hydrogens."""
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    while len(st) > 1:
        del st[1]
    st.remove_ligands_and_waters()
    st.remove_hydrogens()
    st.remove_empty_chains()
    return st


def chain_sequence(chain: gemmi.Chain) -> str:
    # Upper-casing maps modified residues (e.g. MSE -> m) onto their parent residue letter
    return gemmi.one_letter_code([res.name for res in chain]).upper()


def mean_plddt(chain: gemmi.Chain) -> float:
    """Mean CA B-factor; AF2 and AF3 both store per-residue pLDDT there."""
    _, bfactors = ca_coords_and_bfactors(chain)
    return float(bfactors.mean())


def ca_coords_and_bfactors(chain: gemmi.Chain) -> tuple[np.ndarray, np.ndarray]:
    cas = [res.find_atom("CA", "*") for res in chain]
    cas = [atom for atom in cas if atom is not None]
    coords = np.array([atom.pos.tolist() for atom in cas])
    bfactors = np.array([atom.b_iso for atom in cas])
    return coords, bfactors


def protein_chains(path: Path, st: gemmi.Structure) -> list[ChainRecord]:
    """Protein chains of a structure; other polymers (DNA/RNA) are skipped with a warning."""
    records = []
    for chain in st[0]:
        polymer_type = chain.get_polymer().check_polymer_type()
        if polymer_type not in PROTEIN_TYPES:
            log.warning("%s: skipping non-protein chain %s (%s)", path.name, chain.name, polymer_type.name)
            continue
        records.append(ChainRecord(path, chain.name, chain_sequence(chain), mean_plddt(chain)))
    return records


def write_chain_pdb(chain: gemmi.Chain, chain_id: str, out_path: Path) -> None:
    """Write a single chain as PDB (the C++ binaries only read PDB), renamed and renumbered from 1."""
    new_chain = chain.clone()
    new_chain.name = chain_id
    for i, res in enumerate(new_chain, start=1):
        res.seqid = gemmi.SeqId(i, " ")
    model = gemmi.Model(1)
    model.add_chain(new_chain)
    st = gemmi.Structure()
    st.add_model(model)
    st.setup_entities()
    st.write_pdb(str(out_path), gemmi.PdbWriteOptions(minimal=True))
