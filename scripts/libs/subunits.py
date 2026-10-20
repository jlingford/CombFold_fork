"""Auto-detection of subunits (unique protein chains) across overlapping sub-complex predictions.

Current scope: heteromers only. Each protein occurs once in the complex and once per input file;
inputs overlap by sharing whole chains with identical sequences.
"""

import json
import logging
import re
import string
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import polars as pl

from libs.structures import ChainRecord

log = logging.getLogger(__name__)

# The assembler stores chain IDs as a single char
CHAIN_IDS = string.ascii_uppercase + string.ascii_lowercase + string.digits
# Subunit names end up in file names that the assembler splits on "_plus_"
VALID_NAME = re.compile(r"^[A-Za-z0-9.\-]+$")


@dataclass
class Subunit:
    name: str
    chain_id: str
    sequence: str
    occurrences: list[ChainRecord] = field(default_factory=list)

    @property
    def representative(self) -> ChainRecord:
        """The occurrence with the highest mean pLDDT is used as the rigid body for assembly."""
        return max(self.occurrences, key=lambda r: r.mean_plddt)


def check_no_homomers(records: list[ChainRecord]) -> None:
    seen: dict[tuple[Path, str], str] = {}
    for record in records:
        key = (record.file, record.sequence)
        if key in seen:
            raise ValueError(f"{record.file}: chains {seen[key]} and {record.chain} have identical sequences. "
                             f"Homo-oligomers are not supported yet.")
        seen[key] = record.chain


def read_names_fasta(path: Path) -> dict[str, str]:
    """Map sequence -> name from a FASTA file (name = first word of the header)."""
    names = {}
    for entry in gemmi.read_pir_or_fasta(path.read_text()):
        name = entry.header.split()[0]
        if not VALID_NAME.match(name):
            raise ValueError(f"Invalid subunit name '{name}' in {path}: use only letters, digits, '.' and '-'")
        names[entry.seq.replace("\n", "").upper()] = name
    return names


def define_subunits(records: list[ChainRecord], names: dict[str, str] | None = None) -> list[Subunit]:
    """One subunit per unique sequence, ordered by first appearance (records are sorted by file)."""
    check_no_homomers(records)
    names = names or {}
    by_sequence: dict[str, Subunit] = {}
    for record in records:
        if record.sequence not in by_sequence:
            i = len(by_sequence)
            if i >= len(CHAIN_IDS):
                raise ValueError(f"More than {len(CHAIN_IDS)} unique chains; the assembler needs 1-char chain IDs")
            name = names.get(record.sequence, f"S{i + 1}")
            by_sequence[record.sequence] = Subunit(name, CHAIN_IDS[i], record.sequence)
        by_sequence[record.sequence].occurrences.append(record)

    subunits = list(by_sequence.values())
    duplicated = {s.name for s in subunits if [t.name for t in subunits].count(s.name) > 1}
    if duplicated:
        raise ValueError(f"Duplicated subunit names: {sorted(duplicated)}")
    log.info("Detected %d subunits in %d chains", len(subunits), len(records))
    return subunits


def assign_occurrences(subunits: list[Subunit], records: list[ChainRecord]) -> list[Subunit]:
    """Attach input chains to predefined subunits (e.g. read from subunits.json) by exact sequence."""
    check_no_homomers(records)
    by_sequence = {s.sequence: s for s in subunits}
    for subunit in subunits:
        subunit.occurrences = []
    for record in records:
        if record.sequence not in by_sequence:
            log.warning("%s chain %s matches no subunit, ignoring it", record.file.name, record.chain)
            continue
        by_sequence[record.sequence].occurrences.append(record)
    missing = [s.name for s in subunits if not s.occurrences]
    if missing:
        raise ValueError(f"Subunits not found in any input structure: {missing}")
    return subunits


def connected_components(nodes: list[str], edges: list[tuple[str, str]]) -> list[set[str]]:
    neighbours: dict[str, set[str]] = defaultdict(set)
    for a, b in edges:
        neighbours[a].add(b)
        neighbours[b].add(a)
    components, seen = [], set()
    for node in nodes:
        if node in seen:
            continue
        component, stack = set(), [node]
        while stack:
            current = stack.pop()
            if current not in component:
                component.add(current)
                stack.extend(neighbours[current] - component)
        seen |= component
        components.append(component)
    return components


def check_connected(subunits: list[Subunit], edges: list[tuple[str, str]], what: str) -> None:
    """The assembler aborts if the subunit graph is disconnected, so fail early with a readable message."""
    components = connected_components([s.name for s in subunits], edges)
    if len(components) == 1:
        return
    by_name = {s.name: s for s in subunits}
    lines = []
    for i, component in enumerate(components, start=1):
        files = sorted({r.file.name for name in component for r in by_name[name].occurrences})
        lines.append(f"  component {i}: subunits {sorted(component)} from {files}")
    raise ValueError(f"Subunits are not connected by {what}, so they cannot be assembled into one complex:\n"
                     + "\n".join(lines)
                     + "\nCheck that overlapping inputs share a chain with an *identical* sequence.")


def cooccurrence_edges(subunits: list[Subunit]) -> list[tuple[str, str]]:
    names_by_file: dict[Path, list[str]] = defaultdict(list)
    for subunit in subunits:
        for record in subunit.occurrences:
            names_by_file[record.file].append(subunit.name)
    return [(names[0], other) for names in names_by_file.values() for other in names[1:]]


def write_subunits_json(subunits: list[Subunit], path: Path) -> None:
    """Original CombFold subunits.json schema, so the file stays hand-editable."""
    data = {s.name: {"name": s.name, "chain_names": [s.chain_id], "start_res": 1, "sequence": s.sequence}
            for s in subunits}
    path.write_text(json.dumps(data, indent=2) + "\n")


def read_subunits_json(path: Path) -> list[Subunit]:
    subunits = []
    for key, entry in json.loads(path.read_text()).items():
        if key != entry["name"]:
            raise ValueError(f"{path}: key '{key}' does not match name '{entry['name']}'")
        if not VALID_NAME.match(key):
            raise ValueError(f"{path}: invalid subunit name '{key}': use only letters, digits, '.' and '-'")
        if len(entry["chain_names"]) != 1 or entry.get("start_res", 1) != 1:
            raise ValueError(f"{path}: subunit '{key}' has several copies or start_res != 1, "
                             f"which is not supported yet (heteromers of whole chains only)")
        subunits.append(Subunit(key, entry["chain_names"][0], entry["sequence"].upper()))
    chain_ids = [s.chain_id for s in subunits]
    if len(set(chain_ids)) != len(chain_ids) or any(len(c) != 1 for c in chain_ids):
        raise ValueError(f"{path}: chain names must be unique single characters")
    return subunits


def write_subunits_tsv(subunits: list[Subunit], path: Path) -> None:
    rows = [{"name": s.name,
             "chain_id": s.chain_id,
             "length": len(s.sequence),
             "file": str(r.file),
             "source_chain": r.chain,
             "mean_plddt": round(r.mean_plddt, 2),
             "is_representative": r is s.representative}
            for s in subunits for r in s.occurrences]
    pl.DataFrame(rows).write_csv(path, separator="\t")
