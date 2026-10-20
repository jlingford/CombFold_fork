"""Building assembled models from the assembler's transformations."""

import logging
from math import cos, sin
from pathlib import Path

import gemmi

log = logging.getLogger(__name__)


def euler_to_transform(numbers: list[float]) -> gemmi.Transform:
    """Rotation (rx, ry, rz Euler angles) + translation, with the same convention as gamb::RigidTrans."""
    x, y, z, tx, ty, tz = numbers
    cx, cy, cz = cos(x), cos(y), cos(z)
    sx, sy, sz = sin(x), sin(y), sin(z)
    rotation = gemmi.Mat33([
        [cz * cy, -sy * sx * cz - sz * cx, -sy * cx * cz + sz * sx],
        [sz * cy, -sy * sx * sz + cx * cz, -sy * cx * sz - sx * cz],
        [sy, cy * sx, cy * cx],
    ])
    return gemmi.Transform(rotation, gemmi.Vec3(tx, ty, tz))


def build_model(transforms: list[tuple[int, list[float]]], chain_files: list[Path], name: str) -> gemmi.Structure:
    """Place each representative chain (indexed as in chain.list) with its transformation."""
    model = gemmi.Model(1)
    for index, numbers in sorted(transforms):
        st = gemmi.read_structure(str(chain_files[index]))
        st[0].transform_pos_and_adp(euler_to_transform(numbers))
        for chain in st[0]:
            model.add_chain(chain)
    st = gemmi.Structure()
    st.name = name
    st.add_model(model)
    st.setup_entities()
    return st


def write_cif(st: gemmi.Structure, path: Path) -> None:
    st.make_mmcif_document().write_file(str(path))
