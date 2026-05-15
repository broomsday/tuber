import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
import pytest

from tuber.geometry import generate_nanotube
from tuber.structure import build_atom_array
from tuber.writers import write_structure


def _make_atom_array():
    geometry = generate_nanotube(n=3, m=3, units=2)
    return build_atom_array(geometry.coordinates, geometry.elements)


@pytest.mark.parametrize(
    ("file_format", "suffix"),
    [
        ("pdb", ".pdb"),
        ("cif", ".cif"),
    ],
)
def test_write_structure_smoke(tmp_path, file_format: str, suffix: str) -> None:
    atom_array = _make_atom_array()
    output_path = tmp_path / f"tube{suffix}"

    write_structure(atom_array, output_path, file_format)

    assert output_path.exists()
    assert output_path.stat().st_size > 0

    if file_format == "pdb":
        loaded = pdb.PDBFile.read(str(output_path)).get_structure(model=1)
    else:
        loaded = pdbx.get_structure(pdbx.CIFFile.read(str(output_path)), model=1)

    assert loaded.array_length() == atom_array.array_length()


def test_write_structure_rejects_suffix_mismatch(tmp_path) -> None:
    atom_array = _make_atom_array()
    with pytest.raises(ValueError, match="does not match requested format"):
        write_structure(atom_array, tmp_path / "tube.cif", "pdb")
