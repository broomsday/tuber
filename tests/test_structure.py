import numpy as np

from tuber.geometry import generate_nanotube
from tuber.structure import build_atom_array


def test_build_atom_array_preserves_mixed_element_annotations() -> None:
    geometry = generate_nanotube(n=3, m=3, units=2, hydrogen_terminate=True)

    atom_array = build_atom_array(geometry.coordinates, geometry.elements)

    assert atom_array.array_length() == geometry.atom_count
    assert np.array_equal(atom_array.element, geometry.elements)
    assert np.array_equal(atom_array.atom_name, geometry.elements.astype("U6"))
