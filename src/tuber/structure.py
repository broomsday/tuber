from __future__ import annotations

from collections.abc import Sequence

import biotite.structure as struc
import numpy as np
from numpy.typing import NDArray


def build_atom_array(
    coordinates: NDArray[np.float64],
    elements: Sequence[str] | None = None,
) -> struc.AtomArray:
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("coordinates must have shape (N, 3)")

    atom_count = coordinates.shape[0]
    if elements is None:
        element_array = np.full(atom_count, "C", dtype="U2")
    else:
        element_array = np.asarray(list(elements), dtype="U2")
        if element_array.shape != (atom_count,):
            raise ValueError("elements must have one entry per coordinate")

    atom_array = struc.AtomArray(atom_count)
    atom_array.coord = coordinates
    atom_array.chain_id = np.full(atom_count, "A", dtype="U4")
    atom_array.res_id = np.ones(atom_count, dtype=int)
    atom_array.ins_code = np.full(atom_count, "", dtype="U1")
    atom_array.res_name = np.full(atom_count, "CNT", dtype="U5")
    atom_array.hetero = np.ones(atom_count, dtype=bool)
    atom_array.atom_name = np.full(atom_count, "C", dtype="U6")
    atom_array.element = element_array
    atom_array.set_annotation("atom_id", np.arange(1, atom_count + 1, dtype=int))
    atom_array.set_annotation("occupancy", np.ones(atom_count, dtype=float))
    atom_array.set_annotation("b_factor", np.zeros(atom_count, dtype=float))
    return atom_array
