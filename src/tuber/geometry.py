from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from numpy.typing import NDArray

SQRT3 = math.sqrt(3.0)
EPSILON = 1e-9


@dataclass(frozen=True)
class NanotubeParameters:
    n: int
    m: int
    bond_length: float
    chiral_vector: NDArray[np.float64]
    translational_vector: NDArray[np.float64]
    circumference: float
    radius: float
    unit_cell_length: float
    translational_gcd: int
    t1: int
    t2: int
    unit_cell_atom_count: int


@dataclass(frozen=True)
class NanotubeGeometry:
    coordinates: NDArray[np.float64]
    elements: NDArray[np.str_]
    radius: float
    unit_cell_length: float
    total_length: float
    unit_cell_atom_count: int
    atom_count: int
    units: int
    n: int
    m: int
    bond_length: float


def graphene_lattice(
    bond_length: float = 1.42,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    if bond_length <= 0:
        raise ValueError("bond_length must be positive")

    a1 = np.array([SQRT3 * bond_length, 0.0], dtype=float)
    a2 = np.array([SQRT3 * bond_length / 2.0, 3.0 * bond_length / 2.0], dtype=float)
    basis = np.array(
        [
            [0.0, 0.0],
            (a1 + a2) / 3.0,
        ],
        dtype=float,
    )
    return a1, a2, basis


def calculate_nanotube_parameters(
    n: int,
    m: int,
    bond_length: float = 1.42,
) -> NanotubeParameters:
    _validate_chiral_indices(n, m)
    a1, a2, _ = graphene_lattice(bond_length=bond_length)

    translational_gcd = math.gcd(2 * m + n, 2 * n + m)
    t1 = (2 * m + n) // translational_gcd
    t2 = -(2 * n + m) // translational_gcd

    chiral_vector = n * a1 + m * a2
    translational_vector = t1 * a1 + t2 * a2
    circumference = float(np.linalg.norm(chiral_vector))
    unit_cell_length = float(np.linalg.norm(translational_vector))
    radius = circumference / (2.0 * math.pi)
    unit_cell_atom_count = 4 * (n * n + n * m + m * m) // translational_gcd

    return NanotubeParameters(
        n=n,
        m=m,
        bond_length=bond_length,
        chiral_vector=chiral_vector,
        translational_vector=translational_vector,
        circumference=circumference,
        radius=radius,
        unit_cell_length=unit_cell_length,
        translational_gcd=translational_gcd,
        t1=t1,
        t2=t2,
        unit_cell_atom_count=unit_cell_atom_count,
    )


def generate_unit_cell_sheet_coordinates(
    params: NanotubeParameters,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    a1, a2, basis = graphene_lattice(params.bond_length)
    cell_matrix = np.column_stack((params.chiral_vector, params.translational_vector))
    inverse_cell_matrix = np.linalg.inv(cell_matrix)
    corners = np.array(
        [
            [0.0, 0.0],
            [float(params.n), float(params.m)],
            [float(params.t1), float(params.t2)],
            [float(params.n + params.t1), float(params.m + params.t2)],
        ],
        dtype=float,
    )

    min_i = math.floor(float(corners[:, 0].min()) - 1.0)
    max_i = math.ceil(float(corners[:, 0].max()) + 1.0)
    min_j = math.floor(float(corners[:, 1].min()) - 1.0)
    max_j = math.ceil(float(corners[:, 1].max()) + 1.0)

    candidates: list[tuple[float, float, int, NDArray[np.float64], NDArray[np.float64]]] = []
    seen: set[tuple[float, float]] = set()

    for i in range(min_i, max_i + 1):
        for j in range(min_j, max_j + 1):
            lattice_origin = i * a1 + j * a2
            for basis_index, basis_offset in enumerate(basis):
                position = lattice_origin + basis_offset
                fractional = _normalize_fractional_coordinates(
                    inverse_cell_matrix @ position
                )
                if _is_inside_half_open_cell(fractional):
                    key = tuple(np.round(position, 10))
                    if key in seen:
                        continue
                    seen.add(key)
                    candidates.append(
                        (
                            float(fractional[1]),
                            float(fractional[0]),
                            basis_index,
                            position.copy(),
                            fractional.copy(),
                        )
                    )

    candidates.sort(key=lambda item: (item[0], item[1], item[2]))
    sheet_coordinates = np.array([item[3] for item in candidates], dtype=float)
    fractional_coordinates = np.array([item[4] for item in candidates], dtype=float)

    if len(sheet_coordinates) != params.unit_cell_atom_count:
        raise RuntimeError(
            "Failed to enumerate the expected number of unit-cell atoms: "
            f"expected {params.unit_cell_atom_count}, got {len(sheet_coordinates)}"
        )

    return sheet_coordinates, fractional_coordinates


def map_sheet_to_cylinder(
    fractional_coordinates: NDArray[np.float64],
    radius: float,
    unit_cell_length: float,
) -> NDArray[np.float64]:
    circumferential_fraction = fractional_coordinates[:, 0]
    axial_fraction = fractional_coordinates[:, 1]

    theta = 2.0 * math.pi * circumferential_fraction
    z = axial_fraction * unit_cell_length
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    return np.column_stack((x, y, z))


def generate_nanotube(
    n: int,
    m: int,
    units: int,
    bond_length: float = 1.42,
    center_z: bool = True,
) -> NanotubeGeometry:
    _validate_positive_int(units, "units")
    params = calculate_nanotube_parameters(n=n, m=m, bond_length=bond_length)
    _, fractional_coordinates = generate_unit_cell_sheet_coordinates(params)
    unit_cell_coordinates = map_sheet_to_cylinder(
        fractional_coordinates=fractional_coordinates,
        radius=params.radius,
        unit_cell_length=params.unit_cell_length,
    )

    coordinates = np.tile(unit_cell_coordinates, (units, 1))
    z_offsets = np.repeat(
        np.arange(units, dtype=float) * params.unit_cell_length,
        unit_cell_coordinates.shape[0],
    )
    coordinates[:, 2] += z_offsets

    if center_z:
        z_center = 0.5 * (float(coordinates[:, 2].min()) + float(coordinates[:, 2].max()))
        coordinates[:, 2] -= z_center

    elements = np.full(coordinates.shape[0], "C", dtype="U2")
    total_length = units * params.unit_cell_length

    return NanotubeGeometry(
        coordinates=coordinates,
        elements=elements,
        radius=params.radius,
        unit_cell_length=params.unit_cell_length,
        total_length=total_length,
        unit_cell_atom_count=params.unit_cell_atom_count,
        atom_count=coordinates.shape[0],
        units=units,
        n=n,
        m=m,
        bond_length=bond_length,
    )


def _normalize_fractional_coordinates(
    fractional: NDArray[np.float64],
) -> NDArray[np.float64]:
    normalized = fractional.copy()
    normalized[np.abs(normalized) < EPSILON] = 0.0
    normalized[np.abs(normalized - 1.0) < EPSILON] = 1.0
    return normalized


def _is_inside_half_open_cell(fractional: NDArray[np.float64]) -> bool:
    return bool(
        np.all(fractional >= -EPSILON) and np.all(fractional < 1.0 - EPSILON)
    )


def _validate_chiral_indices(n: int, m: int) -> None:
    _validate_nonnegative_int(n, "n")
    _validate_nonnegative_int(m, "m")
    if n == 0 and m == 0:
        raise ValueError("(n, m) cannot both be zero")


def _validate_nonnegative_int(value: int, name: str) -> None:
    if value < 0:
        raise ValueError(f"{name} must be non-negative")


def _validate_positive_int(value: int, name: str) -> None:
    if value < 1:
        raise ValueError(f"{name} must be at least 1")
