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
    carbon_count: int
    hydrogen_count: int
    hydrogen_terminated: bool


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
    sheet_coordinates, fractional_coordinates, _ = _generate_unit_cell_sheet_atom_data(params)
    return sheet_coordinates, fractional_coordinates


def _generate_unit_cell_sheet_atom_data(
    params: NanotubeParameters,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.int_]]:
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
    basis_indices = np.array([item[2] for item in candidates], dtype=int)

    if len(sheet_coordinates) != params.unit_cell_atom_count:
        raise RuntimeError(
            "Failed to enumerate the expected number of unit-cell atoms: "
            f"expected {params.unit_cell_atom_count}, got {len(sheet_coordinates)}"
        )

    return sheet_coordinates, fractional_coordinates, basis_indices


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
    hydrogen_terminate: bool = False,
    hydrogen_bond_length: float = 1.09,
) -> NanotubeGeometry:
    _validate_positive_int(units, "units")
    params = calculate_nanotube_parameters(n=n, m=m, bond_length=bond_length)
    unit_cell_sheet_coordinates, unit_cell_fractional_coordinates, unit_cell_basis_indices = (
        _generate_unit_cell_sheet_atom_data(params)
    )

    unit_cell_atom_count = unit_cell_sheet_coordinates.shape[0]
    unit_offsets = np.repeat(np.arange(units, dtype=float), unit_cell_atom_count)
    sheet_offsets = unit_offsets[:, np.newaxis] * params.translational_vector
    carbon_sheet_coordinates = np.tile(unit_cell_sheet_coordinates, (units, 1)) + sheet_offsets

    carbon_fractional_coordinates = np.tile(unit_cell_fractional_coordinates, (units, 1))
    carbon_fractional_coordinates[:, 1] += unit_offsets
    carbon_basis_indices = np.tile(unit_cell_basis_indices, units)

    carbon_coordinates = map_sheet_to_cylinder(
        fractional_coordinates=carbon_fractional_coordinates,
        radius=params.radius,
        unit_cell_length=params.unit_cell_length,
    )

    coordinates = carbon_coordinates.copy()
    elements = np.full(carbon_coordinates.shape[0], "C", dtype="U2")
    hydrogen_count = 0

    if hydrogen_terminate:
        _validate_positive_float(hydrogen_bond_length, "hydrogen_bond_length")
        hydrogen_coordinates = _generate_terminal_hydrogen_coordinates(
            carbon_coordinates=carbon_coordinates,
            carbon_sheet_coordinates=carbon_sheet_coordinates,
            carbon_fractional_coordinates=carbon_fractional_coordinates,
            carbon_basis_indices=carbon_basis_indices,
            params=params,
            hydrogen_bond_length=hydrogen_bond_length,
        )
        if hydrogen_coordinates.size > 0:
            coordinates = np.vstack((coordinates, hydrogen_coordinates))
            elements = np.concatenate(
                (
                    elements,
                    np.full(hydrogen_coordinates.shape[0], "H", dtype="U2"),
                )
            )
            hydrogen_count = int(hydrogen_coordinates.shape[0])

    if center_z:
        z_center = 0.5 * (float(coordinates[:, 2].min()) + float(coordinates[:, 2].max()))
        coordinates[:, 2] -= z_center

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
        carbon_count=int(carbon_coordinates.shape[0]),
        hydrogen_count=hydrogen_count,
        hydrogen_terminated=hydrogen_terminate,
    )


def _generate_terminal_hydrogen_coordinates(
    carbon_coordinates: NDArray[np.float64],
    carbon_sheet_coordinates: NDArray[np.float64],
    carbon_fractional_coordinates: NDArray[np.float64],
    carbon_basis_indices: NDArray[np.int_],
    params: NanotubeParameters,
    hydrogen_bond_length: float,
) -> NDArray[np.float64]:
    a1, a2, basis = graphene_lattice(params.bond_length)
    inverse_cell_matrix = np.linalg.inv(
        np.column_stack((params.chiral_vector, params.translational_vector))
    )
    present_positions = {
        _fractional_position_key(fractional_position)
        for fractional_position in carbon_fractional_coordinates
    }

    hydrogen_coordinates: list[NDArray[np.float64]] = []
    for atom_coordinate, sheet_position, basis_index, current_fractional in zip(
        carbon_coordinates,
        carbon_sheet_coordinates,
        carbon_basis_indices,
        carbon_fractional_coordinates,
        strict=True,
    ):
        theta = 2.0 * math.pi * float(current_fractional[0])
        circumferential_direction = np.array(
            [-math.sin(theta), math.cos(theta), 0.0],
            dtype=float,
        )
        axial_direction = np.array([0.0, 0.0, 1.0], dtype=float)
        neighbor_positions = _graphene_neighbor_positions(
            position=sheet_position,
            basis_index=int(basis_index),
            a1=a1,
            a2=a2,
            basis=basis,
        )

        missing_neighbor_positions = [
            neighbor_position
            for neighbor_position in neighbor_positions
            if _fractional_position_key(inverse_cell_matrix @ neighbor_position)
            not in present_positions
        ]

        for missing_neighbor_position in missing_neighbor_positions:
            sheet_displacement = missing_neighbor_position - sheet_position
            fractional_displacement = inverse_cell_matrix @ sheet_displacement
            hydrogen_direction = (
                fractional_displacement[0] * params.circumference * circumferential_direction
                + fractional_displacement[1] * params.unit_cell_length * axial_direction
            )
            direction_norm = float(np.linalg.norm(hydrogen_direction))
            if direction_norm <= EPSILON:
                raise RuntimeError(
                    "Failed to determine a stable hydrogen direction for a terminal carbon"
                )

            hydrogen_coordinates.append(
                atom_coordinate + hydrogen_bond_length * (hydrogen_direction / direction_norm)
            )

    if not hydrogen_coordinates:
        return np.empty((0, 3), dtype=float)

    return np.asarray(hydrogen_coordinates, dtype=float)


def _graphene_neighbor_positions(
    position: NDArray[np.float64],
    basis_index: int,
    a1: NDArray[np.float64],
    a2: NDArray[np.float64],
    basis: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    if basis_index == 0:
        basis_offset = basis[1] - basis[0]
        return (
            position + basis_offset,
            position + basis_offset - a1,
            position + basis_offset - a2,
        )

    if basis_index == 1:
        basis_offset = basis[0] - basis[1]
        return (
            position + basis_offset,
            position + basis_offset + a1,
            position + basis_offset + a2,
        )

    raise ValueError(f"Unsupported graphene basis index {basis_index}")


def _fractional_position_key(position: NDArray[np.float64]) -> tuple[float, float]:
    normalized = position.copy()
    normalized[0] = normalized[0] % 1.0
    if abs(normalized[0] - 1.0) < EPSILON:
        normalized[0] = 0.0
    if abs(normalized[0]) < EPSILON:
        normalized[0] = 0.0
    if abs(normalized[1]) < EPSILON:
        normalized[1] = 0.0
    return tuple(np.round(normalized, 10))


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


def _validate_positive_float(value: float, name: str) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be positive")
