import math

import numpy as np
import pytest

from tuber.geometry import calculate_nanotube_parameters, generate_nanotube


@pytest.mark.parametrize(
    ("n", "m", "expected_atom_count"),
    [
        (3, 3, 12),
        (5, 0, 20),
        (6, 4, 152),
    ],
)
def test_unit_cell_atom_count_matches_analytical_formula(
    n: int,
    m: int,
    expected_atom_count: int,
) -> None:
    params = calculate_nanotube_parameters(n=n, m=m)
    assert params.unit_cell_atom_count == expected_atom_count


@pytest.mark.parametrize(
    ("n", "m", "units"),
    [
        (3, 3, 2),
        (5, 0, 3),
        (6, 4, 2),
    ],
)
def test_generated_nanotube_has_constant_radius_and_z_axis_replication(
    n: int,
    m: int,
    units: int,
) -> None:
    geometry = generate_nanotube(n=n, m=m, units=units, center_z=False)
    radii = np.linalg.norm(geometry.coordinates[:, :2], axis=1)
    assert np.allclose(radii, geometry.radius, atol=1e-10)

    unit_cell_atom_count = geometry.unit_cell_atom_count
    replicated_delta = (
        geometry.coordinates[unit_cell_atom_count : 2 * unit_cell_atom_count]
        - geometry.coordinates[:unit_cell_atom_count]
    )
    assert np.allclose(replicated_delta[:, :2], 0.0, atol=1e-10)
    assert np.allclose(
        replicated_delta[:, 2],
        geometry.unit_cell_length,
        atol=1e-10,
    )


def test_total_length_metadata_matches_units_times_translational_repeat() -> None:
    geometry = generate_nanotube(n=6, m=4, units=3)
    assert math.isclose(
        geometry.total_length,
        geometry.units * geometry.unit_cell_length,
        rel_tol=1e-12,
    )
    assert geometry.atom_count == geometry.units * geometry.unit_cell_atom_count


def test_center_z_controls_axial_placement() -> None:
    centered = generate_nanotube(n=5, m=0, units=4, center_z=True)
    uncentered = generate_nanotube(n=5, m=0, units=4, center_z=False)

    assert centered.coordinates[:, 2].min() == pytest.approx(
        -centered.coordinates[:, 2].max(),
        abs=1e-10,
    )
    assert uncentered.coordinates[:, 2].min() == pytest.approx(0.0, abs=1e-10)


@pytest.mark.parametrize(
    ("n", "m", "units", "expected_hydrogen_count"),
    [
        (3, 3, 2, 12),
        (5, 0, 1, 20),
        (6, 4, 1, 24),
    ],
)
def test_hydrogen_termination_adds_expected_end_hydrogens(
    n: int,
    m: int,
    units: int,
    expected_hydrogen_count: int,
) -> None:
    geometry = generate_nanotube(
        n=n,
        m=m,
        units=units,
        hydrogen_terminate=True,
        center_z=False,
    )

    assert geometry.carbon_count == geometry.unit_cell_atom_count * units
    assert geometry.hydrogen_count == expected_hydrogen_count
    assert geometry.atom_count == geometry.carbon_count + geometry.hydrogen_count
    assert geometry.hydrogen_terminated is True

    carbon_coordinates = geometry.coordinates[geometry.elements == "C"]
    hydrogen_coordinates = geometry.coordinates[geometry.elements == "H"]
    carbon_hydrogen_distances = np.linalg.norm(
        hydrogen_coordinates[:, np.newaxis, :] - carbon_coordinates[np.newaxis, :, :],
        axis=2,
    )

    assert np.allclose(np.min(carbon_hydrogen_distances, axis=1), 1.09, atol=1e-10)


def test_default_generation_remains_carbon_only() -> None:
    geometry = generate_nanotube(n=3, m=3, units=2)

    assert geometry.hydrogen_count == 0
    assert geometry.carbon_count == geometry.atom_count
    assert geometry.hydrogen_terminated is False
    assert np.array_equal(np.unique(geometry.elements), np.array(["C"]))


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"n": -1, "m": 1, "units": 1}, "n must be non-negative"),
        ({"n": 0, "m": 0, "units": 1}, "cannot both be zero"),
        ({"n": 1, "m": 1, "units": 0}, "units must be at least 1"),
        ({"n": 1, "m": 1, "units": 1, "bond_length": 0.0}, "bond_length must be positive"),
        (
            {"n": 1, "m": 1, "units": 1, "hydrogen_terminate": True, "hydrogen_bond_length": 0.0},
            "hydrogen_bond_length must be positive",
        ),
    ],
)
def test_invalid_inputs_raise_value_error(kwargs: dict[str, int | float], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        generate_nanotube(**kwargs)
