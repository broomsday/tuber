# Carbon Nanotube CLI Plan

## Goal

Turn this repository into a Python CLI application that generates finite carbon nanotube structures from chiral indices `(n, m)` and a tube length expressed in nanotube translational unit cells, with the tube axis aligned to the global `Z` axis.

Inputs:

- `n`: first chiral index
- `m`: second chiral index
- `units`: number of nanotube translational unit cells along `Z`
- `format`: output format, `pdb` or `cif`
- `output`: destination file path

Outputs:

- finite, open-ended, carbon-only nanotube structure
- coordinates aligned so the tube axis is parallel to `Z`
- serialization as either PDB or CIF

Example:

```bash
tuber generate --n 10 --m 5 --units 8 --format pdb --output cnt_10_5_8.pdb
```

## Locked Decisions

The following design choices are now fixed:

- `units` means nanotube translational unit cells, not rings or arbitrary axial repeats
- geometry should be built as NumPy coordinate arrays
- structure export should use Biotite `AtomArray`
- the structure is carbon-only
- the export is finite and non-periodic; the output should behave like a normal molecular structure file rather than a periodic crystal model

## High-Level Architecture

Use a small `src/`-layout Python package:

```text
pyproject.toml
src/
  tuber/
    __init__.py
    cli.py
    geometry.py
    structure.py
    writers.py
tests/
README.md
PLAN.md
```

Module responsibilities:

- `cli.py`: parse command-line arguments and drive generation/export
- `geometry.py`: graphene and nanotube math, coordinate generation, validation
- `structure.py`: build Biotite `AtomArray` objects from generated coordinates
- `writers.py`: dispatch Biotite-based save logic for PDB and CIF

## Core Geometry Strategy

Implement the nanotube from the standard graphene-to-cylinder construction.

### 1. Graphene lattice model

Use a standard carbon-carbon bond length, likely `1.42 A`, and define the graphene primitive lattice vectors and two-atom basis.

### 2. Chiral and translational vectors

For a requested `(n, m)`:

- compute the chiral vector `Ch = n a1 + m a2`
- compute the translational vector `T`
- derive the tube circumference, radius, and translational unit-cell length
- derive the number of atoms in one nanotube translational unit cell

This gives the exact finite repeat to replicate along `Z`.

### 3. 2D fundamental-domain atom generation

Generate all graphene basis atoms that fall inside one translational nanotube unit cell in the unwrapped 2D sheet.

Requirements:

- robust inclusion test for boundary points
- deterministic handling of edge duplicates
- stable atom ordering so output is reproducible

### 4. Cylindrical mapping

Map the 2D unwrapped coordinates to 3D Cartesian coordinates:

- circumferential position -> angular coordinate `theta`
- axial position -> `z`
- `x = R cos(theta)`
- `y = R sin(theta)`

Construct the coordinates so the nanotube axis is aligned to `Z` by construction rather than by a later rotation step.

### 5. Axial replication

Replicate the translational unit cell `units` times along `Z`.

Default placement convention:

- center the final tube around `z = 0` for symmetry

## Data Model and Export

### NumPy-first representation

The geometry layer should return plain NumPy arrays first:

- coordinates: shape `(N, 3)`
- element symbols: shape `(N,)`, all `"C"`
- optional metadata: radius, unit-cell length, total length, atom count

This keeps the math layer independent from file formats.

### Biotite conversion

Build a Biotite `AtomArray` from the NumPy arrays.

Populate at least:

- Cartesian coordinates
- element name
- atom name, e.g. `C`
- residue or molecule identifiers sufficient for valid PDB/CIF export

### File writing

Use Biotite's structure I/O to save:

- PDB
- CIF

Since this is a finite, non-periodic molecule export, avoid introducing periodic-crystal semantics into the output model.

## CLI Design

Proposed command:

```bash
tuber generate --n <int> --m <int> --units <int> --format <pdb|cif> --output <path>
```

Recommended options:

- `--bond-length`: optional override, default `1.42`
- `--center-z / --no-center-z`: optional control over axial placement
- `--overwrite`: optional explicit overwrite behavior

Validation rules:

- reject `(n, m) = (0, 0)`
- require `n >= 0`, `m >= 0`
- require `units >= 1`
- require supported output format
- fail clearly on invalid output suffix or unwritable path

## Implementation Phases

### Phase 1: Project scaffolding

- add `pyproject.toml`
- define package metadata and CLI entry point
- add dependencies: `numpy`, `biotite`
- create package skeleton under `src/tuber`

### Phase 2: Geometry core

- implement graphene lattice constants and basis
- implement `(n, m)` nanotube parameter calculations
- implement unit-cell atom enumeration in the unwrapped sheet
- implement cylindrical mapping to Cartesian coordinates
- implement axial replication and centering

### Phase 3: Structure conversion and export

- build Biotite `AtomArray` objects from NumPy arrays
- implement PDB writing
- implement CIF writing
- ensure both outputs are valid for finite molecular structures

### Phase 4: CLI integration

- parse command-line arguments
- call geometry generation
- convert to `AtomArray`
- write the requested file
- print a short summary: atom count, radius, total length, output path

### Phase 5: Tests and verification

- unit tests for nanotube geometry formulas
- regression tests for atom counts
- export smoke tests for both formats
- CLI tests for argument validation and file generation

## Test Strategy

Test small representative tubes such as:

- armchair: `(3, 3)`
- zigzag: `(5, 0)`
- chiral: `(6, 4)`

Core assertions:

- atom count matches the analytical unit-cell formula times `units`
- all atoms lie on the same cylinder radius within tolerance
- total tube length matches `units * |T|`
- output coordinates are aligned to `Z`
- PDB and CIF exports can be written successfully by Biotite

## Documentation Requirements

Update `README.md` to document:

- what the tool does
- installation steps
- CLI usage examples
- exact meaning of `units`
- coordinate convention: tube axis along `Z`
- current scope: finite, open-ended, carbon-only nanotubes

## Suggested Initial Deliverable

The first implementation milestone should produce:

- working `tuber generate` command
- correct geometry for finite carbon nanotubes
- export to both PDB and CIF through Biotite
- basic tests covering one armchair, one zigzag, and one chiral example

Status on 2026-05-14: complete.

Delivered:

- working `tuber generate` command
- finite nanotube geometry aligned to `Z`
- PDB and CIF export through Biotite
- tests covering armchair `(3, 3)`, zigzag `(5, 0)`, and chiral `(6, 4)` cases

## Future Extensions

These are explicitly out of scope for the first version:

- hydrogen termination
- capped nanotubes
- periodic crystal export
- multi-wall nanotubes
- defect insertion or doping
- alternate coordinate conventions beyond simple axial centering
