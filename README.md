# Tuber

`tuber` generates idealized, finite carbon nanotube structures from chiral
indices `(n, m)` and a tube length expressed as an integer number of nanotube
translational unit cells.

The generated coordinates are aligned to the global `Z` axis. Output can remain
carbon-only or include terminal hydrogens at the two open ends, and structures
can be written as `PDB` or `CIF`.

## Purpose

Use `tuber` when you want a simple, scriptable way to build open-ended carbon
nanotube geometries for visualization, inspection, or downstream workflows that
can consume `PDB` or `CIF`.

The current codebase is focused on:

- finite, non-periodic nanotubes
- chiral-index-driven geometry generation
- optional hydrogen termination of dangling end bonds
- deterministic coordinate generation aligned to `Z`
- file export through Biotite

## Images

Add renders or screenshots here later.

- `[placeholder]` overview render of a generated nanotube
- `[placeholder]` close-up of an open end with hydrogen termination enabled
- `[placeholder]` comparison showing centered vs. uncentered `Z` placement

## Installation

From the project root:

```bash
uv sync
```

For development and tests:

```bash
uv sync --extra dev
```

## Example Commands

Show the CLI:

```bash
uv run tuber --help
```

Generate a carbon-only nanotube as PDB:

```bash
uv run tuber --n 10 --m 5 --units 8 --format pdb --output cnt_10_5_8.pdb
```

Generate a nanotube as CIF:

```bash
uv run tuber --n 6 --m 4 --units 3 --format cif --output cnt_6_4_3.cif
```

Generate a hydrogen-terminated nanotube:

```bash
uv run tuber --n 5 --m 0 --units 1 --format cif --output cnt_5_0_1_h.cif --hydrogen-terminate
```

Keep the tube uncentered along `Z` so the minimum axial coordinate stays near
`0` instead of centering the structure around `z = 0`:

```bash
uv run tuber --n 5 --m 0 --units 4 --format pdb --output cnt_5_0_4_uncentered.pdb --no-center-z
```

Overwrite an existing file:

```bash
uv run tuber --n 3 --m 3 --units 2 --format pdb --output tube.pdb --overwrite
```

The CLI also accepts a legacy `generate` token before the normal options:

```bash
uv run tuber generate --n 3 --m 3 --units 2 --format pdb --output tube.pdb
```

## CLI Notes

- Required inputs are `--n`, `--m`, `--units`, `--format`, and `--output`.
- Supported formats are `pdb` and `cif`.
- `--center-z` is enabled by default; use `--no-center-z` to disable centering.
- `--hydrogen-terminate` is disabled by default.
- `--hydrogen-bond-length` only matters when hydrogen termination is enabled.
- The output filename suffix must match the selected format.
- Existing files are not replaced unless `--overwrite` is passed.
- After writing a file, the CLI prints a summary including total atom count,
  carbon count, hydrogen count, radius, total length, and the output path.

## Python API

The package also exposes a small programmatic surface:

```python
from tuber import generate_nanotube

geometry = generate_nanotube(
    n=6,
    m=4,
    units=3,
    hydrogen_terminate=True,
)

print(geometry.atom_count)
print(geometry.radius)
print(geometry.total_length)
print(geometry.elements[:5])
```

## Implementation Notes

- `src/tuber/geometry.py` builds the graphene lattice, derives nanotube
  parameters from `(n, m)`, enumerates atoms in a half-open unit-cell sheet
  region to avoid duplicates, maps the sheet onto a cylinder, and optionally
  adds terminal hydrogens by detecting missing graphene neighbors at the open
  ends.
- `src/tuber/structure.py` converts generated coordinates and element labels
  into a Biotite `AtomArray`. The current annotations use chain `A`, residue
  name `CNT`, and one residue for the whole structure.
- `src/tuber/writers.py` writes either `PDB` or `CIF` and validates suffix
  compatibility, parent-directory existence, and overwrite behavior.
- `src/tuber/cli.py` is the main command entry point and is also exposed as the
  `tuber` console script through `pyproject.toml`.

## Assumptions And Limitations

- `units` means nanotube translational unit cells along the tube axis, not
  rings or an arbitrary physical length.
- Generated structures are finite, open-ended, and non-periodic.
- Inputs must satisfy `n >= 0`, `m >= 0`, `(n, m) != (0, 0)`, and `units >= 1`.
- The geometry is idealized from analytical lattice construction with fixed bond
  lengths; there is no force-field relaxation or energy minimization step.
- Hydrogen termination adds one hydrogen per undercoordinated end carbon using a
  fixed `C-H` bond length.
- The nanotube axis is always aligned to global `Z`; only the choice of whether
  the final coordinates are centered around `z = 0` is configurable.

## Development

Run the test suite with:

```bash
uv run --extra dev pytest
```
