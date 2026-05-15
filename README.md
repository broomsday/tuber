# Tuber

`tuber` generates finite, open-ended carbon nanotube structures from chiral
indices `(n, m)` and a tube length expressed in nanotube translational unit
cells.

The generated coordinates are aligned so the nanotube axis is parallel to the
global `Z` axis.

By default the output is carbon-only. Optional terminal hydrogen termination is
also supported.

## Installation

```bash
python -m pip install -e .
```

For development and tests:

```bash
python -m pip install -e .[dev]
```

## Usage

Generate a nanotube and write it as PDB:

```bash
tuber generate --n 10 --m 5 --units 8 --format pdb --output cnt_10_5_8.pdb
```

Generate a nanotube and write it as CIF:

```bash
tuber generate --n 6 --m 4 --units 3 --format cif --output cnt_6_4_3.cif
```

Generate a hydrogen-terminated nanotube:

```bash
tuber generate --n 5 --m 0 --units 1 --format cif --output cnt_5_0_1_h.cif --hydrogen-terminate
```

## Notes

- `units` means nanotube translational unit cells along the tube axis, not
  rings or arbitrary axial repeats.
- The output is finite and non-periodic.
- The default output is open-ended and carbon-only; `--hydrogen-terminate`
  adds terminal hydrogens to saturate dangling graphene bonds at the tube ends.
