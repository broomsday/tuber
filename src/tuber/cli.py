from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path
import sys

if __name__ == "__main__" and __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "tuber"

from .geometry import generate_nanotube


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tuber",
        description=(
            "Generate finite carbon nanotube structures aligned to the global Z axis."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--n", type=int, required=True, help="First chiral index")
    parser.add_argument("--m", type=int, required=True, help="Second chiral index")
    parser.add_argument(
        "--units",
        type=int,
        required=True,
        help="Number of nanotube translational unit cells along Z",
    )
    parser.add_argument(
        "--format",
        dest="file_format",
        choices=("pdb", "cif"),
        required=True,
        help="Output format",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Destination file path",
    )
    parser.add_argument(
        "--bond-length",
        type=float,
        default=1.42,
        help="Carbon-carbon bond length in angstroms",
    )
    parser.add_argument(
        "--center-z",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Center the final tube around z = 0",
    )
    parser.add_argument(
        "--hydrogen-terminate",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Add one terminal hydrogen to each undercoordinated end carbon",
    )
    parser.add_argument(
        "--hydrogen-bond-length",
        type=float,
        default=1.09,
        help="Carbon-hydrogen bond length in angstroms when hydrogen termination is enabled",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow overwriting an existing output file",
    )

    return parser


def _normalize_argv(argv: Sequence[str] | None) -> list[str]:
    normalized_argv = list(sys.argv[1:] if argv is None else argv)
    if normalized_argv and normalized_argv[0] == "generate":
        return normalized_argv[1:]
    return normalized_argv


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    normalized_argv = _normalize_argv(argv)
    try:
        args = parser.parse_args(normalized_argv)
    except SystemExit as error:
        return int(error.code)
    return _run_generate(args)


def _run_generate(args: argparse.Namespace) -> int:
    try:
        geometry = generate_nanotube(
            n=args.n,
            m=args.m,
            units=args.units,
            bond_length=args.bond_length,
            center_z=args.center_z,
            hydrogen_terminate=args.hydrogen_terminate,
            hydrogen_bond_length=args.hydrogen_bond_length,
        )

        from .structure import build_atom_array
        from .writers import write_structure

        atom_array = build_atom_array(
            coordinates=geometry.coordinates,
            elements=geometry.elements,
        )
        output_path = write_structure(
            atom_array=atom_array,
            output_path=args.output,
            file_format=args.file_format,
            overwrite=args.overwrite,
        )
    except (OSError, RuntimeError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2

    print(
        "Wrote "
        f"{geometry.atom_count} atoms; "
        f"C={geometry.carbon_count}; "
        f"H={geometry.hydrogen_count}; "
        f"radius={geometry.radius:.3f} A; "
        f"total_length={geometry.total_length:.3f} A; "
        f"output={output_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
