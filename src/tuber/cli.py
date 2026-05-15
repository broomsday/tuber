from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path
import sys

from .geometry import generate_nanotube


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tuber",
        description=(
            "Generate finite carbon nanotube structures aligned to the global Z axis."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    generate_parser = subparsers.add_parser(
        "generate",
        help="Generate a finite carbon nanotube",
    )
    generate_parser.add_argument("--n", type=int, required=True, help="First chiral index")
    generate_parser.add_argument("--m", type=int, required=True, help="Second chiral index")
    generate_parser.add_argument(
        "--units",
        type=int,
        required=True,
        help="Number of nanotube translational unit cells along Z",
    )
    generate_parser.add_argument(
        "--format",
        dest="file_format",
        choices=("pdb", "cif"),
        required=True,
        help="Output format",
    )
    generate_parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Destination file path",
    )
    generate_parser.add_argument(
        "--bond-length",
        type=float,
        default=1.42,
        help="Carbon-carbon bond length in angstroms",
    )
    generate_parser.add_argument(
        "--center-z",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Center the final tube around z = 0",
    )
    generate_parser.add_argument(
        "--hydrogen-terminate",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Add one terminal hydrogen to each undercoordinated end carbon",
    )
    generate_parser.add_argument(
        "--hydrogen-bond-length",
        type=float,
        default=1.09,
        help="Carbon-hydrogen bond length in angstroms when hydrogen termination is enabled",
    )
    generate_parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow overwriting an existing output file",
    )

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as error:
        return int(error.code)

    if args.command == "generate":
        return _run_generate(args)

    parser.error(f"Unknown command: {args.command}")
    return 2


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
