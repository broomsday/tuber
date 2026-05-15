from __future__ import annotations

from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx

SUPPORTED_SUFFIXES = {
    "pdb": {".pdb"},
    "cif": {".cif"},
}


def write_structure(
    atom_array: struc.AtomArray,
    output_path: str | Path,
    file_format: str,
    overwrite: bool = False,
) -> Path:
    normalized_format = file_format.lower()
    if normalized_format not in SUPPORTED_SUFFIXES:
        raise ValueError(
            f"Unsupported format {file_format!r}; expected one of "
            f"{', '.join(sorted(SUPPORTED_SUFFIXES))}"
        )

    output_path = Path(output_path)
    _validate_output_path(
        output_path=output_path,
        file_format=normalized_format,
        overwrite=overwrite,
    )

    if normalized_format == "pdb":
        pdb_file = pdb.PDBFile()
        pdb_file.set_structure(atom_array)
        pdb_file.write(str(output_path))
    else:
        cif_file = pdbx.CIFFile()
        pdbx.set_structure(cif_file, atom_array, data_block="TUBER")
        cif_file.write(str(output_path))

    return output_path


def _validate_output_path(
    output_path: Path,
    file_format: str,
    overwrite: bool,
) -> None:
    if output_path.suffix.lower() not in SUPPORTED_SUFFIXES[file_format]:
        expected_suffixes = ", ".join(sorted(SUPPORTED_SUFFIXES[file_format]))
        raise ValueError(
            f"Output path {output_path} does not match requested format "
            f"{file_format!r}; expected suffix {expected_suffixes}"
        )

    if output_path.exists() and output_path.is_dir():
        raise IsADirectoryError(f"Output path {output_path} is a directory")

    if output_path.exists() and not overwrite:
        raise FileExistsError(
            f"Output path {output_path} already exists; use --overwrite to replace it"
        )

    if not output_path.parent.exists():
        raise FileNotFoundError(
            f"Parent directory does not exist for output path {output_path}"
        )
