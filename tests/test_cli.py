from pathlib import Path
import subprocess
import sys

import pytest

from tuber.cli import main


def test_cli_help_returns_zero(capsys) -> None:
    exit_code = main(["--help"])

    captured = capsys.readouterr()
    assert exit_code == 0
    assert "usage:" in captured.out
    assert "--n" in captured.out
    assert "--format" in captured.out


def test_cli_accepts_legacy_generate_alias(capsys) -> None:
    exit_code = main(["generate", "--help"])

    captured = capsys.readouterr()
    assert exit_code == 0
    assert "--n" in captured.out


def test_cli_rejects_invalid_format_choice(capsys, tmp_path) -> None:
    output_path = tmp_path / "invalid.xyz"

    exit_code = main(
        [
            "--n",
            "3",
            "--m",
            "3",
            "--units",
            "1",
            "--format",
            "xyz",
            "--output",
            str(output_path),
        ]
    )

    captured = capsys.readouterr()
    assert exit_code == 2
    assert "invalid choice" in captured.err


def test_cli_rejects_invalid_chiral_indices(capsys, tmp_path) -> None:
    output_path = tmp_path / "invalid.pdb"
    exit_code = main(
        [
            "--n",
            "0",
            "--m",
            "0",
            "--units",
            "1",
            "--format",
            "pdb",
            "--output",
            str(output_path),
        ]
    )

    captured = capsys.readouterr()
    assert exit_code == 2
    assert "cannot both be zero" in captured.err


def test_cli_generate_writes_output_file(tmp_path, capsys) -> None:
    pytest.importorskip("biotite")
    output_path = tmp_path / "tube.pdb"

    exit_code = main(
        [
            "--n",
            "3",
            "--m",
            "3",
            "--units",
            "2",
            "--format",
            "pdb",
            "--output",
            str(output_path),
        ]
    )

    captured = capsys.readouterr()
    assert exit_code == 0
    assert output_path.exists()
    assert "Wrote 24 atoms" in captured.out
    assert "C=24; H=0" in captured.out
    assert str(output_path) in captured.out


def test_cli_requires_overwrite_for_existing_file(tmp_path, capsys) -> None:
    pytest.importorskip("biotite")
    output_path = tmp_path / "tube.pdb"
    output_path.write_text("existing file", encoding="utf-8")

    exit_code = main(
        [
            "--n",
            "3",
            "--m",
            "3",
            "--units",
            "1",
            "--format",
            "pdb",
            "--output",
            str(output_path),
        ]
    )

    captured = capsys.readouterr()
    assert exit_code == 2
    assert "already exists" in captured.err


def test_cli_generate_can_add_terminal_hydrogens(tmp_path, capsys) -> None:
    pytest.importorskip("biotite")
    output_path = tmp_path / "tube.cif"

    exit_code = main(
        [
            "--n",
            "5",
            "--m",
            "0",
            "--units",
            "1",
            "--format",
            "cif",
            "--output",
            str(output_path),
            "--hydrogen-terminate",
        ]
    )

    captured = capsys.readouterr()
    assert exit_code == 0
    assert output_path.exists()
    assert "Wrote 40 atoms" in captured.out
    assert "C=20; H=20" in captured.out


def test_cli_rejects_invalid_hydrogen_bond_length(tmp_path, capsys) -> None:
    output_path = tmp_path / "tube.pdb"

    exit_code = main(
        [
            "--n",
            "3",
            "--m",
            "3",
            "--units",
            "1",
            "--format",
            "pdb",
            "--output",
            str(output_path),
            "--hydrogen-terminate",
            "--hydrogen-bond-length",
            "0",
        ]
    )

    captured = capsys.readouterr()
    assert exit_code == 2
    assert "hydrogen_bond_length must be positive" in captured.err


def test_cli_script_runs_from_repo_root() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [sys.executable, "src/tuber/cli.py", "--help"],
        capture_output=True,
        cwd=repo_root,
        text=True,
        check=False,
    )

    assert result.returncode == 0
    assert "--n" in result.stdout
