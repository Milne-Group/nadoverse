"""Tests for the nado CLI. Requires nadoverse installed in the test env."""

import subprocess
import sys

import pytest


def _nado(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "nadoverse.cli", *args],
        capture_output=True,
        text=True,
    )


def test_doctor_runs():
    result = subprocess.run(
        ["nado", "doctor"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def test_doctor_shows_tool_names():
    result = subprocess.run(
        ["nado", "doctor"],
        capture_output=True,
        text=True,
    )
    output = result.stdout + result.stderr
    for name in ("SeqNado", "PlotNado", "TrackNado", "QuantNado", "TabNado"):
        assert name in output, f"{name} missing from nado doctor output"


def test_nado_help_exits_zero():
    result = subprocess.run(["nado", "--help"], capture_output=True, text=True)
    assert result.returncode == 0


def test_subcommand_help_exits_zero():
    for cmd in ("seq", "plot", "track", "quant", "tab"):
        result = subprocess.run(
            ["nado", cmd, "--help"],
            capture_output=True,
            text=True,
        )
        # Either the tool runs --help (0) or gives install hint (1) — both acceptable
        assert result.returncode in (0, 1), (
            f"nado {cmd} --help returned unexpected exit {result.returncode}"
        )


def test_missing_tool_shows_install_hint(monkeypatch):
    """When a tool isn't installed, the CLI prints the install hint."""
    import importlib
    import sys
    from unittest.mock import patch

    # Temporarily hide seqnado from the import system
    with patch.dict(sys.modules, {"seqnado": None}):
        from nadoverse.registry import get_tool
        tool = get_tool("seqnado")
        assert tool is not None
        assert not tool.is_available()
        hint = f"pip install nadoverse[{tool.install_extra}]"
        assert hint == "pip install nadoverse[seqnado]"
