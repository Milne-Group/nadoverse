import importlib
import sys
import types

import pytest


class _FakeFastMCP:
    def __init__(self, *args, **kwargs):
        pass

    def tool(self, *args, **kwargs):
        def decorator(func):
            return func

        return decorator

    def run(self, *args, **kwargs):
        pass


def _import_mcp_server(monkeypatch):
    monkeypatch.setitem(
        sys.modules,
        "fastmcp",
        types.SimpleNamespace(FastMCP=_FakeFastMCP),
    )
    monkeypatch.setitem(sys.modules, "seqnado", types.SimpleNamespace(Assay=object))
    sys.modules.pop("nadoverse.mcp_server", None)
    return importlib.import_module("nadoverse.mcp_server")


def test_seqnado_generate_config_uses_non_interactive_cli(monkeypatch, tmp_path):
    mcp_server = _import_mcp_server(monkeypatch)
    output_yaml = tmp_path / "config_atac.yaml"
    calls = []

    def fake_run(args, capture_output, text):
        calls.append(args)
        output_yaml.write_text("genome:\n  name: hg38\n")
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(mcp_server.subprocess, "run", fake_run)

    result = mcp_server.seqnado_generate_config("atac", str(output_yaml))

    assert calls == [
        [
            "seqnado",
            "config",
            "atac",
            "--no-interactive",
            "--no-make-dirs",
            "--render-options",
            "-o",
            str(output_yaml),
        ]
    ]
    assert output_yaml.exists()
    assert "non-interactive SeqNado defaults" in result


def test_seqnado_generate_config_rejects_prompt_only_modes(monkeypatch, tmp_path):
    mcp_server = _import_mcp_server(monkeypatch)

    with pytest.raises(ValueError, match="does not support 'mcc'"):
        mcp_server.seqnado_generate_config("mcc", str(tmp_path / "config_mcc.yaml"))


def test_bamnado_bam_coverage_invokes_cli(monkeypatch, tmp_path):
    mcp_server = _import_mcp_server(monkeypatch)
    monkeypatch.setitem(sys.modules, "bamnado", types.SimpleNamespace())
    calls = []

    def fake_run(args):
        calls.append(args)
        return ""

    monkeypatch.setattr(mcp_server, "_run", fake_run)
    output = tmp_path / "coverage.bedgraph"

    result = mcp_server.bamnado_bam_coverage(
        "/data/sample.bam",
        str(output),
        bin_size=50,
        normalize="raw",
        fragment_counts=True,
        min_mapq=30,
    )

    assert calls == [
        [
            "bamnado",
            "bam-coverage",
            "--bam",
            "/data/sample.bam",
            "--output",
            str(output),
            "--normalize",
            "raw",
            "--threads",
            "6",
            "--strand",
            "both",
            "--min-mapq",
            "30",
            "--min-length",
            "20",
            "--max-length",
            "1000",
            "--bin-size",
            "50",
            "--fragment-counts",
        ]
    ]
    assert result == f"Coverage written to {output}"
