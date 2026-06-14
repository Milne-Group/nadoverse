"""
Local stdio MCP server for the *Nado bioinformatics toolkit.

Exposes each *Nado tool with a CLI entrypoint as an MCP tool.
Library-only tools (BamNado, MCCNado) appear in list_nado_tools() only.

Install the mcp extra then add to Claude Code:
    pip install nadoverse[mcp]

    # ~/.claude/mcp.json
    { "mcpServers": { "nadoverse": { "command": "nadoverse-mcp" } } }
"""

from __future__ import annotations

import json
import subprocess
import sys
from typing import Optional

from nadoverse.registry import REGISTRY, NadoTool, all_tools

try:
    import fastmcp
    from fastmcp import FastMCP
except ImportError:
    raise ImportError(
        "MCP support requires fastmcp. Run: pip install nadoverse[mcp]"
    )

mcp = FastMCP(
    name="nadoverse",
    instructions=(
        "Tools for running *Nado bioinformatics pipelines. "
        "Always call list_nado_tools first to check what is installed. "
        "Use dry_run=True (default) to preview commands before execution."
    ),
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _tool_unavailable_msg(tool: NadoTool) -> str:
    return (
        f"{tool.name} is not installed. "
        f"Run: pip install nadoverse[{tool.install_extra}]"
    )


def _run(cmd: list[str], dry_run: bool) -> str:
    if dry_run:
        return f"[dry-run] would execute: {' '.join(cmd)}"
    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stdout + result.stderr
    return output.strip() if output.strip() else f"[exit {result.returncode}]"


def _check_elicitation_support(ctx) -> bool:
    """Return True if the MCP client supports elicitation."""
    try:
        caps = ctx.client_capabilities or {}
        return bool(caps.get("elicitation"))
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Tools
# ---------------------------------------------------------------------------

@mcp.tool(
    description=(
        "List all *Nado tools: installation status, version, CLI entrypoints, "
        "input/output types. Run this first to know what is available."
    ),
    annotations={"readOnlyHint": True, "destructiveHint": False},
)
def list_nado_tools() -> str:
    rows = []
    for tool in all_tools():
        rows.append({
            "name": tool.name,
            "installed": tool.is_available(),
            "version": tool.installed_version(),
            "cli_command": tool.cli_command,
            "install_hint": f"pip install nadoverse[{tool.install_extra}]",
            "input_types": tool.input_types,
            "output_types": tool.output_types,
            "python_compatible": tool.python_compatible(),
            "description": tool.description,
        })
    return json.dumps(rows, indent=2)


@mcp.tool(
    description=(
        "Run a SeqNado genomics pipeline (ATAC-seq, ChIP-seq, CUT&RUN/TAG, "
        "RNA-seq, WGS, Methylation, CRISPR, Micro-Capture-C). "
        "Set dry_run=True to preview the command without executing."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def run_seqnado(
    args: list[str],
    dry_run: bool = True,
) -> str:
    from nadoverse.registry import get_tool
    tool = get_tool("seqnado")
    if not tool.is_available():
        return _tool_unavailable_msg(tool)
    return _run(["seqnado", *args], dry_run)


@mcp.tool(
    description=(
        "Run PlotNado to generate genomic track visualisations "
        "(UCSC-style browser tracks, publication-quality plots from bigWig/bigBed). "
        "Set dry_run=True to preview."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def run_plotnado(
    args: list[str],
    dry_run: bool = True,
) -> str:
    from nadoverse.registry import get_tool
    tool = get_tool("plotnado")
    if not tool.is_available():
        return _tool_unavailable_msg(tool)
    return _run(["plotnado", *args], dry_run)


@mcp.tool(
    description=(
        "Run TrackNado to generate UCSC trackhubs from bigWig/bigBed/BED files. "
        "Set dry_run=True to preview."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def run_tracknado(
    args: list[str],
    dry_run: bool = True,
) -> str:
    from nadoverse.registry import get_tool
    tool = get_tool("tracknado")
    if not tool.is_available():
        return _tool_unavailable_msg(tool)
    return _run(["tracknado", *args], dry_run)


@mcp.tool(
    description=(
        "Run QuantNado for genomic signal quantification: BAM/bigWig ingestion, "
        "Zarr storage, feature counting, dimensionality reduction, peak calling. "
        "Set dry_run=True to preview."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def run_quantnado(
    subcommand: str,
    args: list[str],
    dry_run: bool = True,
) -> str:
    """
    subcommand: one of 'run', 'make-zarr', 'combine-metadata'
    """
    from nadoverse.registry import get_tool
    tool = get_tool("quantnado")
    if not tool.is_available():
        return _tool_unavailable_msg(tool)
    cmd_map = {
        "run": "quantnado",
        "make-zarr": "quantnado-make-zarr",
        "combine-metadata": "quantnado-combine-metadata",
    }
    if subcommand not in cmd_map:
        return f"Unknown subcommand '{subcommand}'. Choose from: {list(cmd_map)}"
    return _run([cmd_map[subcommand], *args], dry_run)


@mcp.tool(
    description=(
        "Run TabNado TF-binding prediction from epigenomic cofactors. "
        "Requires Python <3.13. Subcommand: init, run, data, sweep, train, evaluate, or shap. "
        "Set dry_run=True to preview."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def run_tabnado(
    subcommand: str,
    args: list[str],
    dry_run: bool = True,
) -> str:
    from nadoverse.registry import get_tool
    tool = get_tool("tabnado")
    if not tool.is_available():
        return _tool_unavailable_msg(tool)
    if not tool.python_compatible():
        return (
            f"TabNado requires Python <3.13 "
            f"(running {sys.version_info.major}.{sys.version_info.minor}). "
            "See https://github.com/Milne-Group/TabNado for updates."
        )
    cmd_map = {
        "init": "tabnado-init",
        "run": "tabnado-run",
        "data": "tabnado-data",
        "sweep": "tabnado-sweep",
        "train": "tabnado-train",
        "evaluate": "tabnado-evaluate",
        "shap": "tabnado-shap",
    }
    if subcommand not in cmd_map:
        return f"Unknown subcommand '{subcommand}'. Choose from: {list(cmd_map)}"
    return _run([cmd_map[subcommand], *args], dry_run)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
