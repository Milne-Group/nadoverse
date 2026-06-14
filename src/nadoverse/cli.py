from __future__ import annotations

import subprocess
import sys
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

from nadoverse.registry import REGISTRY, NadoTool

app = typer.Typer(
    name="nado",
    help="Unified CLI for the *Nado bioinformatics toolkit.",
    no_args_is_help=True,
    add_completion=False,
)

console = Console()


def _require_tool(tool: NadoTool, extra_hint: Optional[str] = None) -> None:
    """Exit with an install hint if the tool's package is not importable."""
    if not tool.is_available():
        extra = extra_hint or tool.install_extra
        console.print(
            f"[bold red]{tool.name}[/bold red] is not installed.\n"
            f"Run: [bold]pip install nadoverse\\[{extra}][/bold]",
            highlight=False,
        )
        raise typer.Exit(1)


def _dispatch(tool: NadoTool, args: list[str], subcommand: Optional[str] = None) -> None:
    """Check availability then exec the tool's CLI entrypoint."""
    _require_tool(tool)
    cmd = [subcommand or tool.cli_command, *args]
    result = subprocess.run(cmd)
    raise typer.Exit(result.returncode)


# ---------------------------------------------------------------------------
# Per-tool subcommands
# ---------------------------------------------------------------------------

@app.command(
    "seq",
    help="Run SeqNado genomics pipelines.",
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
)
def seq(ctx: typer.Context) -> None:
    from nadoverse.registry import get_tool
    _dispatch(get_tool("seqnado"), ctx.args)


@app.command(
    "plot",
    help="Run PlotNado genomic track visualisation.",
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
)
def plot(ctx: typer.Context) -> None:
    from nadoverse.registry import get_tool
    _dispatch(get_tool("plotnado"), ctx.args)


@app.command(
    "track",
    help="Run TrackNado UCSC trackhub generator.",
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
)
def track(ctx: typer.Context) -> None:
    from nadoverse.registry import get_tool
    _dispatch(get_tool("tracknado"), ctx.args)


@app.command(
    "quant",
    help="Run QuantNado genomic signal quantification.",
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
)
def quant(ctx: typer.Context) -> None:
    from nadoverse.registry import get_tool
    _dispatch(get_tool("quantnado"), ctx.args)


@app.command(
    "tab",
    help="Run TabNado TF-binding prediction (requires Python <3.13).",
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
)
def tab(ctx: typer.Context) -> None:
    from nadoverse.registry import get_tool
    tool = get_tool("tabnado")
    if not tool.python_compatible():
        console.print(
            f"[bold yellow]Warning:[/bold yellow] TabNado requires Python <3.13 "
            f"(running {sys.version_info.major}.{sys.version_info.minor}).\n"
            "See https://github.com/Milne-Group/TabNado for updates.",
        )
        raise typer.Exit(1)
    _dispatch(tool, ctx.args)


# ---------------------------------------------------------------------------
# nado doctor
# ---------------------------------------------------------------------------

@app.command("doctor")
def doctor() -> None:
    """Show installed *Nado tools, versions, and compatibility status."""
    table = Table(title="nadoverse — installed tools", show_lines=False)
    table.add_column("Tool", style="bold")
    table.add_column("Extra")
    table.add_column("Installed")
    table.add_column("Version")
    table.add_column("Status")

    for tool in REGISTRY:
        available = tool.is_available()
        version = tool.installed_version() or "—"
        compat = tool.python_compatible()

        if not available:
            status = f"pip install nadoverse[{tool.install_extra}]"
            marker = "[red]✗[/red]"
            version = "—"
        elif not compat:
            status = (
                f"Requires Python <3.13 "
                f"(running {sys.version_info.major}.{sys.version_info.minor})"
            )
            marker = "[yellow]⚠[/yellow]"
        elif tool.cli_command is None:
            status = "Library-only (no CLI)"
            marker = "[green]✓[/green]"
        else:
            status = "OK"
            marker = "[green]✓[/green]"

        table.add_row(tool.name, tool.install_extra, marker, version, status)

    console.print(table)
