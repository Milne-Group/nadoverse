# nadoverse

Unified installer, CLI, and MCP server for the **\*Nado** bioinformatics toolkit.

```bash
uv venv .venv --python 3.12
source .venv/bin/activate
uv pip install "nadoverse[seqnado,plotnado,quantnado]"
nado doctor
```

## Included Tools

| Tool | Extra | Description | CLI |
|------|-------|-------------|-----|
| [SeqNado](https://github.com/Milne-Group/SeqNado) | `seqnado` | Genomics pipelines for ATAC-seq, ChIP-seq, CUT&RUN/TAG, RNA-seq, WGS, methylation, CRISPR, and MCC | `nado seq` |
| [PlotNado](https://github.com/Milne-Group/PlotNado) | `plotnado` | Genomic track visualisation from bigWig, bigBed, BED, and config files | `nado plot` |
| [TrackNado](https://github.com/Milne-Group/TrackNado) | `tracknado` | UCSC track hub generation | `nado track` |
| [BamNado](https://github.com/Milne-Group/BamNado) | `bamnado` | High-performance BAM processing with Rust/PyO3 bindings | `bamnado` |
| [QuantNado](https://github.com/Milne-Group/QuantNado) | `quantnado` | Zarr-backed genomic signal quantification and peak calling | `nado quant` |
| [MCCNado](https://github.com/alsmith151/MCCNado) | `mccnado` | Micro-Capture-C processing utilities | library/MCP |
| [TabNado](https://github.com/Milne-Group/TabNado) | `tabnado` | TF-binding prediction from epigenomic cofactors | `nado tab` |

## Requirements

- Python 3.12 or newer for `nadoverse`.
- Python 3.12 is the safest environment for all current extras.
- TabNado currently requires Python `<3.13`; `nado doctor` reports this compatibility status.

## Install With uv

Install one or more tools into the active environment:

```bash
uv pip install "nadoverse[seqnado]"
uv pip install "nadoverse[seqnado,plotnado,quantnado]"
```

Install the main published toolset:

```bash
uv pip install "nadoverse[all]"
```

`all` excludes TabNado because TabNado currently requires Python `<3.13` and its published dependency constraints are unsatisfiable on macOS x86_64. Install it separately on a compatible Python 3.12 environment:

```bash
uv pip install "nadoverse[tabnado]"
```

Install the MCP server support:

```bash
uv pip install "nadoverse[mcp]"
```

### Fresh venv

```bash
uv venv .venv --python 3.12
source .venv/bin/activate
uv pip install "nadoverse[seqnado,plotnado,quantnado]"
nado doctor
```

### Existing conda env

```bash
conda create -n nadoverse python=3.12
conda activate nadoverse
uv pip install "nadoverse[seqnado,plotnado,quantnado]"
nado doctor
```

If `uv` is not already available inside the conda environment, install it first:

```bash
conda install -c conda-forge uv
```

## CLI

```bash
nado doctor               # show tools, installed versions, and compatibility
nado seq   [args...]      # SeqNado
nado plot  [args...]      # PlotNado
nado track [args...]      # TrackNado
nado quant [args...]      # QuantNado
nado tab   [args...]      # TabNado
```

Each `nado` subcommand passes arguments through to the native tool CLI. If a tool is not installed, nadoverse prints the matching `uv pip install` command:

```text
SeqNado is not installed.
Run: uv pip install "nadoverse[seqnado]"
```

Example:

```bash
uv pip install "nadoverse[seqnado]"
nado seq --help
nado doctor
```

## MCP Server

The `nadoverse-mcp` server exposes *Nado tools for MCP hosts such as Claude Code.

```bash
uv pip install "nadoverse[mcp]"
```

Add this server to `~/.claude/mcp.json`:

```json
{
  "mcpServers": {
    "nadoverse": {
      "command": "nadoverse-mcp"
    }
  }
}
```

Core MCP tools:

| Tool | Description |
|------|-------------|
| `list_nado_tools` | Show install status, versions, CLI entrypoints, and input/output types |
| `bamnado_get_signal` | Extract per-bin coverage for one chromosome from a BAM |
| `bamnado_bam_coverage` | Generate bedGraph or BigWig coverage from one BAM |
| `bamnado_multi_bam_coverage` | Merge coverage from multiple BAM files |
| `bamnado_bigwig_compare` | Compare two BigWig files bin by bin |
| `plotnado_from_template` | Render genomic figures from PlotNado templates |
| `seqnado_generate_config` | Generate non-interactive SeqNado config files where supported |

The MCP server imports tool packages lazily. This means `nadoverse-mcp` can start when only a subset of extras is installed, and missing tools report their install command when called.

## Development

```bash
git clone https://github.com/alsmith151/nadoverse
cd nadoverse
uv sync --extra dev --extra mcp
uv run pytest
```

Run a focused check while editing:

```bash
uv run pytest tests/test_registry.py tests/test_cli.py tests/test_mcp_server.py
```

## Registry

`src/nadoverse/registry.py` is the single source of truth for tool metadata. The CLI, MCP server, and future web/API integrations use the same `NadoTool` entries for names, install extras, input/output types, repository links, and compatibility checks.

Planned additions include ScNado, GeoNado, WigNado, PeakNado, and PsudobulkNado once their public packages are available.
