# nadoverse

Unified installer and CLI for the **\*Nado** bioinformatics toolkit.

```
pip install nadoverse[seqnado,plotnado,quantnado]
nado doctor
```

---

## Included tools

| Tool | Extra | Description | PyPI |
|------|-------|-------------|------|
| [SeqNado](https://github.com/Milne-Group/SeqNado) | `seqnado` | Genomics pipelines: ATAC-seq, ChIP-seq, CUT&RUN/TAG, RNA-seq, WGS, Methylation, CRISPR, MCC | ✅ |
| [PlotNado](https://github.com/Milne-Group/PlotNado) | `plotnado` | Genomic track visualisation from bigWig/bigBed | ✅ |
| [TrackNado](https://github.com/Milne-Group/TrackNado) | `tracknado` | UCSC trackhub generator | ✅ |
| [BamNado](https://github.com/Milne-Group/BamNado) | `bamnado` | High-performance BAM processing (Rust + PyO3, library only) | ✅ |
| [QuantNado](https://github.com/Milne-Group/QuantNado) | `quantnado` | Zarr-backed genomic signal quantification and peak calling | ✅ |
| [MCCNado](https://github.com/alsmith151/MCCNado) | `mccnado` | Micro-Capture-C processing utilities (Rust, library only) | ✅ |
| [TabNado](https://github.com/Milne-Group/TabNado) | `tabnado` | TF-binding prediction from epigenomic cofactors ⚠️ Python <3.13 | ✅ |

---

## Install

```bash
# Individual tools
pip install nadoverse[seqnado]
pip install nadoverse[seqnado,plotnado,quantnado]

# Everything
pip install nadoverse[all]

# MCP server (Claude Code integration)
pip install nadoverse[mcp]
```

Requires Python ≥ 3.12. TabNado additionally requires Python < 3.13; `nado doctor` warns if your version exceeds that cap.

---

## CLI — `nado`

```
nado doctor               # show all tools, versions, and status
nado seq   [args...]      # SeqNado
nado plot  [args...]      # PlotNado
nado track [args...]      # TrackNado
nado quant [args...]      # QuantNado
nado tab   [args...]      # TabNado
```

Each subcommand passes all arguments through to the native tool CLI. If a tool is not installed you get an install hint instead of an error:

```
SeqNado is not installed.
Run: pip install nadoverse[seqnado]
```

### Example

```bash
pip install nadoverse[seqnado]
nado seq --help            # same as: seqnado --help
nado doctor                # confirm everything is wired up
```

---

## MCP server — Claude Code integration

The `nadoverse-mcp` server exposes each tool as an MCP tool for use in
Claude Code agentic workflows.

```bash
pip install nadoverse[mcp]
```

Add to `~/.claude/mcp.json`:

```json
{
  "mcpServers": {
    "nadoverse": {
      "command": "nadoverse-mcp"
    }
  }
}
```

Available MCP tools:

| Tool | Description |
|------|-------------|
| `list_nado_tools` | Show all tools, install status, versions |
| `run_seqnado` | Run SeqNado pipeline |
| `run_plotnado` | Generate genomic track plots |
| `run_tracknado` | Generate UCSC trackhubs |
| `run_quantnado` | Quantify genomic signal |
| `run_tabnado` | Run TF-binding prediction |

All execution tools default to `dry_run=True` — Claude previews the
command and you confirm before anything runs.

---

## Roadmap

> Notes only — not yet built.

### Registry as MCP tool source of truth

`src/nadoverse/registry.py` is the single source of truth for all tool
metadata. Each `NadoTool` entry maps 1:1 to an MCP tool definition via
`to_dict()`. Future work:

- **MCP server (hosted):** A remote HTTP FastMCP server backed by the
  registry, enabling nadoverse tools in Claude.ai and other MCP hosts
  (not just Claude Code).
- **FastAPI backend:** `GET /tools` returns
  `[tool.to_dict() for tool in REGISTRY]` — ready for a web dashboard.
- **JS web UI:** Form definitions auto-generated from `input_types` and
  `output_types` fields; no schema redesign needed.

The `to_dict()` / `to_json()` methods are stable by design so the same
registry can feed the CLI, the MCP server, FastAPI responses, and UI
form schemas without duplication.

### Planned additions

- ScNado (single-cell CUT&TAG + RNA pipeline) — pending PyPI publication
- GeoNado (GEO/SRA fastq download) — pending PyPI publication
- WigNado, peaknado, psudobulknado — pending public release

---

## Development

```bash
git clone https://github.com/alsmith151/nadoverse
cd nadoverse
pip install -e ".[dev]"
pytest
```
