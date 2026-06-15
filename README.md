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
| [ReguloNado](https://github.com/alsmith151/ReguloNado) | `regulonado` | Dataset building and training for large sequence-to-function models from FASTA, BED, and BigWig inputs | `nado regulo` |

## Requirements

- Python 3.12 or newer for `nadoverse`.
- Python 3.12 is the safest environment for all current extras.
- TabNado currently requires Python `<3.13`; `nado doctor` reports this compatibility status.
- ReguloNado currently requires Python `<3.13`; `nado doctor` reports this compatibility status.

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

`all` excludes TabNado and ReguloNado because they currently require Python `<3.13`; ReguloNado is also installed directly from GitHub and pulls large model-training dependencies. Install them separately on a compatible Python 3.12 environment:

```bash
uv pip install "nadoverse[tabnado]"
uv pip install "nadoverse[regulonado]"
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
nado seq    [args...]     # SeqNado
nado plot   [args...]     # PlotNado
nado track  [args...]     # TrackNado
nado quant  [args...]     # QuantNado
nado tab    [args...]     # TabNado
nado regulo [args...]     # ReguloNado
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

Available MCP tools:

| Tool | Description |
|------|-------------|
| `list_nado_tools` | Show install status, versions, CLI entrypoints, and input/output types |

### BamNado MCP Tools

| Tool | Description |
|------|-------------|
| `bamnado_get_signal` | Extract per-bin coverage for one chromosome from a BAM |
| `bamnado_bam_coverage` | Generate bedGraph or BigWig coverage from one BAM |
| `bamnado_multi_bam_coverage` | Merge coverage from multiple BAM files |
| `bamnado_bigwig_compare` | Compare two BigWig files bin by bin |
| `bamnado_bigwig_aggregate` | Aggregate multiple BigWig files with sum, mean, median, max, or min |
| `bamnado_collapse_bedgraph` | Collapse adjacent equal-score bedGraph bins |
| `bamnado_split` | Split or filter a BAM with BamNado read filters |
| `bamnado_split_exogenous` | Split endogenous and exogenous reads by reference-name prefix |
| `bamnado_modify` | Filter and adjust reads, including optional Tn5 shifting |
| `bamnado_bigwig_infer_scale` | Infer scaling factors and library size from CPM/RPKM-normalised BigWig files |

### MCCNado MCP Tools

| Tool | Description |
|------|-------------|
| `mccnado_deduplicate_fastq` | Deduplicate paired-end FASTQ reads before alignment |
| `mccnado_annotate_bam` | Annotate aligned MCC BAM files with viewpoint/on-capture/reporter tags |
| `mccnado_deduplicate_bam` | Remove PCR duplicates from annotated MCC BAM files |
| `mccnado_split_viewpoint_reads` | Split deduplicated MCC BAM files by viewpoint |
| `mccnado_identify_ligation_junctions` | Write per-viewpoint cooler files from MCC ligation junctions |
| `mccnado_extract_ligation_stats` | Generate MCC ligation-junction statistics |
| `mccnado_combine_coolers` | Merge per-viewpoint cooler files into one multi-viewpoint cooler |

### PlotNado MCP Tools

| Tool | Description |
|------|-------------|
| `plotnado_from_template` | Render genomic figures from PlotNado templates |
| `plotnado_from_igv_session` | Render figures from IGV session XML files |
| `plotnado_build` | Build a plot directly from track inputs and region arguments |
| `plotnado_init_template` | Generate a starter PlotNado YAML template |

### QuantNado MCP Tools

| Tool | Description |
|------|-------------|
| `quantnado_create_dataset` | Create a QuantNado dataset from BAM or BigWig files |
| `quantnado_dataset_info` | Inspect dataset metadata, samples, and intervals |
| `quantnado_extract_region` | Extract signal values for a genomic region |
| `quantnado_call_peaks` | Call peaks with quantile or SEACR-style methods |
| `quantnado_metaplot` | Generate aggregate signal profiles around genomic features |
| `quantnado_pca` | Run PCA over a QuantNado dataset |
| `quantnado_locus_plot` | Plot signal over one locus |

### ReguloNado MCP Tools

| Tool | Description |
|------|-------------|
| `regulonado_build_dataset` | Build Arrow training datasets from BED, FASTA, and BigWig inputs |
| `regulonado_train` | Launch or dry-run sequence-to-function model training |
| `regulonado_scale_bigwigs` | Infer scale factors for BigWig files |
| `regulonado_calculate_original_scaling` | Compute RPKM-to-raw-count scale factors from dataset metadata |
| `regulonado_calculate_tmm_scaling` | Compute TMM normalisation factors |
| `regulonado_enrich_metadata` | Write scaling fields into ReguloNado metadata |
| `regulonado_recompress_dataset` | Rechunk/recompress saved Arrow datasets |

### TrackNado MCP Tools

| Tool | Description |
|------|-------------|
| `tracknado_create` | Create a UCSC track hub from sequencing output files |

### SeqNado MCP Tools

| Tool | Description |
|------|-------------|
| `seqnado_list_assays` | List supported SeqNado assay types |
| `seqnado_list_genomes` | List available genomes for a SeqNado assay |
| `seqnado_generate_design` | Generate a design CSV template |
| `seqnado_generate_config` | Generate non-interactive SeqNado config files where supported |
| `seqnado_run_pipeline` | Run or dry-run a SeqNado pipeline command |
| `seqnado_pipeline_status` | Inspect output directories for common pipeline outputs |
| `seqnado_validate_design` | Validate a SeqNado design CSV |
| `seqnado_download` | Download public sequencing data with SeqNado helpers |
| `seqnado_build_genome` | Build a genome resource for SeqNado workflows |

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
