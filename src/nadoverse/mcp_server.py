"""
Local stdio (or remote HTTP) MCP server for the *Nado bioinformatics toolkit.

Each tool imports its library lazily so the server starts even when only a
subset of extras are installed. Missing packages produce a clear install hint.

Usage:
    pip install nadoverse[mcp]      # adds fastmcp

    # stdio (default) — add to ~/.claude/mcp.json:
    { "mcpServers": { "nadoverse": { "command": "nadoverse-mcp" } } }

    # HTTP — set MCP_TRANSPORT=http before starting
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import Any, Optional

from nadoverse.registry import all_tools

try:
    from fastmcp import FastMCP
except ImportError:
    raise ImportError("Run: pip install nadoverse[mcp]")

mcp = FastMCP(
    name="nadoverse",
    instructions=(
        "Bioinformatics tools from the Milne group: genomics pipelines, "
        "signal quantification, genomic visualisation, track hub generation, "
        "and Micro-Capture-C analysis. Call list_nado_tools first to confirm "
        "which extras are installed."
    ),
)


def _install_hint(pypi_name: str, extra: str) -> str:
    return f"{pypi_name} is not installed. Run: pip install nadoverse[{extra}]"


def _run(cmd: list[str]) -> str:
    """Run a CLI command, return stdout, raise RuntimeError on non-zero exit."""
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or f"Command failed: {' '.join(cmd)}")
    return result.stdout.strip()


# ─── Registry ────────────────────────────────────────────────────────────────

@mcp.tool(
    description=(
        "List all *Nado tools: install status, version, CLI entrypoints, "
        "input/output types. Call this first to know what is available."
    ),
    annotations={"readOnlyHint": True, "destructiveHint": False},
)
def list_nado_tools() -> str:
    rows = [
        {
            "name": t.name,
            "installed": t.is_available(),
            "version": t.installed_version(),
            "cli_command": t.cli_command,
            "install_hint": f"pip install nadoverse[{t.install_extra}]",
            "input_types": t.input_types,
            "output_types": t.output_types,
            "python_compatible": t.python_compatible(),
            "description": t.description,
        }
        for t in all_tools()
    ]
    return json.dumps(rows, indent=2)


# ─── BamNado ─────────────────────────────────────────────────────────────────

@mcp.tool(
    description="Extract per-bin coverage signal for one chromosome from a BAM file.",
    annotations={"readOnlyHint": True, "destructiveHint": False},
)
def bamnado_get_signal(
    bam_path: str,
    chromosome: str,
    bin_size: int = 50,
    scale_factor: float = 1.0,
    use_fragment: bool = False,
    ignore_scaffold_chromosomes: bool = True,
) -> dict[str, Any]:
    """Returns dict with 'chromosome', 'bin_size', 'n_bins', and 'signal' (list[float])."""
    try:
        from bamnado import get_signal_for_chromosome
    except ImportError:
        raise RuntimeError(_install_hint("bamnado", "bamnado"))

    arr = get_signal_for_chromosome(
        bam_path=bam_path,
        chromosome_name=chromosome,
        bin_size=bin_size,
        scale_factor=scale_factor,
        use_fragment=use_fragment,
        ignore_scaffold_chromosomes=ignore_scaffold_chromosomes,
    )
    return {
        "chromosome": chromosome,
        "bin_size": bin_size,
        "n_bins": len(arr),
        "signal": arr.tolist(),
    }


# ─── MCCNado ─────────────────────────────────────────────────────────────────

@mcp.tool(
    description="Annotate a Micro-Capture-C BAM file with viewpoint and fragment information.",
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def mccnado_annotate_bam(bam: str, output: str) -> str:
    """bam must be sorted by query name."""
    return _run(["mccnado", "annotate-bam", bam, output])


@mcp.tool(
    description="Remove PCR duplicates from an annotated MCC BAM file.",
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def mccnado_deduplicate_bam(bam: str, output: str) -> str:
    return _run(["mccnado", "deduplicate-bam", bam, output])


@mcp.tool(
    description="Split MCC BAM reads by viewpoint into separate BAM files.",
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def mccnado_split_viewpoint_reads(bam: str, output: str) -> str:
    return _run(["mccnado", "split-viewpoint-reads", bam, output])


@mcp.tool(
    description="Compute ligation junction statistics from an MCC BAM file.",
    annotations={"readOnlyHint": True, "destructiveHint": False},
)
def mccnado_extract_ligation_stats(bam: str, stats_output: str) -> str:
    return _run(["mccnado", "extract-ligation-stats", bam, stats_output])


@mcp.tool(
    description="Identify ligation junctions in an MCC BAM and write per-viewpoint output.",
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def mccnado_identify_ligation_junctions(bam: str, outdir: str) -> str:
    return _run(["mccnado", "identify-ligation-junctions", bam, outdir])


@mcp.tool(
    description="Deduplicate paired-end FASTQ files before MCC alignment.",
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def mccnado_deduplicate_fastq(
    fastq_1: str,
    fastq_2: str,
    output_1: str,
    output_2: str,
    stats_prefix: str,
) -> str:
    return _run([
        "mccnado", "deduplicate-fastq",
        fastq_1, fastq_2, output_1, output_2,
        "--stats-prefix", stats_prefix,
    ])


@mcp.tool(
    description="Merge multiple per-viewpoint ligation junction cooler files into one.",
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def mccnado_combine_coolers(cooler_files: list[str], output: str) -> str:
    return _run(
        ["mccnado", "combine-ligation-junction-coolers"] + cooler_files + ["--output", output]
    )


# ─── PlotNado ────────────────────────────────────────────────────────────────

@mcp.tool(
    description="Render a genomic figure from a plotnado YAML template file.",
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def plotnado_from_template(
    template_yaml: str,
    region: str,
    output_path: str,
    dpi: int = 300,
    theme: str = "publication",
) -> str:
    """region: 'chrN:start-end'. output_path extension sets format (PNG/SVG/PDF)."""
    try:
        from plotnado.figure import GenomicFigure
    except ImportError:
        raise RuntimeError(_install_hint("plotnado", "plotnado"))

    fig = GenomicFigure.from_template(template_yaml, theme=theme)
    fig.save(output_path, region=region, dpi=dpi)
    return f"Saved to {output_path}"


@mcp.tool(
    description="Render a genomic figure from a saved IGV session XML file.",
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def plotnado_from_igv_session(
    session_xml: str,
    output_path: str,
    region: Optional[str] = None,
    dpi: int = 300,
    theme: str = "publication",
) -> str:
    """Uses the session's stored locus if region is not provided."""
    try:
        from plotnado.figure import GenomicFigure
    except ImportError:
        raise RuntimeError(_install_hint("plotnado", "plotnado"))

    fig, session_locus = GenomicFigure.from_igv_session(session_xml, theme=theme)
    locus = region or session_locus
    if locus is None:
        raise ValueError("No region in session and none provided")
    fig.save(output_path, region=locus, dpi=dpi)
    return f"Saved to {output_path}"


@mcp.tool(
    description=(
        "Build and render a genomic figure from a list of track specifications. "
        "Each track is a dict with a 'type' key plus track-specific params. "
        "Common types: 'bigwig', 'bed', 'narrowpeak', 'genes', 'axis', 'scalebar'."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def plotnado_build(
    tracks: list[dict[str, Any]],
    region: str,
    output_path: str,
    dpi: int = 300,
    theme: str = "publication",
    width: float = 12,
    highlight_regions: Optional[list[str]] = None,
) -> str:
    try:
        from plotnado.figure import GenomicFigure
    except ImportError:
        raise RuntimeError(_install_hint("plotnado", "plotnado"))

    fig = GenomicFigure(width=width, theme=theme)
    for track_spec in tracks:
        spec = dict(track_spec)
        track_type = spec.pop("type")
        fig.add_track(track_type, **spec)
    if highlight_regions:
        for hr in highlight_regions:
            fig.highlight(hr)
    fig.save(output_path, region=region, dpi=dpi)
    return f"Saved to {output_path}"


@mcp.tool(
    description="Infer a plotnado YAML template from a list of track files.",
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def plotnado_init_template(
    input_files: list[str],
    output_yaml: str,
    genome: str = "hg38",
) -> str:
    """Detects file types (BigWig, BED, narrowPeak, etc.) and writes a ready-to-edit template."""
    return _run(
        ["plotnado", "init", "--output", output_yaml, "--genome", genome] + input_files
    )


# ─── QuantNado ───────────────────────────────────────────────────────────────

@mcp.tool(
    description=(
        "Create a QuantNado Zarr dataset from BAM files. "
        "Quantifies binned coverage and writes a Zarr v3 store."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def quantnado_create_dataset(
    bam_files: list[str],
    output_zarr: str,
    chrom_sizes: Optional[str] = None,
    bin_size: int = 200,
    n_workers: int = 4,
) -> str:
    try:
        from quantnado import QuantNado
    except ImportError:
        raise RuntimeError(_install_hint("quantnado", "quantnado"))

    qn = QuantNado.from_bam_files(
        bam_files=bam_files,
        output_path=output_zarr,
        chrom_sizes=chrom_sizes,
        bin_size=bin_size,
        n_workers=n_workers,
    )
    return f"Dataset created at {output_zarr} with {len(qn.samples)} samples: {qn.samples}"


@mcp.tool(
    description="Return metadata about an existing QuantNado Zarr dataset.",
    annotations={"readOnlyHint": True, "destructiveHint": False},
)
def quantnado_dataset_info(dataset_path: str) -> dict[str, Any]:
    try:
        from quantnado import QuantNado
    except ImportError:
        raise RuntimeError(_install_hint("quantnado", "quantnado"))

    qn = QuantNado.open_dataset(dataset_path)
    return {
        "samples": qn.samples,
        "chromosomes": qn.chromosomes,
        "modalities": qn.modalities,
        "n_completed": qn.n_completed,
        "chromsizes": qn.chromsizes,
    }


@mcp.tool(
    description="Extract coverage signal for a genomic region and save to CSV.",
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def quantnado_extract_region(
    dataset_path: str,
    region: str,
    output_csv: str,
    normalise: Optional[str] = None,
    samples: Optional[list[str]] = None,
) -> str:
    """region: 'chrN:start-end'. normalise: 'cpm' or 'rpkm'."""
    try:
        from quantnado import QuantNado
    except ImportError:
        raise RuntimeError(_install_hint("quantnado", "quantnado"))

    qn = QuantNado.open_dataset(dataset_path)
    arr = qn.extract_region(region=region, normalise=normalise, samples=samples, as_xarray=True)
    df = arr.to_pandas()
    df.to_csv(output_csv)
    return f"Extracted {region} ({arr.shape}) to {output_csv}"


@mcp.tool(
    description="Call peaks from a QuantNado dataset and write results to BED.",
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def quantnado_call_peaks(
    dataset_path: str,
    output_bed: str,
    method: str = "quantile",
    sample_name: Optional[str] = None,
    blacklist_file: Optional[str] = None,
) -> str:
    """method: 'quantile', 'seacr', 'lanceotron', 'macs3', 'macs3_broad', 'unet', or 'resnet'."""
    try:
        from quantnado import QuantNado
    except ImportError:
        raise RuntimeError(_install_hint("quantnado", "quantnado"))

    qn = QuantNado.open_dataset(dataset_path)
    peaks = qn.call_peaks(
        method=method,
        sample_name=sample_name,
        blacklist_file=Path(blacklist_file) if blacklist_file else None,
    )
    peaks.to_bed(output_bed)
    return f"Called peaks with {method}: {len(peaks)} peaks → {output_bed}"


@mcp.tool(
    description="Generate a metaplot (average signal over genomic intervals) and save to PNG.",
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def quantnado_metaplot(
    dataset_path: str,
    intervals_path: str,
    output_png: str,
    normalise: Optional[str] = "cpm",
    samples: Optional[list[str]] = None,
) -> str:
    try:
        from quantnado import QuantNado
    except ImportError:
        raise RuntimeError(_install_hint("quantnado", "quantnado"))

    qn = QuantNado.open_dataset(dataset_path)
    fig = qn.metaplot(intervals_path=intervals_path, normalise=normalise, samples=samples)
    fig.savefig(output_png, dpi=150, bbox_inches="tight")
    return f"Metaplot saved to {output_png}"


@mcp.tool(
    description="Run PCA on sample coverage and save a scatter plot to PNG.",
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def quantnado_pca(
    dataset_path: str,
    output_png: str,
    normalise: Optional[str] = "cpm",
    samples: Optional[list[str]] = None,
) -> str:
    try:
        from quantnado import QuantNado
    except ImportError:
        raise RuntimeError(_install_hint("quantnado", "quantnado"))

    qn = QuantNado.open_dataset(dataset_path)
    fig = qn.pca(normalise=normalise, samples=samples)
    fig.savefig(output_png, dpi=150, bbox_inches="tight")
    return f"PCA plot saved to {output_png}"


@mcp.tool(
    description="Plot coverage signal at a locus for all samples and save to PNG.",
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def quantnado_locus_plot(
    dataset_path: str,
    region: str,
    output_png: str,
    normalise: Optional[str] = "cpm",
    samples: Optional[list[str]] = None,
) -> str:
    try:
        from quantnado import QuantNado
    except ImportError:
        raise RuntimeError(_install_hint("quantnado", "quantnado"))

    qn = QuantNado.open_dataset(dataset_path)
    fig = qn.locus_plot(region=region, normalise=normalise, samples=samples)
    fig.savefig(output_png, dpi=150, bbox_inches="tight")
    return f"Locus plot saved to {output_png}"


# ─── TrackNado ───────────────────────────────────────────────────────────────

@mcp.tool(
    description=(
        "Generate a UCSC track hub from track files and optional metadata CSV. "
        "Either input_files or metadata_csv must be provided."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def tracknado_create(
    output_dir: str,
    input_files: Optional[list[str]] = None,
    metadata_csv: Optional[str] = None,
    genome: str = "hg38",
    hub_name: str = "HUB",
    hub_email: str = "alastair.smith@ndcls.ox.ac.uk",
    color_by: Optional[str] = None,
    supergroup_by: Optional[list[str]] = None,
    subgroup_by: Optional[list[str]] = None,
    url_prefix: str = "https://userweb.molbiol.ox.ac.uk",
    seqnado_layout: bool = False,
) -> str:
    args = [
        "tracknado", "create",
        "--output", output_dir,
        "--genome-name", genome,
        "--hub-name", hub_name,
        "--hub-email", hub_email,
        "--url-prefix", url_prefix,
    ]
    if input_files:
        for f in input_files:
            args += ["--input-files", f]
    if metadata_csv:
        args += ["--metadata", metadata_csv]
    if color_by:
        args += ["--color-by", color_by]
    if supergroup_by:
        for g in supergroup_by:
            args += ["--supergroup-by", g]
    if subgroup_by:
        for g in subgroup_by:
            args += ["--subgroup-by", g]
    if seqnado_layout:
        args.append("--seqnado")
    return _run(args)


# ─── SeqNado ─────────────────────────────────────────────────────────────────

@mcp.tool(
    description="List all supported SeqNado assay types and pipeline options.",
    annotations={"readOnlyHint": True, "destructiveHint": False},
)
def seqnado_list_assays() -> dict[str, Any]:
    try:
        from seqnado import (
            Assay, AssaysWithPeakCalling, AssaysWithSpikein, AssaysWithHeatmaps,
            PileupMethod, PeakCallingMethod, PCRDuplicateHandling, DataScalingTechnique,
        )
    except ImportError:
        raise RuntimeError(_install_hint("seqnado", "seqnado"))

    return {
        "assays": {
            a.clean_name: {
                "full_name": a.value,
                "peak_calling": a in AssaysWithPeakCalling,
                "spikein": a in AssaysWithSpikein,
                "heatmaps": a in AssaysWithHeatmaps,
                "ip_required": a in Assay.ip_assays(),
            }
            for a in Assay
            if a != Assay.MULTIOMICS
        },
        "pileup_methods": [m.value for m in PileupMethod],
        "peak_calling_methods": [m.value for m in PeakCallingMethod],
        "pcr_duplicate_handling": [m.value for m in PCRDuplicateHandling],
        "scaling_techniques": [m.value for m in DataScalingTechnique],
    }


@mcp.tool(
    description="List all genome configurations registered in SeqNado for a given assay.",
    annotations={"readOnlyHint": True, "destructiveHint": False},
)
def seqnado_list_genomes(assay: str = "chip") -> dict[str, Any]:
    try:
        from seqnado import Assay
        from seqnado.config.user_input import load_genome_configs
    except ImportError:
        raise RuntimeError(_install_hint("seqnado", "seqnado"))

    assay_enum = Assay.from_clean_name(assay.lower())
    genomes = load_genome_configs(assay_enum)
    return {
        name: {k: str(v) if v is not None else None for k, v in cfg.model_dump().items()}
        for name, cfg in genomes.items()
    }


@mcp.tool(
    description=(
        "Generate a SeqNado design CSV from FASTQ file paths. "
        "Parses Illumina-style filenames to extract sample names, replicates, and read pairs."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def seqnado_generate_design(
    assay: str,
    fastq_files: list[str],
    output_csv: str,
    ip_to_control: Optional[str] = None,
    group_by: Optional[str] = None,
    deseq2_pattern: Optional[str] = None,
) -> str:
    """assay: 'chip', 'atac', 'rna', 'cat', 'meth', 'snp', 'mcc', 'crispr'."""
    try:
        from seqnado import Assay  # noqa: F401 — import check only
    except ImportError:
        raise RuntimeError(_install_hint("seqnado", "seqnado"))

    args = [
        "seqnado", "design", assay,
        "--output", output_csv,
        "--no-interactive",
        "--accept-all-defaults",
    ] + fastq_files
    if ip_to_control:
        args += ["--ip-to-control", ip_to_control]
    if group_by:
        args += ["--group-by", group_by]
    if deseq2_pattern:
        args += ["--deseq2-pattern", deseq2_pattern]
    return _run(args)


@mcp.tool(
    description=(
        "Generate a SeqNado workflow config YAML using the Python API. "
        "The output file can be edited before passing to seqnado_run_pipeline."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def seqnado_generate_config(
    assay: str,
    genome: str,
    output_yaml: str,
    project_name: str = "seqnado_project",
    pileup_method: str = "deeptools",
    peak_calling_method: str = "lanceotron",
    pcr_duplicates: str = "markdup",
    bin_size: int = 10,
    create_heatmaps: bool = False,
    spikein_genome: Optional[str] = None,
    run_deseq2: bool = False,
    strandedness: int = 0,
    ucsc_hub_dir: str = "seqnado_output/hub/",
    ucsc_hub_email: str = "alastair.smith@ndcls.ox.ac.uk",
    ucsc_hub_url: str = "https://userweb.molbiol.ox.ac.uk",
) -> str:
    try:
        from importlib import resources as _resources
        from seqnado import (
            Assay, PCRDuplicateHandling, PeakCallingMethod,
            PileupMethod, QuantificationMethod,
        )
        from seqnado.config.user_input import (
            build_default_assay_config, load_genome_configs, render_config,
        )
        from seqnado.config.configs import (
            BigwigConfig, PCRDuplicatesConfig, PeakCallingConfig, ProjectConfig,
            RNAQuantificationConfig, SpikeInConfig, UCSCHubConfig,
        )
    except ImportError:
        raise RuntimeError(_install_hint("seqnado", "seqnado"))

    assay_enum = Assay.from_clean_name(assay.lower())
    genomes = load_genome_configs(assay_enum)
    if genome not in genomes:
        raise ValueError(
            f"Genome '{genome}' not found. Available: {list(genomes)}. "
            "Run seqnado_list_genomes to see options."
        )
    genome_cfg = genomes[genome]
    assay_cfg = build_default_assay_config(assay_enum, genome_cfg)
    if assay_cfg is None:
        raise ValueError(f"Could not build default config for assay '{assay}'")

    assay_cfg.bigwigs = BigwigConfig(
        pileup_method=[PileupMethod(pileup_method)], binsize=bin_size,
    )
    assay_cfg.create_heatmaps = create_heatmaps
    assay_cfg.ucsc_hub = UCSCHubConfig(
        directory=ucsc_hub_dir, genome=genome, email=ucsc_hub_email,
        genome_name=genome, url=ucsc_hub_url,
    )
    if hasattr(assay_cfg, "peak_calling") and assay_cfg.peak_calling is not None:
        assay_cfg.peak_calling = PeakCallingConfig(
            method=[PeakCallingMethod(peak_calling_method)], consensus_counts=False,
        )
    if hasattr(assay_cfg, "spikein") and spikein_genome:
        spikein_genomes = load_genome_configs(assay_enum)
        if spikein_genome in spikein_genomes:
            assay_cfg.spikein = SpikeInConfig(genome=spikein_genomes[spikein_genome])
    if hasattr(assay_cfg, "rna_quantification") and assay_cfg.rna_quantification is not None:
        assay_cfg.rna_quantification = RNAQuantificationConfig(
            method=QuantificationMethod.FEATURE_COUNTS,
            run_deseq2=run_deseq2,
            strandedness=strandedness,
        )

    from seqnado._version import __version__ as seqnado_version
    pkg_root = _resources.files("seqnado")
    template_path = pkg_root / "config" / "templates" / f"config_{assay.lower()}.yaml"
    pcr_dup_cfg = PCRDuplicatesConfig(strategy=PCRDuplicateHandling(pcr_duplicates))
    from seqnado.config.user_input import SeqnadoConfig  # type: ignore[attr-defined]
    workflow_cfg = SeqnadoConfig(
        project=ProjectConfig(project_name=project_name),
        pcr_duplicates=pcr_dup_cfg,
        **assay_cfg.model_dump(),
    )
    render_config(
        template=Path(str(template_path)),
        workflow_config=workflow_cfg,
        outfile=Path(output_yaml),
        seqnado_version=seqnado_version,
    )
    return f"Config written to {output_yaml}"


@mcp.tool(
    description=(
        "Run (or dry-run) a SeqNado Snakemake pipeline. "
        "Always dry_run=True first to preview the job graph before committing."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def seqnado_run_pipeline(
    assay: str,
    config_yaml: str,
    working_dir: str,
    cores: int = 4,
    dry_run: bool = True,
    targets: Optional[list[str]] = None,
    profile: Optional[str] = None,
    extra_snakemake_args: Optional[list[str]] = None,
) -> str:
    """dry_run defaults to True — call with dry_run=False only after reviewing the plan."""
    try:
        from seqnado import Assay  # noqa: F401
    except ImportError:
        raise RuntimeError(_install_hint("seqnado", "seqnado"))

    args = [
        "seqnado", assay,
        "--configfile", config_yaml,
        "--cores", str(cores),
        "--directory", working_dir,
    ]
    if dry_run:
        args.append("--dry-run")
    if profile:
        args += ["--profile", profile]
    if targets:
        args += targets
    if extra_snakemake_args:
        args += extra_snakemake_args

    result = subprocess.run(args, capture_output=True, text=True, cwd=working_dir)
    if result.returncode != 0:
        raise RuntimeError(
            f"SeqNado failed (exit {result.returncode}):\n"
            f"STDOUT:\n{result.stdout[-3000:]}\n"
            f"STDERR:\n{result.stderr[-3000:]}"
        )
    return result.stdout.strip() or "Pipeline completed successfully."


@mcp.tool(
    description=(
        "Check the completion status of a SeqNado pipeline run. "
        "Inspects seqnado_output/ for BAMs, BigWigs, peaks, QC reports, and trackhubs."
    ),
    annotations={"readOnlyHint": True, "destructiveHint": False},
)
def seqnado_pipeline_status(working_dir: str) -> dict[str, Any]:
    output_dir = Path(working_dir) / "seqnado_output"
    if not output_dir.exists():
        return {"status": "not_started", "output_dir": str(output_dir)}

    checks = {
        "bams": list(output_dir.glob("**/*.bam")),
        "bigwigs": list(output_dir.glob("**/*.bw")),
        "peaks": list(output_dir.glob("**/*.bed")) + list(output_dir.glob("**/*.narrowPeak")),
        "qc_multiqc": list(output_dir.glob("**/multiqc_report.html")),
        "hub": list((output_dir / "hub").glob("hub.txt")) if (output_dir / "hub").exists() else [],
        "logs": (
            list((Path(working_dir) / ".snakemake" / "log").glob("*.log"))
            if (Path(working_dir) / ".snakemake" / "log").exists() else []
        ),
    }
    steps = {
        k: {"count": len(v), "files": [f.name for f in v[:10]]}
        for k, v in checks.items()
    }
    total = sum(s["count"] for s in steps.values())
    return {
        "output_dir": str(output_dir),
        "status": "completed" if total > 0 else "empty",
        "steps": steps,
    }


@mcp.tool(
    description="Validate a SeqNado design CSV against the schema for the given assay.",
    annotations={"readOnlyHint": True, "destructiveHint": False},
)
def seqnado_validate_design(design_csv: str, assay: str) -> dict[str, Any]:
    try:
        import pandas as pd
        from seqnado import Assay
    except ImportError:
        raise RuntimeError(_install_hint("seqnado", "seqnado"))

    assay_enum = Assay.from_clean_name(assay.lower())
    df = pd.read_csv(design_csv)
    errors: list[str] = []
    required = ["sample_id", "fastq_1", "fastq_2"]
    if assay_enum in Assay.ip_assays():
        required += ["ip", "control"]
    for col in required:
        if col not in df.columns:
            errors.append(f"Missing required column: '{col}'")
    for col in ["fastq_1", "fastq_2"]:
        if col in df.columns:
            missing = [str(p) for p in df[col] if not Path(str(p)).exists()]
            if missing:
                errors.append(f"Missing FASTQ files in '{col}': {missing[:5]}")
    return {
        "valid": len(errors) == 0,
        "errors": errors,
        "columns": list(df.columns),
        "n_samples": len(df),
        "samples": df["sample_id"].tolist() if "sample_id" in df.columns else [],
    }


# ─── Entry point ─────────────────────────────────────────────────────────────

def main() -> None:
    transport = os.getenv("MCP_TRANSPORT", "stdio")
    if transport == "http":
        mcp.run(
            transport="streamable-http",
            host="0.0.0.0",
            port=int(os.getenv("PORT", "8080")),
        )
    else:
        mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
