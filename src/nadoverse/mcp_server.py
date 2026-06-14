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
import shutil
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


def _require_cli(command: str, pypi_name: str, extra: str) -> None:
    """Raise a clear install hint when a required CLI command is unavailable."""
    if shutil.which(command) is None:
        raise RuntimeError(_install_hint(pypi_name, extra))


def _add_read_filter_args(
    args: list[str],
    *,
    strand: str = "both",
    proper_pairs: bool = False,
    min_mapq: int = 20,
    min_length: int = 20,
    max_length: int = 1000,
    min_fragment_len: Optional[int] = None,
    max_fragment_len: Optional[int] = None,
    blacklist: Optional[str] = None,
    barcode_allowlist: Optional[str] = None,
    read_group: Optional[str] = None,
    tag: Optional[str] = None,
    tag_value: Optional[str] = None,
) -> None:
    args += [
        "--strand", strand,
        "--min-mapq", str(min_mapq),
        "--min-length", str(min_length),
        "--max-length", str(max_length),
    ]
    if proper_pairs:
        args.append("--proper-pairs")
    if min_fragment_len is not None:
        args += ["--min-fragment-len", str(min_fragment_len)]
    if max_fragment_len is not None:
        args += ["--max-fragment-len", str(max_fragment_len)]
    if blacklist:
        args += ["--blacklist", blacklist]
    if barcode_allowlist:
        args += ["--barcode-allowlist", barcode_allowlist]
    if read_group:
        args += ["--read-group", read_group]
    if tag:
        args += ["--tag", tag]
    if tag_value:
        args += ["--tag-value", tag_value]


def _add_coverage_args(
    args: list[str],
    *,
    bin_size: Optional[int] = 50,
    normalize: str = "raw",
    scale_factor: Optional[float] = None,
    fragment_counts: bool = False,
    shift: str = "0,0,0,0",
    truncate: Optional[str] = None,
    ignore_scaffolds: bool = False,
    threads: int = 6,
) -> None:
    args += ["--normalize", normalize, "--threads", str(threads)]
    if bin_size is not None:
        args += ["--bin-size", str(bin_size)]
    if scale_factor is not None:
        args += ["--scale-factor", str(scale_factor)]
    if fragment_counts:
        args.append("--fragment-counts")
    if shift != "0,0,0,0":
        args += ["--shift", shift]
    if truncate:
        args += ["--truncate", truncate]
    if ignore_scaffolds:
        args.append("--ignore-scaffolds")


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
    description=(
        "Extract per-bin coverage signal for one chromosome from a sorted, indexed BAM file. "
        "Returns a float32 array of length chrom_size // bin_size. "
        "Use scale_factor = 1e6 / total_mapped_reads for CPM normalisation. "
        "Call list_nado_tools first to confirm bamnado is installed."
    ),
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
    """
    Args:
        bam_path: Absolute path to a sorted, indexed BAM file. A corresponding
            .bai index file must exist alongside it (run `samtools index` first).
        chromosome: Exact contig name as it appears in the BAM header. Use
            `samtools view -H sample.bam | grep '^@SQ'` to list valid names.
            UCSC-style uses 'chr1'; Ensembl-style uses '1'.
        bin_size: Coverage bin width in base pairs. Smaller values give finer
            resolution but produce a larger array. 50 for fine ChIP signal,
            200 for broad chromatin domains.
        scale_factor: Linear multiplier applied to every bin value after counting.
            Use 1e6 / total_mapped_reads for CPM normalisation; default 1.0
            returns raw read counts.
        use_fragment: If True, count paired-end fragment spans (insert size)
            rather than individual read alignments. Use True for ATAC-seq
            fragment-mode coverage; False for single-end or read-level counting.
        ignore_scaffold_chromosomes: If True, skip non-standard contigs such as
            chrUn_*, *_random, and EBV. Recommended True for standard analyses.

    Returns:
        Dict with 'chromosome' (str), 'bin_size' (int), 'n_bins' (int),
        and 'signal' (list[float]) — one float per bin across the chromosome.
    """
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


@mcp.tool(
    description=(
        "Generate a bedGraph or BigWig coverage track from one sorted, indexed BAM file "
        "using `bamnado bam-coverage`. Output format is inferred from output_path: "
        ".bedgraph/.bdg writes bedGraph; .bw/.bigwig writes BigWig. "
        "For ATAC-seq, set fragment_counts=True and optionally use "
        "min_fragment_len/max_fragment_len to restrict fragment classes."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def bamnado_bam_coverage(
    bam_path: str,
    output_path: str,
    bin_size: Optional[int] = 50,
    normalize: str = "raw",
    scale_factor: Optional[float] = None,
    fragment_counts: bool = False,
    shift: str = "0,0,0,0",
    truncate: Optional[str] = None,
    ignore_scaffolds: bool = False,
    threads: int = 6,
    strand: str = "both",
    proper_pairs: bool = False,
    min_mapq: int = 20,
    min_length: int = 20,
    max_length: int = 1000,
    min_fragment_len: Optional[int] = None,
    max_fragment_len: Optional[int] = None,
    blacklist: Optional[str] = None,
    barcode_allowlist: Optional[str] = None,
    read_group: Optional[str] = None,
    tag: Optional[str] = None,
    tag_value: Optional[str] = None,
) -> str:
    """
    Args:
        bam_path: Path to a sorted, indexed BAM file.
        output_path: Destination bedGraph or BigWig path. Extension controls format.
        bin_size: Coverage bin width in bp. Use 10-50 for fine signal, 100+ for quick tests.
        normalize: Signal normalisation: raw, cpm, or rpkm.
        scale_factor: Optional multiplier applied to final signal values.
        fragment_counts: Count paired-end fragment spans instead of individual reads.
        shift: Read/fragment end shifts as L,R,FL,FR. Default is no shift.
        truncate: Trim read/fragment ends as L,R,FL,FR.
        ignore_scaffolds: Skip scaffold/unplaced chromosomes.
        threads: Threads for BigWig writing.
        strand: both, forward, or reverse.
        proper_pairs: Keep only properly paired reads.
        min_mapq: Minimum mapping quality.
        min_length: Minimum read length.
        max_length: Maximum read length.
        min_fragment_len: Minimum paired-end fragment length.
        max_fragment_len: Maximum paired-end fragment length.
        blacklist: Optional BED file of regions to exclude.
        barcode_allowlist: Optional file of barcodes to retain, one per line.
        read_group: Optional read group to retain.
        tag: Optional SAM tag to filter on.
        tag_value: Required tag value when tag is set.
    """
    _require_cli("bamnado", "bamnado", "bamnado")
    args = [
        "bamnado", "bam-coverage",
        "--bam", bam_path,
        "--output", output_path,
    ]
    _add_coverage_args(
        args,
        bin_size=bin_size,
        normalize=normalize,
        scale_factor=scale_factor,
        fragment_counts=fragment_counts,
        shift=shift,
        truncate=truncate,
        ignore_scaffolds=ignore_scaffolds,
        threads=threads,
    )
    _add_read_filter_args(
        args,
        strand=strand,
        proper_pairs=proper_pairs,
        min_mapq=min_mapq,
        min_length=min_length,
        max_length=max_length,
        min_fragment_len=min_fragment_len,
        max_fragment_len=max_fragment_len,
        blacklist=blacklist,
        barcode_allowlist=barcode_allowlist,
        read_group=read_group,
        tag=tag,
        tag_value=tag_value,
    )

    output = _run(args)
    return output or f"Coverage written to {output_path}"


@mcp.tool(
    description=(
        "Merge coverage from multiple sorted, indexed BAM files into one bedGraph "
        "or BigWig using `bamnado multi-bam-coverage`. Output format is inferred "
        "from output_path."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def bamnado_multi_bam_coverage(
    bam_paths: list[str],
    output_path: str,
    bin_size: Optional[int] = 50,
    normalize: str = "raw",
    scale_factor: Optional[float] = None,
    fragment_counts: bool = False,
    shift: str = "0,0,0,0",
    truncate: Optional[str] = None,
    ignore_scaffolds: bool = False,
    threads: int = 6,
    strand: str = "both",
    proper_pairs: bool = False,
    min_mapq: int = 20,
    min_length: int = 20,
    max_length: int = 1000,
    min_fragment_len: Optional[int] = None,
    max_fragment_len: Optional[int] = None,
    blacklist: Optional[str] = None,
    barcode_allowlist: Optional[str] = None,
    read_group: Optional[str] = None,
    tag: Optional[str] = None,
    tag_value: Optional[str] = None,
) -> str:
    """
    Args:
        bam_paths: Input BAM files to merge into one coverage track.
        output_path: Destination bedGraph or BigWig path.
        bin_size: Coverage bin width in bp.
        normalize: Signal normalisation: raw, cpm, or rpkm.
        scale_factor: Optional multiplier applied to final signal values.
        fragment_counts: Count paired-end fragment spans instead of individual reads.
        shift: Read/fragment end shifts as L,R,FL,FR.
        truncate: Trim read/fragment ends as L,R,FL,FR.
        ignore_scaffolds: Skip scaffold/unplaced chromosomes.
        threads: Threads for BigWig writing.
        strand: both, forward, or reverse.
        proper_pairs: Keep only properly paired reads.
        min_mapq: Minimum mapping quality.
        min_length: Minimum read length.
        max_length: Maximum read length.
        min_fragment_len: Minimum paired-end fragment length.
        max_fragment_len: Maximum paired-end fragment length.
        blacklist: Optional BED file of regions to exclude.
        barcode_allowlist: Optional file of barcodes to retain, one per line.
        read_group: Optional read group to retain.
        tag: Optional SAM tag to filter on.
        tag_value: Required tag value when tag is set.
    """
    _require_cli("bamnado", "bamnado", "bamnado")
    args = ["bamnado", "multi-bam-coverage", "--output", output_path]
    for bam in bam_paths:
        args += ["--bams", bam]
    _add_coverage_args(
        args,
        bin_size=bin_size,
        normalize=normalize,
        scale_factor=scale_factor,
        fragment_counts=fragment_counts,
        shift=shift,
        truncate=truncate,
        ignore_scaffolds=ignore_scaffolds,
        threads=threads,
    )
    _add_read_filter_args(
        args,
        strand=strand,
        proper_pairs=proper_pairs,
        min_mapq=min_mapq,
        min_length=min_length,
        max_length=max_length,
        min_fragment_len=min_fragment_len,
        max_fragment_len=max_fragment_len,
        blacklist=blacklist,
        barcode_allowlist=barcode_allowlist,
        read_group=read_group,
        tag=tag,
        tag_value=tag_value,
    )
    output = _run(args)
    return output or f"Coverage written to {output_path}"


@mcp.tool(
    description=(
        "Compare two BigWig files bin by bin using `bamnado bigwig-compare`. "
        "comparison must be subtraction, ratio, or log-ratio."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def bamnado_bigwig_compare(
    bw1: str,
    bw2: str,
    output_path: str,
    comparison: str,
    bin_size: int = 50,
    chunk_size: Optional[int] = None,
    pseudocount: Optional[float] = None,
    scale_factor_bw1: Optional[float] = None,
    scale_factor_bw2: Optional[float] = None,
    threads: int = 6,
) -> str:
    """
    Args:
        bw1: First BigWig input.
        bw2: Second BigWig input.
        output_path: Output BigWig path.
        comparison: subtraction, ratio, or log-ratio.
        bin_size: Bin size for comparison.
        chunk_size: Optional processing chunk size.
        pseudocount: Value added before ratio/log-ratio comparisons.
        scale_factor_bw1: Scale factor applied to bw1 before comparison.
        scale_factor_bw2: Scale factor applied to bw2 before comparison.
        threads: Threads for BigWig writing.
    """
    _require_cli("bamnado", "bamnado", "bamnado")
    args = [
        "bamnado", "bigwig-compare",
        "--bw1", bw1,
        "--bw2", bw2,
        "--output", output_path,
        "--comparison", comparison,
        "--bin-size", str(bin_size),
        "--threads", str(threads),
    ]
    if chunk_size is not None:
        args += ["--chunk-size", str(chunk_size)]
    if pseudocount is not None:
        args += ["--pseudocount", str(pseudocount)]
    if scale_factor_bw1 is not None:
        args += ["--scale-factor-bw1", str(scale_factor_bw1)]
    if scale_factor_bw2 is not None:
        args += ["--scale-factor-bw2", str(scale_factor_bw2)]
    output = _run(args)
    return output or f"BigWig comparison written to {output_path}"


@mcp.tool(
    description=(
        "Aggregate multiple BigWig files into one BigWig using "
        "`bamnado bigwig-aggregate`. method must be sum, mean, median, max, or min."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def bamnado_bigwig_aggregate(
    bigwig_paths: list[str],
    output_path: str,
    method: str,
    bin_size: int = 50,
    pseudocount: Optional[float] = None,
    scale_factors: Optional[list[float]] = None,
    threads: int = 6,
) -> str:
    """
    Args:
        bigwig_paths: BigWig files to aggregate.
        output_path: Output BigWig path.
        method: sum, mean, median, max, or min.
        bin_size: Bin size for aggregation.
        pseudocount: Value added to all inputs before aggregation.
        scale_factors: Per-file scale factors, in the same order as bigwig_paths.
        threads: Threads for BigWig writing.
    """
    _require_cli("bamnado", "bamnado", "bamnado")
    args = [
        "bamnado", "bigwig-aggregate",
        "--output", output_path,
        "--method", method,
        "--bin-size", str(bin_size),
        "--threads", str(threads),
    ]
    for bigwig in bigwig_paths:
        args += ["--bigwigs", bigwig]
    if pseudocount is not None:
        args += ["--pseudocount", str(pseudocount)]
    if scale_factors is not None:
        args.append("--scale-factors")
        args += [str(scale_factor) for scale_factor in scale_factors]
    output = _run(args)
    return output or f"BigWig aggregate written to {output_path}"


@mcp.tool(
    description=(
        "Collapse adjacent equal-score bins in a bedGraph using "
        "`bamnado collapse-bedgraph`. Provide input_path/output_path for file IO."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def bamnado_collapse_bedgraph(
    input_path: Optional[str] = None,
    output_path: Optional[str] = None,
) -> str:
    """
    Args:
        input_path: Input bedGraph path. If omitted, bamnado reads stdin.
        output_path: Output bedGraph path. If omitted, bamnado writes stdout.
    """
    _require_cli("bamnado", "bamnado", "bamnado")
    args = ["bamnado", "collapse-bedgraph"]
    if input_path:
        args += ["--input", input_path]
    if output_path:
        args += ["--output", output_path]
    output = _run(args)
    return output or (f"Collapsed bedGraph written to {output_path}" if output_path else "")


@mcp.tool(
    description=(
        "Split a BAM file using BamNado read filters. Writes filtered reads to "
        "files named from output_prefix."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def bamnado_split(
    input_path: str,
    output_prefix: str,
    strand: str = "both",
    proper_pairs: bool = False,
    min_mapq: int = 20,
    min_length: int = 20,
    max_length: int = 1000,
    min_fragment_len: Optional[int] = None,
    max_fragment_len: Optional[int] = None,
    blacklist: Optional[str] = None,
    barcode_allowlist: Optional[str] = None,
    read_group: Optional[str] = None,
    tag: Optional[str] = None,
    tag_value: Optional[str] = None,
) -> str:
    """
    Args:
        input_path: Input BAM file.
        output_prefix: Output prefix.
        strand: both, forward, or reverse.
        proper_pairs: Keep only properly paired reads.
        min_mapq: Minimum mapping quality.
        min_length: Minimum read length.
        max_length: Maximum read length.
        min_fragment_len: Minimum paired-end fragment length.
        max_fragment_len: Maximum paired-end fragment length.
        blacklist: Optional BED file of regions to exclude.
        barcode_allowlist: Optional file of barcodes to retain, one per line.
        read_group: Optional read group to retain.
        tag: Optional SAM tag to filter on.
        tag_value: Required tag value when tag is set.
    """
    _require_cli("bamnado", "bamnado", "bamnado")
    args = ["bamnado", "split", "--input", input_path, "--output", output_prefix]
    _add_read_filter_args(
        args,
        strand=strand,
        proper_pairs=proper_pairs,
        min_mapq=min_mapq,
        min_length=min_length,
        max_length=max_length,
        min_fragment_len=min_fragment_len,
        max_fragment_len=max_fragment_len,
        blacklist=blacklist,
        barcode_allowlist=barcode_allowlist,
        read_group=read_group,
        tag=tag,
        tag_value=tag_value,
    )
    output = _run(args)
    return output or f"Split BAM written with prefix {output_prefix}"


@mcp.tool(
    description=(
        "Split a BAM into endogenous and exogenous reads using "
        "`bamnado split-exogenous`. Exogenous reads are identified by reference "
        "name prefix."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def bamnado_split_exogenous(
    input_path: str,
    output_prefix: str,
    exogenous_prefix: str,
    stats_path: Optional[str] = None,
    allow_unknown_mapq: bool = False,
    strand: str = "both",
    proper_pairs: bool = False,
    min_mapq: int = 20,
    min_length: int = 20,
    max_length: int = 1000,
    min_fragment_len: Optional[int] = None,
    max_fragment_len: Optional[int] = None,
    blacklist: Optional[str] = None,
    barcode_allowlist: Optional[str] = None,
    read_group: Optional[str] = None,
    tag: Optional[str] = None,
    tag_value: Optional[str] = None,
) -> str:
    """
    Args:
        input_path: Input BAM file.
        output_prefix: Output prefix.
        exogenous_prefix: Reference-name prefix identifying exogenous contigs.
        stats_path: Optional summary statistics output path.
        allow_unknown_mapq: Allow reads with MAPQ 255, common in STAR output.
        strand: both, forward, or reverse.
        proper_pairs: Keep only properly paired reads.
        min_mapq: Minimum mapping quality.
        min_length: Minimum read length.
        max_length: Maximum read length.
        min_fragment_len: Minimum paired-end fragment length.
        max_fragment_len: Maximum paired-end fragment length.
        blacklist: Optional BED file of regions to exclude.
        barcode_allowlist: Optional file of barcodes to retain, one per line.
        read_group: Optional read group to retain.
        tag: Optional SAM tag to filter on.
        tag_value: Required tag value when tag is set.
    """
    _require_cli("bamnado", "bamnado", "bamnado")
    args = [
        "bamnado", "split-exogenous",
        "--input", input_path,
        "--output", output_prefix,
        "--exogenous-prefix", exogenous_prefix,
    ]
    if stats_path:
        args += ["--stats", stats_path]
    if allow_unknown_mapq:
        args.append("--allow-unknown-mapq")
    _add_read_filter_args(
        args,
        strand=strand,
        proper_pairs=proper_pairs,
        min_mapq=min_mapq,
        min_length=min_length,
        max_length=max_length,
        min_fragment_len=min_fragment_len,
        max_fragment_len=max_fragment_len,
        blacklist=blacklist,
        barcode_allowlist=barcode_allowlist,
        read_group=read_group,
        tag=tag,
        tag_value=tag_value,
    )
    output = _run(args)
    return output or f"Endogenous/exogenous BAMs written with prefix {output_prefix}"


@mcp.tool(
    description=(
        "Filter and/or adjust reads in a BAM file using `bamnado modify`. "
        "Set tn5_shift=True to apply the standard Tn5 offset."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def bamnado_modify(
    input_path: str,
    output_prefix: str,
    tn5_shift: bool = False,
    strand: str = "both",
    proper_pairs: bool = False,
    min_mapq: int = 20,
    min_length: int = 20,
    max_length: int = 1000,
    min_fragment_len: Optional[int] = None,
    max_fragment_len: Optional[int] = None,
    blacklist: Optional[str] = None,
    barcode_allowlist: Optional[str] = None,
    read_group: Optional[str] = None,
    tag: Optional[str] = None,
    tag_value: Optional[str] = None,
) -> str:
    """
    Args:
        input_path: Input BAM file.
        output_prefix: Output prefix.
        tn5_shift: Apply the standard Tn5 offset.
        strand: both, forward, or reverse.
        proper_pairs: Keep only properly paired reads.
        min_mapq: Minimum mapping quality.
        min_length: Minimum read length.
        max_length: Maximum read length.
        min_fragment_len: Minimum paired-end fragment length.
        max_fragment_len: Maximum paired-end fragment length.
        blacklist: Optional BED file of regions to exclude.
        barcode_allowlist: Optional file of barcodes to retain, one per line.
        read_group: Optional read group to retain.
        tag: Optional SAM tag to filter on.
        tag_value: Required tag value when tag is set.
    """
    _require_cli("bamnado", "bamnado", "bamnado")
    args = ["bamnado", "modify", "--input", input_path, "--output", output_prefix]
    _add_read_filter_args(
        args,
        strand=strand,
        proper_pairs=proper_pairs,
        min_mapq=min_mapq,
        min_length=min_length,
        max_length=max_length,
        min_fragment_len=min_fragment_len,
        max_fragment_len=max_fragment_len,
        blacklist=blacklist,
        barcode_allowlist=barcode_allowlist,
        read_group=read_group,
        tag=tag,
        tag_value=tag_value,
    )
    if tn5_shift:
        args.append("--tn5-shift")
    output = _run(args)
    return output or f"Modified BAM written with prefix {output_prefix}"


@mcp.tool(
    description=(
        "Infer scaling factor and library size from a CPM/RPKM-normalised BigWig "
        "using `bamnado bigwig-infer-scale`."
    ),
    annotations={"readOnlyHint": True, "destructiveHint": False},
)
def bamnado_bigwig_infer_scale(
    bigwig_path: str,
    format: str = "table",
) -> str:
    """
    Args:
        bigwig_path: Input CPM/RPKM-normalised BigWig file.
        format: Output format: table or tsv.
    """
    _require_cli("bamnado", "bamnado", "bamnado")
    return _run([
        "bamnado", "bigwig-infer-scale",
        "--bigwig", bigwig_path,
        "--format", format,
    ])


# ─── MCCNado ─────────────────────────────────────────────────────────────────
#
# MCC pipeline order:
#   Step 1  mccnado_deduplicate_fastq      — PCR dedup on raw FASTQs (pre-alignment)
#   [external]  align with BWA-MEM2
#   Step 2  mccnado_annotate_bam           — add VP/OC/RT tags (BAM sorted by name)
#   Step 3  mccnado_deduplicate_bam        — coordinate dedup using annotation tags
#   Step 4  mccnado_split_viewpoint_reads  — split into per-viewpoint BAMs
#   Step 5  mccnado_identify_ligation_junctions — per-viewpoint cooler files
#   Step 6  mccnado_extract_ligation_stats — JSON cis/trans stats (parallel with 5)
#   Step 7  mccnado_combine_coolers        — merge per-viewpoint coolers

@mcp.tool(
    description=(
        "MCC Step 1 of 7: Deduplicate paired-end FASTQ files before alignment. "
        "Removes PCR duplicates from raw reads based on sequence identity. "
        "Run before alignment with BWA-MEM2. "
        "Next step: align deduplicated FASTQs, then call mccnado_annotate_bam."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def mccnado_deduplicate_fastq(
    fastq_1: str,
    fastq_2: str,
    output_1: str,
    output_2: str,
    stats_prefix: str,
) -> str:
    """
    Args:
        fastq_1: Path to R1 FASTQ file (gzip-compressed OK).
        fastq_2: Path to R2 FASTQ file (gzip-compressed OK).
        output_1: Output path for deduplicated R1 FASTQ.
        output_2: Output path for deduplicated R2 FASTQ.
        stats_prefix: Filename prefix for statistics output files
            (e.g. 'dedup_stats/sample1' writes 'sample1.json' and 'sample1.txt').
    """
    return _run([
        "mccnado", "deduplicate-fastq",
        fastq_1, fastq_2, output_1, output_2,
        "--stats-prefix", stats_prefix,
    ])


@mcp.tool(
    description=(
        "MCC Step 2 of 7: Annotate an aligned MCC BAM with viewpoint (VP), "
        "on-capture (OC), and reporter (RT) tags needed for all downstream steps. "
        "BAM must be sorted by query name (samtools sort -n) before calling this. "
        "Next step: mccnado_deduplicate_bam."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def mccnado_annotate_bam(bam: str, output: str) -> str:
    """
    Args:
        bam: Path to aligned BAM sorted by query name (samtools sort -n).
            Coordinate-sorted BAMs will fail — sort by name first.
        output: Path for the annotated output BAM.
    """
    return _run(["mccnado", "annotate-bam", bam, output])


@mcp.tool(
    description=(
        "MCC Step 3 of 7: Remove PCR duplicates from an annotated MCC BAM. "
        "Uses coordinate positions together with VP/OC/RT annotation tags "
        "for accurate duplicate detection in MCC data. "
        "Input must be the output of mccnado_annotate_bam. "
        "Next step: mccnado_split_viewpoint_reads."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def mccnado_deduplicate_bam(bam: str, output: str) -> str:
    """
    Args:
        bam: Annotated BAM from mccnado_annotate_bam (coordinate-sorted is fine here).
        output: Path for the deduplicated output BAM.
    """
    return _run(["mccnado", "deduplicate-bam", bam, output])


@mcp.tool(
    description=(
        "MCC Step 4 of 7: Split a deduplicated MCC BAM by viewpoint into separate BAM files. "
        "Each output BAM contains reads from a single viewpoint for independent analysis. "
        "Input must be the output of mccnado_deduplicate_bam. "
        "Next step: mccnado_identify_ligation_junctions on each viewpoint BAM."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def mccnado_split_viewpoint_reads(bam: str, output: str) -> str:
    """
    Args:
        bam: Deduplicated, annotated BAM from mccnado_deduplicate_bam.
        output: Output directory where per-viewpoint BAM files will be written.
    """
    return _run(["mccnado", "split-viewpoint-reads", bam, output])


@mcp.tool(
    description=(
        "MCC Step 5 of 7: Identify ligation junctions in an MCC BAM and write "
        "per-viewpoint cooler files to the output directory. "
        "Can be run in parallel with mccnado_extract_ligation_stats (Step 6). "
        "Next step: mccnado_combine_coolers to merge the per-viewpoint coolers."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def mccnado_identify_ligation_junctions(bam: str, outdir: str) -> str:
    """
    Args:
        bam: Annotated, deduplicated BAM (output of mccnado_deduplicate_bam
            or a per-viewpoint BAM from mccnado_split_viewpoint_reads).
        outdir: Directory where per-viewpoint .cool files will be written.
            One cooler file is created per viewpoint.
    """
    return _run(["mccnado", "identify-ligation-junctions", bam, outdir])


@mcp.tool(
    description=(
        "MCC Step 6 of 7 (parallel with Step 5): Compute ligation junction statistics "
        "from an MCC BAM and write a JSON/TSV report. "
        "Reports cis- vs. trans-interaction ratios and junction counts per viewpoint. "
        "Can be run concurrently with mccnado_identify_ligation_junctions."
    ),
    annotations={"readOnlyHint": True, "destructiveHint": False},
)
def mccnado_extract_ligation_stats(bam: str, stats_output: str) -> str:
    """
    Args:
        bam: Annotated, deduplicated BAM (output of mccnado_deduplicate_bam).
        stats_output: Path for the output stats file (JSON or TSV depending on extension).
    """
    return _run(["mccnado", "extract-ligation-stats", bam, stats_output])


@mcp.tool(
    description=(
        "MCC Step 7 of 7: Merge multiple per-viewpoint ligation junction cooler files "
        "into a single multi-viewpoint cooler. "
        "Input cooler files are produced by mccnado_identify_ligation_junctions."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def mccnado_combine_coolers(cooler_files: list[str], output: str) -> str:
    """
    Args:
        cooler_files: List of per-viewpoint .cool files from
            mccnado_identify_ligation_junctions. All must share the same
            resolution and reference genome.
        output: Output path for the merged multi-viewpoint cooler file.
    """
    return _run(
        ["mccnado", "combine-ligation-junction-coolers"] + cooler_files + ["--output", output]
    )


# ─── PlotNado ────────────────────────────────────────────────────────────────

_PLOTNADO_YAML_SCHEMA = """\
YAML template schema:
  genome: hg38          # genome build; enables gene guides (hg38, mm10, dm6, ...)
  width: 12.0           # figure width in inches
  track_height: 1.0     # default panel height multiplier
  guides:
    genes: true         # add gene annotation panel (requires genome: set)
  tracks:               # ordered list of panels, rendered top-to-bottom
    - path: sample.bw
      type: bigwig      # track type -- see valid types below
      title: My Signal
      color: "#1f77b4"  # any matplotlib colour string or name
      height: 1.5       # relative panel height (multiplies track_height)
      style: fill       # fill | fragment | scatter | std  (bigwig only)
      group: grp1       # shared autoscale/colour group reference
      options: {}       # extra kwargs forwarded to the figure method
  groups:               # optional shared-scale group definitions
    - name: grp1
      autoscale: true
      autocolor: true

Valid track types (aliases in parentheses):
  bigwig (bw, signal, bedgraph)  -- continuous signal from BigWig file
  bed (annotation)               -- BED/BigBed intervals or peaks
  narrowpeak                     -- ENCODE narrowPeak with score and summit display
  genes (gene)                   -- gene models (requires genome: key to be set)
  links                          -- paired anchors / loops from BEDPE-like file
  overlay (bigwig_overlay)       -- multiple BigWig signals on one shared y-axis
  scalebar (scale)               -- scale bar panel
  axis                           -- y-axis reference panel
  spacer                         -- blank gap between panels\
"""

_PLOTNADO_REGION_NOTE = (
    "region: genomic window as 'chrN:start-end', "
    "e.g. 'chr1:1,000,000-1,100,000' (commas optional)."
)

_PLOTNADO_OUTPUT_NOTE = (
    "output_path: file extension sets format -- .png, .svg, or .pdf."
)

@mcp.tool(
    description=(
        "Render a genomic figure from a plotnado YAML template file. "
        "Generate a starter template with plotnado_init_template, then edit it, "
        "then render it here. "
        + _PLOTNADO_REGION_NOTE + " "
        + _PLOTNADO_OUTPUT_NOTE + "\n\n"
        + _PLOTNADO_YAML_SCHEMA
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def plotnado_from_template(
    template_yaml: str,
    region: str,
    output_path: str,
    dpi: int = 300,
    theme: str = "publication",
) -> str:
    """
    Args:
        template_yaml: Path to a YAML config file matching the schema described
            above. Generate a starter with plotnado_init_template, then edit it
            to adjust colours, heights, and track order before rendering.
        region: Genomic window to render, e.g. 'chr1:1,000,000-1,100,000'.
            Commas in coordinates are optional.
        output_path: Destination file. Extension sets format: .png, .svg, or .pdf.
        dpi: Raster resolution in dots per inch. 300 for publication figures,
            150 for quick preview.
        theme: Visual theme. 'publication' gives clean, minimal axes suitable
            for papers. 'default' uses standard matplotlib styling.
    """
    try:
        from plotnado.figure import GenomicFigure
    except ImportError:
        raise RuntimeError(_install_hint("plotnado", "plotnado"))

    fig = GenomicFigure.from_template(template_yaml, theme=theme)
    fig.save(output_path, region=region, dpi=dpi)
    return f"Saved to {output_path}"


@mcp.tool(
    description=(
        "Render a genomic figure from a saved IGV session XML file. "
        "Uses the session's stored locus if region is not provided. "
        + _PLOTNADO_REGION_NOTE + " "
        + _PLOTNADO_OUTPUT_NOTE
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def plotnado_from_igv_session(
    session_xml: str,
    output_path: str,
    region: Optional[str] = None,
    dpi: int = 300,
    theme: str = "publication",
) -> str:
    """
    Args:
        session_xml: Path to a saved IGV session XML file (File -> Save Session
            in IGV). Track order, colours, and locus are read from the session.
        output_path: Destination file. Extension sets format: .png, .svg, or .pdf.
        region: Override the session's stored locus with a custom window,
            e.g. 'chr1:1,000,000-1,100,000'. If omitted, the session locus is used.
        dpi: Raster resolution. 300 for publication, 150 for preview.
        theme: 'publication' for clean axes; 'default' for standard matplotlib.
    """
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
        "Build and render a genomic figure from a list of track specifications -- "
        "no YAML file needed. Each element in 'tracks' is a dict with a 'type' key "
        "plus type-specific params. Track dict examples:\n"
        '  {"type": "bigwig",     "path": "s.bw",  "title": "Signal", "color": "#1f77b4", "style": "fill"}\n'
        '  {"type": "bed",        "path": "p.bed",  "title": "Peaks",  "show_labels": true}\n'
        '  {"type": "narrowpeak", "path": "p.narrowPeak", "color_by": "score"}\n'
        '  {"type": "genes",      "title": "Genes"}   (requires genome set in GenomicFigure)\n'
        '  {"type": "links",      "path": "loops.bedpe", "color_by_score": true}\n'
        '  {"type": "overlay",    "paths": ["a.bw","b.bw"], "titles": ["A","B"]}\n'
        '  {"type": "scalebar"}\n'
        '  {"type": "axis"}\n'
        '  {"type": "spacer"}\n\n'
        + _PLOTNADO_REGION_NOTE + " "
        + _PLOTNADO_OUTPUT_NOTE
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
    """
    Args:
        tracks: Ordered list of track specification dicts, rendered top-to-bottom.
            Each dict must have a 'type' key (see description for valid types and
            examples). Additional keys are type-specific parameters.
        region: Genomic window to render, e.g. 'chr1:1,000,000-1,100,000'.
        output_path: Destination file. Extension sets format: .png, .svg, or .pdf.
        dpi: Raster resolution. 300 for publication, 150 for preview.
        theme: 'publication' for clean axes; 'default' for standard matplotlib.
        width: Figure width in inches (default 12).
        highlight_regions: Optional list of 'chrN:start-end' windows to shade
            behind all tracks, e.g. to highlight a peak or regulatory element.
    """
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
    description=(
        "Infer a plotnado YAML template from a list of track files. "
        "Detects file types from extensions (.bw/.bigwig -> bigwig, .bed -> bed, "
        ".narrowPeak -> narrowpeak, .bedpe -> links) and writes a ready-to-edit "
        "template. Edit the template to adjust colours, heights, and order, "
        "then pass it to plotnado_from_template."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def plotnado_init_template(
    input_files: list[str],
    output_yaml: str,
    genome: str = "hg38",
) -> str:
    """
    Args:
        input_files: List of data file paths. Types are auto-detected from
            extensions: .bw/.bigwig -> bigwig, .bed -> bed,
            .narrowPeak -> narrowpeak, .bedpe -> links, .bigBed -> bed.
        output_yaml: Path to write the generated YAML template. Edit this file
            before passing it to plotnado_from_template.
        genome: Genome build name for gene annotation guides
            (e.g. 'hg38', 'mm10', 'dm6', 'hg19').
    """
    return _run(
        ["plotnado", "init", "--output", output_yaml, "--genome", genome] + input_files
    )


# ─── QuantNado ───────────────────────────────────────────────────────────────

@mcp.tool(
    description=(
        "Create a QuantNado Zarr v3 dataset from a list of BAM files. "
        "Quantifies binned coverage for all samples and writes a single Zarr store. "
        "Each BAM becomes one sample; .bai index must exist alongside each BAM. "
        "After creation, use quantnado_dataset_info to inspect sample names and "
        "chromosomes, then quantnado_extract_region or quantnado_call_peaks."
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
    """
    Args:
        bam_files: List of absolute paths to sorted, indexed BAM files.
            One sample per BAM. Each file must have a .bai index alongside it.
        output_zarr: Path for the new Zarr v3 store directory (will be created).
        chrom_sizes: Path to a 2-column TSV (chrom\\tsize). If None, chromosome
            sizes are read from the BAM headers. Provide this when BAMs have
            inconsistent headers across samples.
        bin_size: Coverage bin width in bp. 200 is standard for ChIP-seq and
            ATAC-seq. Use 10 for RNA-seq or when fine resolution is needed.
        n_workers: Number of parallel workers for BAM processing. Scale to the
            number of available CPU cores.

    Returns:
        Confirmation string with the output path and list of sample names ingested.
    """
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
    """
    Args:
        dataset_path: Path to an existing QuantNado Zarr store directory.

    Returns:
        Dict with 'samples' (list of sample names), 'chromosomes', 'modalities',
        'n_completed' (samples with complete data), and 'chromsizes' mapping.
    """
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
    description=(
        "Extract coverage signal for a genomic region across samples and save to CSV. "
        "Region format: 'chrN:start-end' (e.g. 'chr1:1000000-2000000'). "
        "Use quantnado_dataset_info to get valid sample names and chromosomes."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def quantnado_extract_region(
    dataset_path: str,
    region: str,
    output_csv: str,
    normalise: Optional[str] = None,
    samples: Optional[list[str]] = None,
) -> str:
    """
    Args:
        dataset_path: Path to an existing QuantNado Zarr store.
        region: Genomic window as 'chrN:start-end', e.g. 'chr1:1000000-2000000'.
        output_csv: Path for the output CSV file (bins x samples).
        normalise: Normalisation method. 'cpm' = counts per million (sequencing
            depth correction), 'rpkm' = reads per kilobase per million (depth +
            feature length correction), None = raw bin counts.
        samples: List of sample names to include (from quantnado_dataset_info).
            None returns all samples in the dataset.
    """
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
    description=(
        "Call peaks from a QuantNado dataset and write results to BED. "
        "Choose method by assay type: "
        "'quantile' -- fast threshold, sharp ChIP/ATAC peaks; "
        "'seacr' -- Sparse Enrichment Analysis, designed for CUT&RUN/CUT&TAG low-background data; "
        "'lanceotron' -- ML caller, handles broad and irregular peaks; "
        "'macs3' -- classical ChIP-seq narrow peak calling; "
        "'macs3_broad' -- broad domains (H3K27me3, H3K9me3); "
        "'unet' / 'resnet' -- deep-learning models (require GPU extras)."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def quantnado_call_peaks(
    dataset_path: str,
    output_bed: str,
    method: str = "quantile",
    sample_name: Optional[str] = None,
    blacklist_file: Optional[str] = None,
) -> str:
    """
    Args:
        dataset_path: Path to an existing QuantNado Zarr store.
        output_bed: Path for the output BED file of called peaks.
        method: Peak calling algorithm. Choose based on assay: 'quantile' or
            'lanceotron' for ChIP/ATAC; 'seacr' for CUT&RUN/CUT&TAG;
            'macs3' or 'macs3_broad' for classical ChIP; 'unet'/'resnet'
            for deep-learning-based calling (GPU required).
        sample_name: Name of the sample to call peaks for (from
            quantnado_dataset_info). None calls peaks across all samples.
        blacklist_file: Path to a BED file of regions to exclude from peak
            calling (e.g. ENCODE blacklist for hg38). Optional but recommended.
    """
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
    return f"Called peaks with {method}: {len(peaks)} peaks -> {output_bed}"


@mcp.tool(
    description=(
        "Generate a metaplot (average signal over genomic intervals) and save to PNG. "
        "Useful for visualising signal enrichment at peaks, promoters, or other features."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def quantnado_metaplot(
    dataset_path: str,
    intervals_path: str,
    output_png: str,
    normalise: Optional[str] = "cpm",
    samples: Optional[list[str]] = None,
) -> str:
    """
    Args:
        dataset_path: Path to an existing QuantNado Zarr store.
        intervals_path: Path to a BED file of genomic intervals to average over
            (e.g. peak calls or TSS regions).
        output_png: Output path for the metaplot PNG image.
        normalise: 'cpm', 'rpkm', or None. Defaults to 'cpm'.
        samples: Sample names to include; None = all samples.
    """
    try:
        from quantnado import QuantNado
    except ImportError:
        raise RuntimeError(_install_hint("quantnado", "quantnado"))

    qn = QuantNado.open_dataset(dataset_path)
    fig = qn.metaplot(intervals_path=intervals_path, normalise=normalise, samples=samples)
    fig.savefig(output_png, dpi=150, bbox_inches="tight")
    return f"Metaplot saved to {output_png}"


@mcp.tool(
    description=(
        "Run PCA on genome-wide sample coverage and save a scatter plot to PNG. "
        "Useful for quality control and sample clustering."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def quantnado_pca(
    dataset_path: str,
    output_png: str,
    normalise: Optional[str] = "cpm",
    samples: Optional[list[str]] = None,
) -> str:
    """
    Args:
        dataset_path: Path to an existing QuantNado Zarr store.
        output_png: Output path for the PCA scatter plot PNG.
        normalise: 'cpm', 'rpkm', or None. Defaults to 'cpm'.
        samples: Sample names to include; None = all samples.
    """
    try:
        from quantnado import QuantNado
    except ImportError:
        raise RuntimeError(_install_hint("quantnado", "quantnado"))

    qn = QuantNado.open_dataset(dataset_path)
    fig = qn.pca(normalise=normalise, samples=samples)
    fig.savefig(output_png, dpi=150, bbox_inches="tight")
    return f"PCA plot saved to {output_png}"


@mcp.tool(
    description=(
        "Plot coverage signal at a genomic locus for all (or selected) samples "
        "and save to PNG. Useful for inspecting individual loci across samples."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def quantnado_locus_plot(
    dataset_path: str,
    region: str,
    output_png: str,
    normalise: Optional[str] = "cpm",
    samples: Optional[list[str]] = None,
) -> str:
    """
    Args:
        dataset_path: Path to an existing QuantNado Zarr store.
        region: Genomic window as 'chrN:start-end', e.g. 'chr1:1000000-1010000'.
        output_png: Output path for the locus plot PNG.
        normalise: 'cpm', 'rpkm', or None. Defaults to 'cpm'.
        samples: Sample names to include; None = all samples.
    """
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
        "Generate a UCSC track hub directory from BigWig/BigBed/BED track files. "
        "Two input modes: (a) pass input_files directly for a flat hub; "
        "(b) pass metadata_csv to drive grouping, colouring, and hierarchy. "
        "The hub must be served over HTTP -- set url_prefix to the public base URL "
        "of output_dir so UCSC can fetch the files. "
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
    """
    Args:
        output_dir: Local directory where hub files will be written. Contains
            hub.txt, genomes.txt, and trackDb.txt files after generation.
        input_files: List of BigWig, BigBed, or BED file paths. Used when no
            metadata_csv is provided. Track types are auto-detected from extensions.
        metadata_csv: CSV file with a column matching track filenames/paths, plus
            additional metadata columns (e.g. 'condition', 'antibody', 'cell_type').
            These columns drive color_by, supergroup_by, and subgroup_by grouping.
        genome: UCSC genome assembly name shown in the hub (e.g. 'hg38', 'mm10',
            'dm6', 'hg19').
        hub_name: Display name for the hub as shown in the UCSC browser track list.
        hub_email: Contact email address included in hub.txt for UCSC registration.
        color_by: Column name in metadata_csv whose values determine track colours.
            Each unique value gets a distinct colour from a palette.
        supergroup_by: List of metadata column names used to create a SuperTrack
            hierarchy (top-level collapsible groups in UCSC).
        subgroup_by: List of metadata column names used to create a CompositeTrack
            matrix with filter controls in UCSC.
        url_prefix: Base HTTP URL where output_dir will be publicly accessible.
            UCSC fetches hub files via HTTP -- local paths will not work.
            Example: 'https://userweb.molbiol.ox.ac.uk/public/username'.
        seqnado_layout: If True, use SeqNado-specific seqnado_output/ directory
            conventions when locating track files.
    """
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
#
# Canonical SeqNado workflow:
#   (optional) seqnado_download           -- fetch public FASTQ from GEO/SRA
#   (optional) seqnado_build_genome       -- build and register a reference genome
#   1. seqnado_list_assays                -- confirm assay name and available methods
#   2. seqnado_list_genomes               -- confirm genome name is registered
#   3. seqnado_generate_design            -- parse FASTQs -> metadata CSV
#      [edit metadata CSV: fix sample_id and ip field names]
#   4. seqnado_generate_config            -- generate config YAML
#      [edit config YAML to tune pipeline settings]
#   5. seqnado_run_pipeline dry_run=True  -- validate job graph
#   6. seqnado_run_pipeline dry_run=False -- execute pipeline
#   7. seqnado_pipeline_status            -- check completion

@mcp.tool(
    description=(
        "List all supported SeqNado assay types and available pipeline methods. "
        "Use this first to confirm the correct assay name and supported methods "
        "before calling seqnado_generate_design or seqnado_generate_config. "
        "Canonical workflow: seqnado_list_assays -> seqnado_list_genomes -> "
        "seqnado_generate_design -> seqnado_generate_config -> "
        "seqnado_run_pipeline(dry_run=True) -> seqnado_run_pipeline(dry_run=False)."
    ),
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
    description=(
        "List all genome configurations registered in SeqNado for a given assay. "
        "If the needed genome is not listed, run seqnado_build_genome to download "
        "and register it. Registered genomes are stored at "
        "~/.config/seqnado/genome_config.json."
    ),
    annotations={"readOnlyHint": True, "destructiveHint": False},
)
def seqnado_list_genomes(assay: str = "chip") -> dict[str, Any]:
    """
    Args:
        assay: Assay type to list genomes for. One of: chip, atac, rna, cat,
            meth, snp, mcc, crispr. Different assays may require different genome
            configurations (e.g. RNA-seq needs a STAR index; ChIP needs Bowtie2).
    """
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
        "SeqNado Step 1: Generate a design CSV from FASTQ file paths. "
        "Parses Illumina-style filenames to extract sample names, replicates, "
        "and read pairs. "
        "IMPORTANT: After generation, edit the CSV directly to rename outputs -- "
        "changing 'sample_id' renames sample directories in all pipeline outputs "
        "without touching the original FASTQ files. For ChIP/CUT&RUN/TAG assays, "
        "changing the 'ip' column renames the antibody label in BigWig and peak "
        "output filenames ({sample_id}_{ip}.*). This is the correct way to clean "
        "up messy sequencer filenames. "
        "Output paths in the pipeline follow seqnado_output/{assay}/{sample_id}/{step}/."
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
    """
    Args:
        assay: Assay type. One of: chip, atac, rna, cat, meth, snp, mcc, crispr.
            Use seqnado_list_assays to see all options.
        fastq_files: List of FASTQ file paths. Illumina naming convention
            '{sample}_{antibody}_R{1,2}_001.fastq.gz' is parsed automatically
            to extract sample names and read pairing.
        output_csv: Output path for the design CSV (convention: metadata_{assay}.csv
            in the project working directory).
        ip_to_control: Explicit antibody-to-control mapping as comma-separated
            'IP1:ctrl1,IP2:ctrl2' pairs (e.g. 'H3K27ac:Input,H3K4me3:Input').
            Use when multiple controls exist and auto-detection is ambiguous.
        group_by: Metadata column name or regex pattern for grouping samples
            into comparison groups.
        deseq2_pattern: Regex to extract DESeq2 group label from filenames
            (RNA-seq only). Example: '-(WT|KO)-' extracts 'WT' or 'KO'.
            Common terms (control, treated, WT, KO, DMSO, vehicle) are
            auto-detected without this argument.

    Design CSV columns produced:
        sample_id     -- unique sample name; EDIT THIS to rename all pipeline outputs
        fastq_1       -- R1 FASTQ path (absolute)
        fastq_2       -- R2 FASTQ path (absolute)
        ip            -- antibody / IP target (ChIP/CAT/MCC); EDIT to rename in outputs
        control       -- control sample name (IP assays only)
        scaling_group -- normalisation group (default: 'default')
        condition     -- biological condition (optional)
        group         -- DESeq2 comparison group (RNA-seq)
        deseq2        -- 0=reference, 1=treatment (RNA-seq DESeq2)
    """
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
        "SeqNado Step 2: Generate a workflow config YAML non-interactively. "
        "This wraps the documented CLI path: "
        "`seqnado config ASSAY --no-interactive --no-make-dirs --render-options -o OUTPUT`. "
        "Use this instead of plain `seqnado config`, because MCP agents cannot "
        "answer terminal prompts. The generated YAML uses SeqNado defaults and "
        "the first registered compatible genome; review and edit the YAML before "
        "running the pipeline."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def seqnado_generate_config(
    assay: str,
    output_yaml: str,
) -> str:
    """
    Args:
        assay: Assay type. One of: chip, atac, rna, cat, meth, snp, mcc, crispr.
        output_yaml: Output path for the config YAML. Convention: config_{assay}.yaml
            in the project working directory.
    """
    try:
        from seqnado import Assay  # noqa: F401
    except ImportError:
        raise RuntimeError(_install_hint("seqnado", "seqnado"))

    assay_clean = assay.lower()
    if assay_clean in {"mcc", "multiomics"}:
        raise ValueError(
            f"SeqNado non-interactive config generation does not support '{assay}'. "
            "Run the interactive CLI outside MCP or create/edit the YAML manually."
        )

    args = [
        "seqnado", "config", assay_clean,
        "--no-interactive",
        "--no-make-dirs",
        "--render-options",
        "-o", output_yaml,
    ]
    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"seqnado config failed (exit {result.returncode}):\n"
            f"STDOUT:\n{result.stdout[-3000:]}\n"
            f"STDERR:\n{result.stderr[-3000:]}"
        )

    output_path = Path(output_yaml)
    if not output_path.exists():
        raise RuntimeError(
            "seqnado config exited successfully but did not create "
            f"the expected file: {output_yaml}"
        )

    messages = [f"Config written to {output_yaml} using non-interactive SeqNado defaults."]
    if result.stdout.strip():
        messages.append(f"stdout:\n{result.stdout.strip()[-1000:]}")
    if result.stderr.strip():
        messages.append(f"stderr:\n{result.stderr.strip()[-1000:]}")
    return "\n\n".join(messages)


@mcp.tool(
    description=(
        "SeqNado Step 3: Run (or dry-run) a SeqNado Snakemake pipeline. "
        "ALWAYS call with dry_run=True first to preview the job graph and validate "
        "inputs before committing to a full run. "
        "The pipeline reads metadata_{assay}.csv and config_{assay}.yaml from "
        "working_dir and writes all outputs to seqnado_output/ there. "
        "Execution profiles: 'le' (local+conda envs, default), "
        "'lc' (local+containers), 'ls' (local+Singularity), "
        "'ss' (SLURM+Singularity for HPC), 't' (testing). "
        "Pass ['--unlock'] in extra_snakemake_args if the pipeline is locked "
        "from a previous interrupted run."
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
    """
    Args:
        assay: Assay type. One of: chip, atac, rna, cat, meth, snp, mcc, crispr.
        config_yaml: Path to the YAML generated by seqnado_generate_config.
            Review and edit it before passing here.
        working_dir: Project directory that contains the design CSV
            (metadata_{assay}.csv) and where seqnado_output/ will be created.
        cores: CPU cores for local execution. Ignored by SLURM profiles ('ss')
            which schedule jobs via the cluster.
        dry_run: If True (default), preview the Snakemake job graph without
            executing any steps. Always run with True first, then set False
            to execute after reviewing the plan.
        targets: List of specific Snakemake rule names or output file paths to
            run, limiting execution to those targets only. None runs all rules.
        profile: Snakemake execution profile. 'le' = local+conda (default),
            'lc' = local+containers, 'ls' = local+Singularity,
            'ss' = SLURM+Singularity (HPC), 't' = testing preset.
        extra_snakemake_args: Additional flags passed verbatim to Snakemake.
            Common examples: ['--unlock'] to release a locked directory,
            ['--rerun-incomplete'] to retry unfinished jobs,
            ['--forceall'] to re-run all steps regardless of completion status.
    """
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
    """
    Args:
        working_dir: Project directory where the pipeline was run. Inspects
            the seqnado_output/ subdirectory for completed output files.
    """
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
    """
    Args:
        design_csv: Path to the design CSV generated by seqnado_generate_design.
        assay: Assay type the design was created for (chip, atac, rna, etc.).
            Used to determine which columns are required (e.g. ip/control for
            IP-based assays like chip, cat).
    """
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


@mcp.tool(
    description=(
        "Download public FASTQ files from GEO/SRA using an ENA metadata TSV. "
        "Obtain the TSV from the ENA Browser: search by project accession "
        "(e.g. PRJNA1234567 or ERP123456) then click 'Download report' as TSV. "
        "Downloaded files are named '{library_name}-{sample_title}_R1.fastq.gz' "
        "(paired-end) or '{library_name}-{sample_title}.fastq.gz' (single-end). "
        "Optionally auto-generates a SeqNado design CSV after download. "
        "Requires sra-tools (prefetch + fasterq-dump) in PATH. "
        "Always use dry_run=True first to preview the download plan."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def seqnado_download(
    metadata_tsv: str,
    output_dir: str = "fastqs",
    assay: Optional[str] = None,
    design_output: Optional[str] = None,
    cores: int = 4,
    preset: str = "le",
    dry_run: bool = True,
) -> str:
    """
    Args:
        metadata_tsv: Path to an ENA file report TSV. Required columns:
            'run_accession' (e.g. SRR123456), 'sample_title' (sample name),
            'library_name' (e.g. GSM identifier), 'library_layout' (PAIRED or SINGLE).
            Download from ENA Browser: search project accession -> Download report -> TSV.
        output_dir: Directory where FASTQ files will be written. Created if it
            does not exist. Defaults to 'fastqs'.
        assay: If provided, auto-generates a SeqNado design CSV after download.
            One of: chip, atac, rna, cat, meth, snp, mcc, crispr.
        design_output: Output path for the auto-generated design CSV.
            Defaults to 'metadata_{assay}.csv' in the current directory.
        cores: Number of parallel download jobs (runs prefetch + fasterq-dump).
        preset: Snakemake execution profile. 'le' = local+conda (default),
            'ls' = local+Singularity, 'ss' = SLURM+Singularity (HPC).
        dry_run: If True (default), preview the download plan without fetching
            any files. Review before setting False.
    """
    try:
        from seqnado import Assay  # noqa: F401
    except ImportError:
        raise RuntimeError(_install_hint("seqnado", "seqnado"))

    args = [
        "seqnado", "download",
        metadata_tsv,
        "--outdir", output_dir,
        "--cores", str(cores),
        "--preset", preset,
    ]
    if assay:
        args += ["--assay", assay]
    if design_output:
        args += ["--design-output", design_output]
    if dry_run:
        args.append("--dry-run")

    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"seqnado download failed (exit {result.returncode}):\n"
            f"STDOUT:\n{result.stdout[-3000:]}\n"
            f"STDERR:\n{result.stderr[-3000:]}"
        )
    return result.stdout.strip() or "Download completed."


@mcp.tool(
    description=(
        "Download a reference genome from UCSC and build all required indices: "
        "Bowtie2, STAR, samtools FAI, chromosome sizes, ENCODE blacklist, RefSeq GTF. "
        "On completion the genome is automatically registered at "
        "~/.config/seqnado/genome_config.json and immediately available via "
        "seqnado_list_genomes. Multiple genomes and spike-in composite builds are "
        "supported in a single call. STAR indexing requires at least 8 GB RAM. "
        "Always use dry_run=True first to preview the build DAG."
    ),
    annotations={"readOnlyHint": False, "destructiveHint": False},
)
def seqnado_build_genome(
    genome_names: list[str],
    output_dir: str,
    spikein: Optional[str] = None,
    cores: int = 4,
    preset: str = "le",
    dry_run: bool = True,
) -> str:
    """
    Args:
        genome_names: List of UCSC genome assembly identifiers to download and
            build. Examples: ['hg38'], ['hg38', 'mm39'], ['hg38', 'mm39', 'dm6'].
        output_dir: Directory where genome files will be written. The path is
            registered in ~/.config/seqnado/genome_config.json after a successful
            build so seqnado pipelines can find it automatically.
        spikein: Optional spike-in genome name to build a composite combined
            index alongside the primary genome (e.g. 'dm6' for Drosophila
            spike-in normalisation). Must be a valid UCSC genome name.
        cores: CPU cores for parallel index building. STAR genome generation
            requires at least 8 GB RAM; scale cores to available memory.
        preset: Snakemake execution profile. 'le' = local+conda (default),
            'ls' = local+Singularity, 'ss' = SLURM+Singularity (HPC).
        dry_run: If True (default), preview the build DAG without downloading
            or indexing anything. Review the plan before setting False.
    """
    try:
        from seqnado import Assay  # noqa: F401
    except ImportError:
        raise RuntimeError(_install_hint("seqnado", "seqnado"))

    args = [
        "seqnado", "genomes", "build",
        "--name", ",".join(genome_names),
        "--outdir", output_dir,
        "--cores", str(cores),
        "--preset", preset,
    ]
    if spikein:
        args += ["--spikein", spikein]
    if dry_run:
        args.append("--dry-run")

    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"seqnado genomes build failed (exit {result.returncode}):\n"
            f"STDOUT:\n{result.stdout[-3000:]}\n"
            f"STDERR:\n{result.stderr[-3000:]}"
        )
    return result.stdout.strip() or f"Genome build completed for {genome_names}."


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
