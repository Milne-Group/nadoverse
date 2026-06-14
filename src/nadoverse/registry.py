from __future__ import annotations

import json
import sys
from dataclasses import asdict, dataclass
from typing import Optional


@dataclass
class NadoTool:
    name: str
    description: str
    pypi_name: str
    install_extra: str
    cli_command: Optional[str]
    input_types: list[str]
    output_types: list[str]
    repo_url: str
    docs_url: Optional[str]
    min_python: str
    container_image: Optional[str] = None

    def to_dict(self) -> dict:
        return asdict(self)

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2)

    def is_available(self) -> bool:
        try:
            __import__(self.pypi_name.lower().replace("-", "_"))
            return True
        except ImportError:
            return False

    def installed_version(self) -> Optional[str]:
        try:
            from importlib.metadata import version
            return version(self.pypi_name)
        except Exception:
            return None

    def python_compatible(self) -> bool:
        """False if running Python version exceeds this tool's known upper cap."""
        caps = {
            "tabnado": (3, 13),  # tabnado requires <3.13
        }
        cap = caps.get(self.pypi_name.lower())
        if cap is None:
            return True
        return sys.version_info < cap


REGISTRY: list[NadoTool] = [
    NadoTool(
        name="SeqNado",
        description=(
            "Unified genomics pipelines: ATAC-seq, ChIP-seq, CUT&RUN/TAG, "
            "RNA-seq, WGS, Methylation (Bisulphite/TAPS), CRISPR screens, "
            "and Micro-Capture-C."
        ),
        pypi_name="seqnado",
        install_extra="seqnado",
        cli_command="seqnado",
        input_types=["FASTQ", "BAM"],
        output_types=["BAM", "bigWig", "bigBed", "peaks", "counts"],
        repo_url="https://github.com/Milne-Group/SeqNado",
        docs_url="https://Milne-Group.github.io/SeqNado",
        min_python="3.10",
    ),
    NadoTool(
        name="PlotNado",
        description=(
            "Genomic track visualisation: generate UCSC-style browser tracks "
            "and publication-quality plots from bigWig/bigBed files."
        ),
        pypi_name="plotnado",
        install_extra="plotnado",
        cli_command="plotnado",
        input_types=["bigWig", "bigBed", "BED", "config TOML"],
        output_types=["PNG", "SVG", "PDF"],
        repo_url="https://github.com/Milne-Group/PlotNado",
        docs_url=None,
        min_python="3.12",
        container_image="ghcr.io/milne-group/plotnado:latest",
    ),
    NadoTool(
        name="TrackNado",
        description=(
            "CLI utility to generate UCSC trackhubs from sequencing output files "
            "(bigWig, bigBed, etc.)."
        ),
        pypi_name="tracknado",
        install_extra="tracknado",
        cli_command="tracknado",
        input_types=["bigWig", "bigBed", "BED"],
        output_types=["UCSC trackhub"],
        repo_url="https://github.com/Milne-Group/TrackNado",
        docs_url=None,
        min_python="3.10",
    ),
    NadoTool(
        name="BamNado",
        description=(
            "High-performance BAM processing (Rust core, PyO3 bindings): "
            "parallel coverage/pileup, flexible read filtering, signal normalisation. "
            "Provides both a standalone CLI (`bamnado`) and Python bindings."
        ),
        pypi_name="bamnado",
        install_extra="bamnado",
        cli_command="bamnado",
        input_types=["BAM"],
        output_types=["bigWig", "coverage arrays"],
        repo_url="https://github.com/Milne-Group/BamNado",
        docs_url=None,
        min_python="3.10",
        container_image="ghcr.io/milne-group/bamnado:latest",
    ),
    NadoTool(
        name="QuantNado",
        description=(
            "Zarr-backed genomic signal storage and analysis: BAM/bigWig ingestion, "
            "signal reduction, feature counting, dimensionality reduction, and "
            "quantile-based peak calling."
        ),
        pypi_name="quantnado",
        install_extra="quantnado",
        cli_command="quantnado",
        input_types=["BAM", "bigWig"],
        output_types=["Zarr store", "peaks", "counts matrix"],
        repo_url="https://github.com/Milne-Group/QuantNado",
        docs_url=None,
        min_python="3.12",
        container_image="ghcr.io/milne-group/quantnado:latest",
    ),
    NadoTool(
        name="MCCNado",
        description=(
            "Rust-based tools for processing Micro-Capture-C data, used internally "
            "by SeqNado. Library only — no standalone CLI."
        ),
        pypi_name="mccnado",
        install_extra="mccnado",
        cli_command=None,
        input_types=["BAM", "pairs"],
        output_types=["cooler", "HDF5"],
        repo_url="https://github.com/alsmith151/MCCNado",
        docs_url=None,
        min_python="3.10",
        container_image="ghcr.io/alsmith151/mccnado:latest",
    ),
    NadoTool(
        name="TabNado",
        description=(
            "Predicts transcription-factor binding from epigenomic cofactors "
            "(ChIP-seq, CUT&TAG, CUT&RUN) using GANDALF or XGBoost. "
            "Requires Python <3.13."
        ),
        pypi_name="tabnado",
        install_extra="tabnado",
        cli_command="tabnado-run",
        input_types=["Zarr store", "BED", "config TOML"],
        output_types=["predictions TSV", "model checkpoint"],
        repo_url="https://github.com/Milne-Group/TabNado",
        docs_url=None,
        min_python="3.12",
        container_image="ghcr.io/milne-group/tabnado:latest",
    ),
]


def get_tool(name: str) -> Optional[NadoTool]:
    """Return first registry entry whose name or pypi_name matches (case-insensitive)."""
    key = name.lower()
    for tool in REGISTRY:
        if tool.name.lower() == key or tool.pypi_name.lower() == key:
            return tool
    return None


def all_tools() -> list[NadoTool]:
    return list(REGISTRY)


def to_json() -> str:
    return json.dumps([t.to_dict() for t in REGISTRY], indent=2)
