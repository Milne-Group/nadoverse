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
    monkeypatch.setattr(mcp_server.shutil, "which", lambda command: f"/bin/{command}")
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
            "--bin-size",
            "50",
            "--fragment-counts",
            "--strand",
            "both",
            "--min-mapq",
            "30",
            "--min-length",
            "20",
            "--max-length",
            "1000",
        ]
    ]
    assert result == f"Coverage written to {output}"


def test_bamnado_multi_bam_coverage_invokes_cli(monkeypatch, tmp_path):
    mcp_server = _import_mcp_server(monkeypatch)
    monkeypatch.setattr(mcp_server.shutil, "which", lambda command: f"/bin/{command}")
    calls = []

    def fake_run(args):
        calls.append(args)
        return ""

    monkeypatch.setattr(mcp_server, "_run", fake_run)
    output = tmp_path / "merged.bw"

    result = mcp_server.bamnado_multi_bam_coverage(
        ["/data/a.bam", "/data/b.bam"],
        str(output),
        normalize="cpm",
        proper_pairs=True,
    )

    assert calls == [
        [
            "bamnado",
            "multi-bam-coverage",
            "--output",
            str(output),
            "--bams",
            "/data/a.bam",
            "--bams",
            "/data/b.bam",
            "--normalize",
            "cpm",
            "--threads",
            "6",
            "--bin-size",
            "50",
            "--strand",
            "both",
            "--min-mapq",
            "20",
            "--min-length",
            "20",
            "--max-length",
            "1000",
            "--proper-pairs",
        ]
    ]
    assert result == f"Coverage written to {output}"


def test_bamnado_bigwig_compare_invokes_cli(monkeypatch, tmp_path):
    mcp_server = _import_mcp_server(monkeypatch)
    monkeypatch.setattr(mcp_server.shutil, "which", lambda command: f"/bin/{command}")
    calls = []

    def fake_run(args):
        calls.append(args)
        return ""

    monkeypatch.setattr(mcp_server, "_run", fake_run)
    output = tmp_path / "compare.bw"

    result = mcp_server.bamnado_bigwig_compare(
        "/data/a.bw",
        "/data/b.bw",
        str(output),
        "log-ratio",
        pseudocount=1.0,
    )

    assert calls == [
        [
            "bamnado",
            "bigwig-compare",
            "--bw1",
            "/data/a.bw",
            "--bw2",
            "/data/b.bw",
            "--output",
            str(output),
            "--comparison",
            "log-ratio",
            "--bin-size",
            "50",
            "--threads",
            "6",
            "--pseudocount",
            "1.0",
        ]
    ]
    assert result == f"BigWig comparison written to {output}"


def test_bamnado_modify_invokes_cli(monkeypatch, tmp_path):
    mcp_server = _import_mcp_server(monkeypatch)
    monkeypatch.setattr(mcp_server.shutil, "which", lambda command: f"/bin/{command}")
    calls = []

    def fake_run(args):
        calls.append(args)
        return ""

    monkeypatch.setattr(mcp_server, "_run", fake_run)
    output_prefix = tmp_path / "shifted"

    result = mcp_server.bamnado_modify(
        "/data/sample.bam",
        str(output_prefix),
        tn5_shift=True,
        tag="CB",
        tag_value="cell-1",
    )

    assert calls == [
        [
            "bamnado",
            "modify",
            "--input",
            "/data/sample.bam",
            "--output",
            str(output_prefix),
            "--strand",
            "both",
            "--min-mapq",
            "20",
            "--min-length",
            "20",
            "--max-length",
            "1000",
            "--tag",
            "CB",
            "--tag-value",
            "cell-1",
            "--tn5-shift",
        ]
    ]
    assert result == f"Modified BAM written with prefix {output_prefix}"


def test_regulonado_build_dataset_invokes_cli(monkeypatch, tmp_path):
    mcp_server = _import_mcp_server(monkeypatch)
    monkeypatch.setattr(mcp_server.shutil, "which", lambda command: f"/bin/{command}")
    calls = []

    def fake_run(args):
        calls.append(args)
        return ""

    monkeypatch.setattr(mcp_server, "_run", fake_run)
    output = tmp_path / "dataset"

    result = mcp_server.regulonado_build_dataset(
        "/data/intervals.bed",
        "/data/genome.fa",
        str(output),
        bigwig_files=["/data/a.bw", "/data/b.bw"],
        splits=["train:fold0,fold1", "validation:fold4"],
        overwrite=True,
        chroms=["chr1"],
    )

    assert calls == [
        [
            "regulonado",
            "build",
            "/data/intervals.bed",
            "/data/genome.fa",
            str(output),
            "--bigwig-glob",
            "*.bw",
            "--context-length",
            "524288",
            "--bin-size",
            "32",
            "--n-pred-bins",
            "6144",
            "--shift-max-bp",
            "0",
            "--io-threads",
            "8",
            "--num-proc",
            "1",
            "--writer-batch-size",
            "500",
            "--dedupe-tracks",
            "none",
            "--n-extract-threads",
            "32",
            "--signal-sample-chunk",
            "8",
            "--signal-track-chunk",
            "128",
            "--arrow-batch-size",
            "8",
            "--arrow-compression",
            "lz4",
            "--strategy",
            "chrom_pass",
            "--bigwig",
            "/data/a.bw",
            "--bigwig",
            "/data/b.bw",
            "--split",
            "train:fold0,fold1",
            "--split",
            "validation:fold4",
            "--overwrite",
            "--chrom",
            "chr1",
        ]
    ]
    assert result == f"ReguloNado dataset written to {output}"


def test_regulonado_build_dataset_requires_one_bigwig_source(monkeypatch, tmp_path):
    mcp_server = _import_mcp_server(monkeypatch)
    monkeypatch.setattr(mcp_server.shutil, "which", lambda command: f"/bin/{command}")

    with pytest.raises(ValueError, match="exactly one"):
        mcp_server.regulonado_build_dataset(
            "/data/intervals.bed",
            "/data/genome.fa",
            str(tmp_path / "dataset"),
        )


def test_regulonado_train_invokes_cli_dry_run(monkeypatch, tmp_path):
    mcp_server = _import_mcp_server(monkeypatch)
    monkeypatch.setattr(mcp_server.shutil, "which", lambda command: f"/bin/{command}")
    calls = []

    def fake_run(args):
        calls.append(args)
        return "python -m regulonado.train +experiment=head_only_borzoi data.path=/data/dataset"

    monkeypatch.setattr(mcp_server, "_run", fake_run)

    result = mcp_server.regulonado_train(
        "/data/dataset",
        output_dir=str(tmp_path / "run"),
        experiment="head_only_borzoi",
        max_steps=10,
        no_wandb=True,
        hydra_overrides=["trainer.max_eval_samples=2"],
    )

    assert calls == [
        [
            "regulonado",
            "train",
            "/data/dataset",
            "--experiment",
            "head_only_borzoi",
            "--nproc-per-node",
            "1",
            "--output-dir",
            str(tmp_path / "run"),
            "--max-steps",
            "10",
            "--no-wandb",
            "--dry-run",
            "trainer.max_eval_samples=2",
        ]
    ]
    assert "regulonado.train" in result


def test_regulonado_scaling_tools_invoke_cli(monkeypatch):
    mcp_server = _import_mcp_server(monkeypatch)
    monkeypatch.setattr(mcp_server.shutil, "which", lambda command: f"/bin/{command}")
    calls = []

    def fake_run(args):
        calls.append(args)
        return "done"

    monkeypatch.setattr(mcp_server, "_run", fake_run)

    assert mcp_server.regulonado_scale_bigwigs("/bw", "/out/scale.parquet") == "done"
    assert (
        mcp_server.regulonado_calculate_original_scaling(
            "/ds/regulonado_metadata.json",
            output_path="/ds/scale.parquet",
            workers=4,
        )
        == "done"
    )
    assert (
        mcp_server.regulonado_calculate_tmm_scaling(
            "/ds/regulonado_metadata.json",
            scale_factors="/ds/scale.parquet",
            output_path="/ds/tmm.parquet",
            split="validation",
        )
        == "done"
    )
    assert (
        mcp_server.regulonado_enrich_metadata(
            "/ds/regulonado_metadata.json",
            "/ds/tmm.parquet",
            fields=["scale_factor", "tmm_factor"],
        )
        == "done"
    )
    assert (
        mcp_server.regulonado_recompress_dataset(
            "/ds/raw",
            "/ds/rechunked",
            max_batch_size=4,
        )
        == "done"
    )

    assert calls == [
        [
            "regulonado",
            "scale",
            "/bw",
            "--output",
            "/out/scale.parquet",
            "--format",
            "parquet",
            "--workers",
            "16",
            "--glob",
            "*.bw",
        ],
        [
            "regulonado",
            "calculate-original-scaling",
            "/ds/regulonado_metadata.json",
            "--format",
            "parquet",
            "--workers",
            "4",
            "--output",
            "/ds/scale.parquet",
        ],
        [
            "regulonado",
            "calculate-tmm-scaling",
            "/ds/regulonado_metadata.json",
            "--format",
            "parquet",
            "--split",
            "validation",
            "--trim-m",
            "0.3",
            "--trim-a",
            "0.05",
            "--min-count",
            "1.0",
            "--scale-factors",
            "/ds/scale.parquet",
            "--output",
            "/ds/tmm.parquet",
        ],
        [
            "regulonado",
            "enrich-metadata",
            "/ds/regulonado_metadata.json",
            "/ds/tmm.parquet",
            "--field",
            "scale_factor",
            "--field",
            "tmm_factor",
        ],
        [
            "regulonado",
            "recompress-dataset",
            "/ds/raw",
            "/ds/rechunked",
            "--level",
            "3",
            "--workers",
            "4",
            "--max-batch-size",
            "4",
        ],
    ]
