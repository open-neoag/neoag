from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "scripts" / "build_longrna_unique_splice_peptides.py"
MERGE = ROOT / "scripts" / "merge_longrna_comprehensive_evidence.py"
PROVENANCE = ROOT / "scripts" / "capture_longrna_provenance.py"
FULL_RUNNER = ROOT / "scripts" / "run_longrna_splice_full.sh"


def write(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def isoquant_dry_run_command(tmp_path: Path, result_dir: Path) -> list[str]:
    reference = write(tmp_path / "reference.fa", ">chr1\nACGT\n")
    annotation = write(tmp_path / "annotation.gtf", "annotation\n")
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    shortread = write(tmp_path / "shortread.tsv", "junction_key\n")
    sqanti_home = tmp_path / "SQANTI3"
    write(sqanti_home / "sqanti3_qc.py", "# fixture\n")
    return [
        "bash",
        str(FULL_RUNNER),
        "--start-from",
        "isoquant",
        "--isoquant-result-dir",
        str(result_dir),
        "--sample-id",
        "CASE-42",
        "--workdir",
        str(tmp_path / "work"),
        "--reference-fasta",
        str(reference),
        "--annotation-gtf",
        str(annotation),
        "--hla",
        str(hla),
        "--shortread-comparisons",
        str(shortread),
        "--sqanti-home",
        str(sqanti_home),
        "--dry-run",
    ]


def candidate_command(tmp_path: Path, *, normal: bool) -> tuple[list[str], Path]:
    classification = write(
        tmp_path / "classification.tsv",
        "isoform\tstructural_category\tCDS_type\tpredicted_NMD\tpsauron_score\tCDS_start\tCDS_end\tassociated_gene\n"
        "iso1\tnovel_not_in_catalog\tcomplete\tFALSE\t0.95\t1\t300\tGENE1\n",
    )
    junctions = write(
        tmp_path / "junctions.tsv",
        "isoform\tchrom\tgenomic_start_coord\tgenomic_end_coord\tstrand\tjunction_category\t"
        "total_coverage_unique\ttotal_coverage_multi\tcanonical\tjunction_number\n"
        "iso1\t1\t200\t300\t+\tnovel\t5\t6\tcanonical\t1\n",
    )
    gtf = write(
        tmp_path / "corrected.gtf",
        'chr1\ttest\texon\t100\t199\t.\t+\t.\ttranscript_id "iso1";\n'
        'chr1\ttest\texon\t301\t500\t.\t+\t.\ttranscript_id "iso1";\n',
    )
    proteins = write(tmp_path / "proteins.faa", ">iso1.p1\n" + "A" * 100 + "\n")
    shortread = write(tmp_path / "shortread.tsv", "junction_key\n")
    hla = write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    outdir = tmp_path / ("with_normal" if normal else "without_normal")
    command = [
        sys.executable,
        str(BUILD),
        "--sample-id",
        "CASE-42",
        "--sqanti-classification",
        str(classification),
        "--sqanti-junctions",
        str(junctions),
        "--corrected-gtf",
        str(gtf),
        "--td2-peptides",
        str(proteins),
        "--shortread-comparisons",
        str(shortread),
        "--hla",
        str(hla),
        "--outdir",
        str(outdir),
    ]
    if normal:
        normal_path = write(tmp_path / "normal.tsv", "junction_key\nchr1:199:301:+\n")
        command.extend(["--normal-junctions", str(normal_path)])
    return command, outdir


def test_candidate_builder_propagates_sample_id(tmp_path: Path) -> None:
    command, outdir = candidate_command(tmp_path, normal=False)
    subprocess.run(command, check=True, capture_output=True, text=True)

    rows = (outdir / "longrna_unique_junction_peptides.tsv").read_text().splitlines()
    assert len(rows) > 1
    assert all(row.startswith("CASE-42\t") for row in rows[1:])
    summary = json.loads((outdir / "longrna_unique_candidate_summary.json").read_text())
    assert summary["normal_background_status"] == "UNASSESSED"
    assert summary["longrna_unique_abnormal_junction_rows"] == 1


def test_candidate_builder_excludes_normal_background(tmp_path: Path) -> None:
    command, outdir = candidate_command(tmp_path, normal=True)
    subprocess.run(command, check=True, capture_output=True, text=True)

    summary = json.loads((outdir / "longrna_unique_candidate_summary.json").read_text())
    assert summary["normal_background_status"] == "ASSESSED"
    assert summary["normal_background_junction_keys"] == 1
    assert summary["longrna_unique_abnormal_junction_rows"] == 0


def test_longrna_shell_entrypoints_have_valid_syntax() -> None:
    for name in (
        "run_longrna_splice_full.sh",
        "run_longrna_downstream.sh",
        "run_longrna_downstream_260F400881011.sh",
    ):
        subprocess.run(["bash", "-n", str(ROOT / "scripts" / name)], check=True)


def test_evidence_merge_accepts_hla_allele_column(tmp_path: Path) -> None:
    provenance = write(
        tmp_path / "provenance.tsv",
        "sample_id\tpeptide\thla_allele\tlongrna_support\nCASE-42\tAAAAAAAA\tHLA-A*02:01\ttrue\n",
    )
    presentation = write(
        tmp_path / "presentation.tsv",
        "sample_id\tpeptide_id\tevent_id\tpeptide\thla_allele\tpresentation_evidence_grade\n"
        "CASE-42\tp1\te1\tAAAAAAAA\tHLA-A*02:01\tA\n",
    )
    output = tmp_path / "merged.tsv"
    summary = tmp_path / "summary.json"
    subprocess.run(
        [
            sys.executable,
            str(MERGE),
            "--provenance",
            str(provenance),
            "--presentation",
            str(presentation),
            "--output",
            str(output),
            "--summary",
            str(summary),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    result = json.loads(summary.read_text())
    assert result["complete_rows"] == 1
    assert result["missing_prediction_rows"] == 0


def test_provenance_retention_controls_hashing(tmp_path: Path) -> None:
    source = write(tmp_path / "input.txt", "auditable input\n")
    for retention, expects_hash in (("standard", False), ("full", True)):
        output = tmp_path / f"{retention}.json"
        subprocess.run(
            [
                sys.executable,
                str(PROVENANCE),
                "--sample-id",
                "CASE-42",
                "--retention",
                retention,
                "--scope",
                "output",
                "--output",
                str(output),
                "--path",
                str(source),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        manifest = json.loads(output.read_text())
        assert manifest["retention"] == retention
        assert ("sha256" in manifest["files"][0]) is expects_hash


def test_downstream_resume_validates_outputs(tmp_path: Path) -> None:
    work = tmp_path / "work"
    isoquant = work / "isoquant" / "CASE-42"
    sqanti_out = work / "sqanti3"
    sqanti_home = tmp_path / "SQANTI3"
    isoquant.mkdir(parents=True)
    sqanti_out.mkdir(parents=True)
    sqanti_home.mkdir()
    write(isoquant / "CASE-42.transcript_models.gtf", "transcript\n")
    write(sqanti_home / "sqanti3_qc.py", "# fixture\n")
    reference = write(tmp_path / "reference.fa", ">chr1\nACGT\n")
    annotation = write(tmp_path / "annotation.gtf", "annotation\n")
    for suffix in ("classification.txt", "junctions.txt", "corrected.gtf", "corrected.fasta"):
        write(sqanti_out / f"CASE-42_{suffix}", "fixture\n")
    protein = write(sqanti_out / "CASE-42_corrected.faa", ">iso1\nAAAA\n")
    write(work / ".sqanti3_complete", "complete\n")

    command = [
        "bash",
        str(ROOT / "scripts" / "run_longrna_downstream.sh"),
        "--sample-id",
        "CASE-42",
        "--workdir",
        str(work),
        "--reference-fasta",
        str(reference),
        "--annotation-gtf",
        str(annotation),
        "--sqanti-home",
        str(sqanti_home),
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    assert (sqanti_out / "protein_output.path").read_text().strip() == str(protein)

    protein.unlink()
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    assert result.returncode != 0
    assert "outputs failed validation" in (work / "longrna_splice_status.tsv").read_text()


def test_full_runner_accepts_existing_isoquant_result_without_fastq(tmp_path: Path) -> None:
    result_dir = tmp_path / "existing_isoquant"
    isoform_gtf = write(
        result_dir / "CASE-42" / "CASE-42.transcript_models.gtf",
        "transcript\n",
    )
    result = subprocess.run(
        isoquant_dry_run_command(tmp_path, result_dir),
        check=True,
        capture_output=True,
        text=True,
    )

    assert "start_from=isoquant" in result.stdout
    assert "isoquant_stage=REUSED" in result.stdout
    assert f"isoform_gtf={isoform_gtf}" in result.stdout
    assert "isoquant_command=" not in result.stdout


def test_full_runner_rejects_ambiguous_isoquant_results(tmp_path: Path) -> None:
    result_dir = tmp_path / "existing_isoquant"
    write(result_dir / "run1" / "one.transcript_models.gtf", "one\n")
    write(result_dir / "run2" / "two.transcript_models.gtf", "two\n")
    result = subprocess.run(
        isoquant_dry_run_command(tmp_path, result_dir),
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "expected exactly one non-empty" in result.stderr


def test_verifier_accepts_external_isoquant_gtf_pointer(tmp_path: Path) -> None:
    work = tmp_path / "work"
    external_gtf = write(tmp_path / "isoquant" / "CASE-42.transcript_models.gtf", "transcript\n")
    write(work / "longrna_splice_status.tsv", "stage\tstatus\tupdated_at\tdetail\n")
    write(work / "README.workflow.txt", "sample_id=CASE-42\n")
    write(work / "isoquant_transcript_model.path", f"{external_gtf}\n")
    write(work / "sqanti3" / "CASE-42_classification.txt", "classification\n")
    write(work / "translation" / "CASE-42.pep", ">iso1\nAAAA\n")
    write(work / "comprehensive_evidence.tsv", "sample_id\n")

    result = subprocess.run(
        ["bash", str(ROOT / "scripts" / "verify_longrna_splice.sh"), str(work)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "workflow_status=PASS" in result.stdout
    assert str(external_gtf) in result.stdout


def test_full_runner_rejects_resume_mode_mismatch(tmp_path: Path) -> None:
    result_dir = tmp_path / "existing_isoquant"
    write(result_dir / "CASE-42" / "CASE-42.transcript_models.gtf", "transcript\n")
    work = tmp_path / "work"
    write(work / "longrna_splice_status.tsv", "stage\tstatus\tupdated_at\tdetail\n")
    write(work / "README.workflow.txt", "sample_id=CASE-42\nstart_from=fastq\n")
    command = isoquant_dry_run_command(tmp_path, result_dir)
    command.remove("--dry-run")
    command.append("--resume")

    result = subprocess.run(command, check=False, capture_output=True, text=True)
    assert result.returncode != 0
    assert "resume mode mismatch" in result.stderr
