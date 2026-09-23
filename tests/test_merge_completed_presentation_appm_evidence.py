from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_merge_overlays_mhcflurry_mt_and_wt_scores(tmp_path: Path) -> None:
    script = Path(__file__).parents[1] / "scripts" / "merge_completed_presentation_appm_evidence.py"
    input_path = tmp_path / "input.tsv"
    mhcflurry_path = tmp_path / "mhcflurry.tsv"
    output_path = tmp_path / "output.tsv"

    _write_tsv(
        input_path,
        ["peptide_id", "peptide", "wildtype_peptide", "hla_allele", "mhcflurry_presentation_score"],
        [{
            "peptide_id": "p1",
            "peptide": "MUTPEPTID",
            "wildtype_peptide": "WILDPEPTI",
            "hla_allele": "HLA-A*02:01",
            "mhcflurry_presentation_score": "0.0",
        }],
    )
    mhcflurry_fields = [
        "peptide",
        "hla_allele",
        "mhcflurry_affinity",
        "mhcflurry_affinity_percentile",
        "mhcflurry_processing_score",
        "mhcflurry_presentation_score",
    ]
    _write_tsv(
        mhcflurry_path,
        mhcflurry_fields,
        [
            {
                "peptide": "MUTPEPTID",
                "hla_allele": "HLA-A*02:01",
                "mhcflurry_affinity": "42",
                "mhcflurry_affinity_percentile": "0.2",
                "mhcflurry_processing_score": "0.7",
                "mhcflurry_presentation_score": "0.8",
            },
            {
                "peptide": "WILDPEPTI",
                "hla_allele": "A0201",
                "mhcflurry_affinity": "420",
                "mhcflurry_affinity_percentile": "2.0",
                "mhcflurry_processing_score": "0.2",
                "mhcflurry_presentation_score": "0.3",
            },
        ],
    )

    empty_inputs = {
        "deepimmuno": ["peptide", "hla_allele", "deepimmuno_score"],
        "netmhcstabpan": ["peptide", "hla_allele", "netmhcstabpan_score"],
        "netchop": ["peptide_id", "netchop_31d_max_score"],
        "appm-modifiers": ["peptide_id", "appm_multiplier"],
        "appm-flags": ["peptide_id", "tap_status"],
    }
    paths: dict[str, Path] = {}
    for argument, fields in empty_inputs.items():
        paths[argument] = tmp_path / f"{argument}.tsv"
        _write_tsv(paths[argument], fields, [])
    summary_path = tmp_path / "appm-summary.tsv"
    _write_tsv(summary_path, ["mhc_i_integrity_status"], [{"mhc_i_integrity_status": "INTACT"}])

    command = [
        sys.executable,
        str(script),
        "--input",
        str(input_path),
        "--mhcflurry",
        str(mhcflurry_path),
        "--appm-summary",
        str(summary_path),
        "--output",
        str(output_path),
    ]
    for argument, path in paths.items():
        command.extend([f"--{argument}", str(path)])
    subprocess.run(command, check=True)

    with output_path.open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["mhcflurry_presentation_score"] == "0.8"
    assert row["mhcflurry_mt_presentation_score"] == "0.8"
    assert row["mhcflurry_wt_presentation_score"] == "0.3"
    assert float(row["mhcflurry_mt_wt_presentation_difference"]) == 0.5
    assert row["mhcflurry_allele_support_status"] == "SUPPORTED"
