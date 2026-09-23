#!/usr/bin/env python3
"""Overlay completed presentation and APPM evidence onto a canonical TSV."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalize_hla(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", value.upper().replace("HLA", ""))


def peptide(row: dict[str, str]) -> str:
    return str(row.get("peptide") or row.get("mutant_peptide") or "").strip().upper()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--deepimmuno", required=True, type=Path)
    parser.add_argument(
        "--mhcflurry",
        type=Path,
        help="Optional normalized MHCflurry TSV keyed by peptide and HLA allele.",
    )
    parser.add_argument("--netmhcstabpan", required=True, type=Path)
    parser.add_argument("--netchop", required=True, type=Path)
    parser.add_argument("--appm-modifiers", required=True, type=Path)
    parser.add_argument("--appm-flags", required=True, type=Path)
    parser.add_argument("--appm-summary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    deep = {
        (r.get("peptide", "").upper(), normalize_hla(r.get("hla_allele", ""))): r
        for r in read_tsv(args.deepimmuno)
        if r.get("deepimmuno_score")
    }
    mhcflurry = (
        {
            (r.get("peptide", "").upper(), normalize_hla(r.get("hla_allele", ""))): r
            for r in read_tsv(args.mhcflurry)
        }
        if args.mhcflurry
        else {}
    )
    stab = {
        (
            r.get("Peptide", r.get("peptide", "")).upper(),
            normalize_hla(r.get("HLA", r.get("hla_allele", ""))),
        ): r
        for r in read_tsv(args.netmhcstabpan)
    }
    netchop = {r.get("peptide_id", ""): r for r in read_tsv(args.netchop)}
    modifiers = {r.get("peptide_id", ""): r for r in read_tsv(args.appm_modifiers)}
    flags = {r.get("peptide_id", ""): r for r in read_tsv(args.appm_flags)}
    summary_rows = read_tsv(args.appm_summary)
    summary = summary_rows[0] if summary_rows else {}
    rows = read_tsv(args.input)

    appm_fields = [
        "appm_multiplier",
        "appm_multiplier_reason",
        "appm_integrity_status",
        "appm_evidence_completeness",
        "appm_review_required",
        "appm_action",
        "appm_call_confidence",
        "appm_call_confidence_score",
        "confidence_reason",
        "functional_validation_status",
        "restricting_locus_expression_status",
    ]
    fields = list(rows[0]) if rows else []
    additions = [
        "deepimmuno_score",
        "mhcflurry_affinity",
        "mhcflurry_affinity_percentile",
        "mhcflurry_processing_score",
        "mhcflurry_presentation_score",
        "mhcflurry_mt_affinity_percentile",
        "mhcflurry_mt_processing_score",
        "mhcflurry_mt_presentation_score",
        "mhcflurry_wt_affinity_percentile",
        "mhcflurry_wt_processing_score",
        "mhcflurry_wt_presentation_score",
        "mhcflurry_mt_wt_presentation_difference",
        "mhcflurry_allele_support_status",
        "netmhcstabpan_score",
        "netmhcstabpan_rank",
        "netmhcstabpan_wt_score",
        "netmhcstabpan_wt_rank",
        "netmhcstabpan_allele_support_status",
        "netchop_31d_max_score",
        "netchop_31d_mean_score",
        "netchop_31d_cterm_score",
        "netchop_31d_cleavage_sites",
        "netchop_processing_status",
        "tap_processing_status",
        "appm_mhc_i_integrity",
        "appm_mhc_ii_integrity",
        "b2m_status",
        "jak_stat_status",
        "nlrc5_status",
        "ciita_status",
        *appm_fields,
    ]
    for name in additions:
        if name not in fields:
            fields.append(name)

    for row in rows:
        pid = row.get("peptide_id", "")
        hla = normalize_hla(row.get("hla_allele", ""))
        mt = peptide(row)
        wt = str(row.get("wildtype_peptide") or row.get("wt_peptide") or "").strip().upper()
        if hit := deep.get((mt, hla)):
            row["deepimmuno_score"] = hit.get("deepimmuno_score", "")
        if hit := mhcflurry.get((mt, hla)):
            row["mhcflurry_affinity"] = hit.get("mhcflurry_affinity", "")
            for suffix in ("affinity_percentile", "processing_score", "presentation_score"):
                value = hit.get(f"mhcflurry_{suffix}", "")
                row[f"mhcflurry_{suffix}"] = value
                row[f"mhcflurry_mt_{suffix}"] = value
            row["mhcflurry_allele_support_status"] = "SUPPORTED"
        if wt and (hit := mhcflurry.get((wt, hla))):
            for suffix in ("affinity_percentile", "processing_score", "presentation_score"):
                row[f"mhcflurry_wt_{suffix}"] = hit.get(f"mhcflurry_{suffix}", "")
        try:
            row["mhcflurry_mt_wt_presentation_difference"] = str(
                float(row["mhcflurry_mt_presentation_score"])
                - float(row["mhcflurry_wt_presentation_score"])
            )
        except (KeyError, TypeError, ValueError):
            pass
        if hit := stab.get((mt, hla)):
            row["netmhcstabpan_score"] = hit.get("score", hit.get("netmhcstabpan_score", ""))
            row["netmhcstabpan_rank"] = hit.get(
                "percentile_rank", hit.get("netmhcstabpan_rank", "")
            )
            row["netmhcstabpan_allele_support_status"] = "SUPPORTED"
        if wt and (hit := stab.get((wt, hla))):
            row["netmhcstabpan_wt_score"] = hit.get(
                "score", hit.get("netmhcstabpan_score", "")
            )
            row["netmhcstabpan_wt_rank"] = hit.get(
                "percentile_rank", hit.get("netmhcstabpan_rank", "")
            )
        if hit := netchop.get(pid):
            for name in (
                "netchop_31d_max_score",
                "netchop_31d_mean_score",
                "netchop_31d_cterm_score",
                "netchop_31d_cleavage_sites",
                "netchop_processing_status",
            ):
                row[name] = hit.get(name, "")
        if hit := modifiers.get(pid):
            for name in appm_fields:
                row[name] = hit.get(name, "")
        if hit := flags.get(pid):
            row["tap_processing_status"] = hit.get(
                "tap_status", hit.get("tap_processing_status", "")
            )

        row["appm_mhc_i_integrity"] = summary.get("mhc_i_integrity_status", "")
        row["appm_mhc_ii_integrity"] = summary.get("mhc_ii_integrity_status", "")
        row["b2m_status"] = summary.get("b2m_risk", "")
        row["jak_stat_status"] = summary.get("ifng_response_status", "")
        row["nlrc5_status"] = summary.get("nlrc5_risk", "")
        row["ciita_status"] = summary.get("ciita_risk", "")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
