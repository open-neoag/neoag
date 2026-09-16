from pathlib import Path

from neoag.utils import write_tsv
from neoag.validation import make_validation_plan
from neoag.validation_design import (
    classify_validation_mode,
    design_validation_row,
    minigene_to_long_peptide,
)


def _pep(**kwargs):
    base = {
        "peptide_id": "P1",
        "event_id": "E1",
        "gene": "GENE1",
        "peptide": "AAAAAAAAA",
        "wildtype_peptide": "AAABAAAAA",
        "hla_allele": "HLA-A*02:01",
        "event_type": "SNV",
        "peptide_consequence": "missense",
        "mhc_class": "I",
        "final_priority": "B",
        "safety_status": "PASS",
        "presentation_evidence_grade": "A",
        "appm_multiplier": "1.0",
        "ccf_multiplier": "1.0",
        "recommended_use": "",
    }
    base.update(kwargs)
    return base


def test_classify_missense_short_pair():
    assert classify_validation_mode(_pep()) == "missense_short_pair"


def test_classify_frameshift_long():
    p = _pep(event_type="InDel", peptide_consequence="frameshift", peptide="YLSFTLETL")
    assert classify_validation_mode(p) == "frameshift_long"


def test_classify_splice_junction_long():
    p = _pep(event_type="InDel", peptide_consequence="splice_junction", peptide="KLKLRRVKK")
    assert classify_validation_mode(p) == "splice_junction_long"


def test_classify_fusion_junction_long():
    p = _pep(event_type="Fusion", peptide_consequence="fusion", crosses_junction="yes")
    assert classify_validation_mode(p) == "fusion_junction_long"


def test_minigene_to_long_peptide():
    assert minigene_to_long_peptide("ABC|DEF|GHI") == "ABCDEFGHI"
    assert minigene_to_long_peptide("FLANK|NOVELTAIL") == "FLANKNOVELTAIL"


def test_design_missense_includes_wt_control():
    row = design_validation_row(_pep())
    assert row["validation_mode"] == "missense_short_pair"
    assert row["short_peptide_mt"] == "AAAAAAAAA"
    assert row["short_peptide_wt"] == "AAABAAAAA"
    assert "WT" in row["recommended_assay"] or "pair" in row["recommended_assay"]


def test_design_frameshift_uses_catalog_minigene(tmp_path):
    catalog = tmp_path / "variant_peptides.tsv"
    write_tsv(catalog, [{
        "variant_key": "NLRP4|chr19:55862081TCA>T",
        "mutant_peptide": "YLSFTLETL",
        "minigene": "GQSVLLFEVLFYQPDLK|YLSFTLET|LS",
        "minigene_nt": "AAA|BBB|CCC",
    }])
    peptide = _pep(
        event_id="NLRP4|chr19:55862081TCA>T",
        event_type="InDel",
        peptide_consequence="frameshift",
        peptide="YLSFTLETL",
        wildtype_peptide="YLSFTLTKL",
    )
    rows = make_validation_plan([peptide], peptide_catalog_tsv=catalog)
    assert rows[0]["validation_mode"] == "frameshift_long"
    assert rows[0]["minigene"] == "GQSVLLFEVLFYQPDLK|YLSFTLET|LS"
    assert rows[0]["long_peptide"] == "GQSVLLFEVLFYQPDLKYLSFTLETLS"
    assert "frameshift minigene" in rows[0]["recommended_assay"].lower()


def test_design_priority_d_do_not_advance():
    row = design_validation_row(_pep(final_priority="D"))
    assert row["validation_mode"] == "do_not_advance"
    assert row["recommended_assay"] == "Do not advance"


def test_consensus_r3_without_hard_failure_is_not_blocked_by_weighted_priority_d():
    row = design_validation_row(_pep(
        final_priority="D",
        evidence_grade="R3",
        event_type="InDel",
        peptide_consequence="frameshift",
        rna_support_state="RNA_UNASSESSED",
        rna_depth="0",
    ))
    assert row["validation_mode"] == "frameshift_long"
    assert row["recommended_assay"] != "Do not advance"
    assert "RNA locus evidence unassessed" in row["validation_notes"]


def test_consensus_r3_hard_failure_remains_blocked():
    row = design_validation_row(_pep(
        final_priority="D",
        evidence_grade="R3",
        hard_failure_codes="HARD_REFERENCE_PROTEOME_MATCH",
    ))
    assert row["validation_mode"] == "do_not_advance"


def test_event_level_best_evidence_grade_also_overrides_stale_weighted_d():
    row = design_validation_row(_pep(
        final_priority="D",
        best_evidence_grade="R3-GAP",
        event_type="InDel",
        peptide_consequence="frameshift",
    ))
    assert row["validation_mode"] == "frameshift_long"
