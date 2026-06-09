import io

import pytest
from fastapi.testclient import TestClient

from two_hit.web_demo import create_app


@pytest.fixture
def client():
    return TestClient(create_app())


def test_index_lists_cases(client):
    r = client.get("/")
    assert r.status_code == 200
    assert "TP53" in r.text
    assert "/case/tp53_loh" in r.text


@pytest.mark.parametrize(
    "name,gene",
    [
        ("tp53_loh", "TP53"),
        ("kras_mut_only", "KRAS"),
        ("egfr_amp_mut", "EGFR"),
        ("brca_compound_het", "BRCA1"),
    ],
)
def test_case_renders_report(client, name, gene):
    r = client.get(f"/case/{name}")
    assert r.status_code == 200
    assert gene in r.text
    assert "<svg" in r.text


def test_unknown_case_404(client):
    r = client.get("/case/does_not_exist")
    assert r.status_code == 404


def test_analyze_rejects_oversized_file(client):
    big = io.BytesIO(b"x" * (1_000_001))
    small = io.BytesIO(b"ID\tchrom\tloc.start\tloc.end\tseg.mean\n")
    r = client.post(
        "/analyze",
        files={
            "maf_file": ("big.maf", big, "text/plain"),
            "seg_file": ("s.seg", small, "text/plain"),
        },
    )
    assert r.status_code == 400
    assert "1 MB" in r.text or "too large" in r.text.lower()


def test_analyze_valid_upload_returns_report(client):
    maf = (
        "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
        "Variant_Classification\tHGVSp_Short\tTumor_Sample_Barcode\t"
        "t_alt_count\tt_ref_count\n"
        "KRAS\t12\t25398284\t25398284\tMissense_Mutation\tp.G12D\tUP\t90\t110\n"
    )
    seg = (
        "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n"
        "UP\t12\t1\t133851895\t450\t0.05\n"
    )
    r = client.post(
        "/analyze",
        files={
            "maf_file": ("u.maf", io.BytesIO(maf.encode()), "text/plain"),
            "seg_file": ("u.seg", io.BytesIO(seg.encode()), "text/plain"),
        },
    )
    assert r.status_code == 200
    assert "KRAS" in r.text


def test_analyze_off_panel_only_shows_note(client):
    # NOTCH1 is not in the panel -> nothing to report
    maf = (
        "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
        "Variant_Classification\tHGVSp_Short\tTumor_Sample_Barcode\t"
        "t_alt_count\tt_ref_count\n"
        "NOTCH1\t9\t139399365\t139399365\tNonsense_Mutation\tp.R1598*\tUP\t40\t60\n"
    )
    seg = "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\nUP\t9\t1\t141213431\t450\t0.0\n"
    r = client.post(
        "/analyze",
        files={
            "maf_file": ("u.maf", io.BytesIO(maf.encode()), "text/plain"),
            "seg_file": ("u.seg", io.BytesIO(seg.encode()), "text/plain"),
        },
    )
    assert r.status_code == 200
    assert "no panel" in r.text.lower()
