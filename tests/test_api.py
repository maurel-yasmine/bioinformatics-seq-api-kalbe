from pathlib import Path
from fastapi.testclient import TestClient
from app.main import app
from app.services.aligner import extract_variants, needleman_wunsch

client = TestClient(app)
DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def test_align_patient_01_identical():
    with open(DATA_DIR / "patient_01.fasta", "rb") as f:
        response = client.post(
            "/api/align", files={"file": ("patient_01.fasta", f, "text/plain")}
        )

    assert response.status_code == 200
    body = response.json()
    assert body["sequence_id"] == "patient_01"
    assert body["reference_id"] == "REF_001"
    assert body["alignment"]["score"] == 300.0
    assert body["alignment"]["identity"] == 1.0
    assert body["variants"] == []


def test_align_patient_02_snvs():
    with open(DATA_DIR / "patient_02.fasta", "rb") as f:
        response = client.post(
            "/api/align", files={"file": ("patient_02.fasta", f, "text/plain")}
        )

    assert response.status_code == 200
    body = response.json()
    assert body["alignment"]["score"] == 294.0
    assert body["alignment"]["identity"] == 0.99
    assert body["variants"] == [
        {"position": 52, "ref": "C", "alt": "A", "type": "SNV"},
        {"position": 177, "ref": "A", "alt": "T", "type": "SNV"},
        {"position": 245, "ref": "T", "alt": "G", "type": "SNV"},
    ]


def test_align_patient_03_deletions():
    with open(DATA_DIR / "patient_03.fasta", "rb") as f:
        response = client.post(
            "/api/align", files={"file": ("patient_03.fasta", f, "text/plain")}
        )

    assert response.status_code == 200
    body = response.json()
    assert body["alignment"]["score"] == 288.0
    assert body["alignment"]["identity"] == 0.986667
    assert body["variants"] == [
        {"position": 106, "ref": "AAC", "alt": "-", "type": "deletion"},
        {"position": 208, "ref": "A", "alt": "-", "type": "deletion"},
    ]


def test_align_patient_04_ambiguous_n():
    with open(DATA_DIR / "patient_04.fasta", "rb") as f:
        response = client.post(
            "/api/align", files={"file": ("patient_04.fasta", f, "text/plain")}
        )

    assert response.status_code == 200
    body = response.json()
    assert body["alignment"]["score"] == 130.0
    assert body["alignment"]["identity"] == 0.433333
    assert body["variants"] == []


def test_align_patient_05_insertions():
    with open(DATA_DIR / "patient_05.fasta", "rb") as f:
        response = client.post(
            "/api/align", files={"file": ("patient_05.fasta", f, "text/plain")}
        )

    assert response.status_code == 200
    body = response.json()
    assert body["alignment"]["score"] == 292.0
    assert body["alignment"]["identity"] == 0.986842
    assert len(body["variants"]) == 3
    assert all(v["type"] == "insertion" for v in body["variants"])


def test_leading_insertion_position_zero():
    score, ref_aligned, query_aligned = needleman_wunsch("ACGT", "TTACGT")
    assert score == 0.0
    variants = extract_variants(ref_aligned, query_aligned)
    assert variants == [{"position": 0, "ref": "-", "alt": "TT", "type": "insertion"}]


def test_translate_patient_02_default_frame():
    with open(DATA_DIR / "patient_02.fasta", "rb") as f:
        response = client.post(
            "/api/translate", files={"file": ("patient_02.fasta", f, "text/plain")}
        )

    assert response.status_code == 200
    body = response.json()
    assert body["sequence_id"] == "patient_02"
    assert body["frame"] == 1
    assert body["length"] == 81
    assert body["has_stop_codon"] is True
    assert "*" not in body["protein"]
    assert (
        body["protein"]
        == "MSSFVPLLYVSRVGTSNRLFSLLTKFTWVSVTRLRNPAEFAACVRHPYLAHSSQSHEHNLALVITSGQSTPDSKGEQNWVS"
    )


def test_translate_reading_frames_and_stop_codons():
    fasta_content = b">SEQ_FRAMES\nAATGAAATAA\n"

    res1 = client.post(
        "/api/translate",
        files={"file": ("frames.fasta", fasta_content, "text/plain")},
        data={"frame": "1"},
    )
    assert res1.status_code == 200
    assert res1.json() == {
        "sequence_id": "SEQ_FRAMES",
        "frame": 1,
        "protein": "NEI",
        "length": 3,
        "has_stop_codon": False,
    }

    res2 = client.post(
        "/api/translate",
        files={"file": ("frames.fasta", fasta_content, "text/plain")},
        data={"frame": "2"},
    )
    assert res2.status_code == 200
    assert res2.json() == {
        "sequence_id": "SEQ_FRAMES",
        "frame": 2,
        "protein": "MK",
        "length": 2,
        "has_stop_codon": True,
    }

    res3 = client.post(
        "/api/translate",
        files={"file": ("frames.fasta", fasta_content, "text/plain")},
        data={"frame": "3"},
    )
    assert res3.status_code == 200
    assert res3.json() == {
        "sequence_id": "SEQ_FRAMES",
        "frame": 3,
        "protein": "",
        "length": 0,
        "has_stop_codon": True,
    }


def test_translate_ambiguous_codons_and_short_sequence():
    ambiguous_fasta = b">AMBIG\nACNGGNNAT\n"
    res_ambig = client.post(
        "/api/translate",
        files={"file": ("ambig.fasta", ambiguous_fasta, "text/plain")},
    )
    assert res_ambig.status_code == 200
    assert res_ambig.json()["protein"] == "XXX"
    assert res_ambig.json()["length"] == 3

    short_fasta = b">SHORT\nAT\n"
    res_short = client.post(
        "/api/translate",
        files={"file": ("short.fasta", short_fasta, "text/plain")},
    )
    assert res_short.status_code == 200
    assert res_short.json() == {
        "sequence_id": "SHORT",
        "frame": 1,
        "protein": "",
        "length": 0,
        "has_stop_codon": False,
    }


def test_fasta_formatting_normalization():
    messy_fasta = b">CUSTOM_ID extra header description\n\natg\ncgt\n"
    res = client.post(
        "/api/translate",
        files={"file": ("messy.fasta", messy_fasta, "text/plain")},
    )
    assert res.status_code == 200
    assert res.json()["sequence_id"] == "CUSTOM_ID"
    assert res.json()["protein"] == "MR"


def test_fasta_validation_errors():
    invalid_cases = [
        b"",
        b"ATGCGT",
        b">\nATGCGT",
        b">SEQ1\nATGC\n>SEQ2\nCGTA",
        b">SPACE\nATG CGT",
        b">RNA\nAUGC",
        b">IUPAC\nATGR",
        b">TOO_LONG\n" + (b"A" * 1001),
    ]

    for payload in invalid_cases:
        res = client.post(
            "/api/align",
            files={"file": ("bad.fasta", payload, "text/plain")},
        )
        assert res.status_code == 400


def test_translate_invalid_frame():
    res = client.post(
        "/api/translate",
        files={"file": ("ok.fasta", b">OK\nATGGCC", "text/plain")},
        data={"frame": "4"},
    )
    assert res.status_code == 422