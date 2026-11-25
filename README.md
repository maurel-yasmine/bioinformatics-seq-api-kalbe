# DNA Sequence Analysis REST API

A modular FastAPI service for **global DNA sequence alignment**, **variant calling (SNV, insertion, deletion)**, **3-frame protein translation**, and **biological variant effect annotation** against a reference coding sequence (`REF_001`).

## Features

- **Strict Single-Record FASTA Parser (`app/services/fasta_parser.py`):**
  - Enforces single-record FASTA structure, 1-based sequence ID extraction (first whitespace-delimited token after `>`), uppercase normalization, and length bounds (`1..1000` bp).
  - Supports multi-line sequences and blank lines while rejecting empty files, multi-FASTA files, embedded sequence whitespace, and non-`A/C/G/T/N` characters with structured HTTP `400` responses.
- **Pure-Python Needleman-Wunsch Global Aligner (`app/services/aligner.py`):**
  - Dynamic programming global alignment (`match = +1`, `mismatch = -1`, `gap = -2`, ambiguous `N = 0`).
  - Deterministic traceback tie-breaking priority: `Diagonal (match/mismatch/N) > Up (deletion) > Left (insertion)`.
  - Identity computed strictly as exact `A/C/G/T` matches divided by total alignment columns (`N/N` excluded from matches).
  - 1-based variant extraction grouping contiguous indels and sorting variants by reference position (`0` for leading insertions).
- **3-Frame DNA-to-Protein Translator (`app/services/translator.py`):**
  - Standard NCBI genetic code table supporting reading frames `1`, `2`, and `3` on the forward strand.
  - Translates any codon containing `N` as `X`, ignores trailing incomplete codons (`< 3` bp), and stops immediately before the first stop codon (`TAA`, `TAG`, `TGA`) without including `*` in the protein sequence.
- **Bonus Endpoints & Web UI:**
  - `POST /api/analyze`: Combines alignment, translation, and coding variant effect classification (`synonymous`, `missense`, `nonsense`, `frameshift`, `in-frame insertion`, `in-frame deletion`).
  - `GET /`: Lightweight browser UI for uploading FASTA files and inspecting JSON responses alongside interactive OpenAPI docs at `/docs`.

## Project Structure

```text
bioinformatics-seq-api-kalbe/
├── app/
│   ├── main.py                  # FastAPI entrypoint, exception handlers & UI route
│   ├── schemas.py               # Pydantic request/response models
│   ├── routers/
│   │   └── bio.py               # /api/align, /api/translate, /api/analyze endpoints
│   ├── services/
│   │   ├── fasta_parser.py      # Strict single-record FASTA validation
│   │   ├── aligner.py           # Needleman-Wunsch alignment & variant caller
│   │   ├── translator.py        # 3-frame genetic code translation
│   │   └── annotator.py         # Coding variant functional effect classifier
│   └── static/
│       └── index.html           # Simple interactive web client
├── data/
│   ├── reference.fasta          # Reference sequence (REF_001, 300 bp)
│   └── patient_01..05.fasta     # Patient test samples
├── tests/
│   └── test_api.py              # Automated pytest suite
├── Dockerfile
├── requirements.txt
└── README.md
```

## Setup & Running Locally

1. **Create and activate a virtual environment:**
   ```bash
   python -m venv .venv
   # Windows PowerShell:
   .venv\Scripts\Activate.ps1
   # Linux/macOS:
   source .venv/bin/activate
   ```

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Start the API server:**
   ```bash
   python -m uvicorn app.main:app --reload
   ```
   - Web UI: `http://127.0.0.1:8000/`
   - Interactive Swagger Docs: `http://127.0.0.1:8000/docs`

4. **Run automated tests:**
   ```bash
   python -m pytest -v
   ```
   Example terminal output (`13 passed`):
   ```text
   ========================= test session starts =========================
   platform win32 -- Python 3.14.0, pytest-9.2.0, pluggy-1.7.0
   collected 13 items

   tests/test_api.py::test_align_patient_01_identical PASSED        [  7%]
   tests/test_api.py::test_align_patient_02_snvs PASSED             [ 15%]
   tests/test_api.py::test_align_patient_03_deletions PASSED        [ 23%]
   tests/test_api.py::test_align_patient_04_ambiguous_n PASSED      [ 30%]
   tests/test_api.py::test_align_patient_05_insertions PASSED       [ 38%]
   tests/test_api.py::test_leading_insertion_position_zero PASSED   [ 46%]
   tests/test_api.py::test_translate_patient_02_default_frame PASSED [ 53%]
   tests/test_api.py::test_translate_reading_frames_and_stop_codons PASSED [ 61%]
   tests/test_api.py::test_translate_ambiguous_codons_and_short_sequence PASSED [ 69%]
   tests/test_api.py::test_fasta_formatting_normalization PASSED    [ 76%]
   tests/test_api.py::test_fasta_validation_errors PASSED           [ 84%]
   tests/test_api.py::test_translate_invalid_frame PASSED           [ 92%]
   tests/test_api.py::test_analyze_patient_02_and_03_variant_effects PASSED [100%]

   ==================== 13 passed, 1 warning in 0.95s ====================
   ```

## Running with Docker

```bash
docker build -t bioinformatics-seq-api .
docker run -p 8000:8000 bioinformatics-seq-api
```

## Example `curl` Requests

```bash
# Global alignment & variant detection
curl -X POST "[http://127.0.0.1:8000/api/align](http://127.0.0.1:8000/api/align)" \
  -F "file=@data/patient_02.fasta"

# Protein translation (frame 1)
curl -X POST "[http://127.0.0.1:8000/api/translate](http://127.0.0.1:8000/api/translate)" \
  -F "file=@data/patient_02.fasta" \
  -F "frame=1"

# Combined alignment, translation & variant effect annotation
curl -X POST "[http://127.0.0.1:8000/api/analyze](http://127.0.0.1:8000/api/analyze)" \
  -F "file=@data/patient_02.fasta" \
  -F "frame=1"
```

## Design Assumptions & Edge-Case Handling

- **Ambiguous Bases (`N`):** Scored as `0.0` during alignment so uncalled sequencer stretches (e.g., `patient_04`) do not trigger false gaps or false SNVs. Excluded from exact `A/C/G/T` match counts in `identity`. During translation, any codon containing `N` yields `X`.
- **Tied Optimal Alignments:** Broken deterministically in favor of `Diagonal > Up (deletion) > Left (insertion)`.
- **Variant Effects (`/api/analyze`):** Evaluated relative to the canonical reading frame (`frame = 1`) of `REF_001`.