from pathlib import Path
from fastapi import APIRouter, File, Form, UploadFile
from app.schemas import (
    AlignResponse,
    AlignmentDetail,
    AnalyzeResponse,
    TranslateResponse,
)
from app.services.aligner import calculate_identity, extract_variants, needleman_wunsch
from app.services.annotator import annotate_variants
from app.services.fasta_parser import FastaValidationError, parse_single_fasta
from app.services.translator import translate_sequence

router = APIRouter(prefix="/api", tags=["Bioinformatics"])

REFERENCE_PATH = Path(__file__).resolve().parents[2] / "data" / "reference.fasta"
REFERENCE_RECORD = parse_single_fasta(REFERENCE_PATH.read_bytes())


@router.post("/align", response_model=AlignResponse)
async def align_sequence(file: UploadFile = File(...)):
    raw_bytes = await file.read()
    query_record = parse_single_fasta(raw_bytes)

    score, ref_aligned, query_aligned = needleman_wunsch(
        REFERENCE_RECORD.sequence, query_record.sequence
    )
    identity = calculate_identity(ref_aligned, query_aligned)
    variants = extract_variants(ref_aligned, query_aligned)

    return AlignResponse(
        sequence_id=query_record.sequence_id,
        reference_id=REFERENCE_RECORD.sequence_id,
        alignment=AlignmentDetail(
            score=score,
            identity=identity,
            reference_aligned=ref_aligned,
            query_aligned=query_aligned,
        ),
        variants=variants,
    )


@router.post("/translate", response_model=TranslateResponse)
async def translate_fasta(
    file: UploadFile = File(...),
    frame: int = Form(1),
):
    if frame not in (1, 2, 3):
        raise FastaValidationError(
            "Invalid frame value. Frame must be 1, 2, or 3.", status_code=422
        )

    raw_bytes = await file.read()
    query_record = parse_single_fasta(raw_bytes)

    protein, length, has_stop_codon = translate_sequence(
        query_record.sequence, frame=frame
    )

    return TranslateResponse(
        sequence_id=query_record.sequence_id,
        frame=frame,
        protein=protein,
        length=length,
        has_stop_codon=has_stop_codon,
    )


@router.post("/analyze", response_model=AnalyzeResponse)
async def analyze_sequence(
    file: UploadFile = File(...),
    frame: int = Form(1),
):
    if frame not in (1, 2, 3):
        raise FastaValidationError(
            "Invalid frame value. Frame must be 1, 2, or 3.", status_code=422
        )

    raw_bytes = await file.read()
    query_record = parse_single_fasta(raw_bytes)

    score, ref_aligned, query_aligned = needleman_wunsch(
        REFERENCE_RECORD.sequence, query_record.sequence
    )
    identity = calculate_identity(ref_aligned, query_aligned)
    variants = extract_variants(ref_aligned, query_aligned)
    annotated = annotate_variants(REFERENCE_RECORD.sequence, variants)

    protein, length, has_stop_codon = translate_sequence(
        query_record.sequence, frame=frame
    )

    return AnalyzeResponse(
        sequence_id=query_record.sequence_id,
        reference_id=REFERENCE_RECORD.sequence_id,
        alignment=AlignmentDetail(
            score=score,
            identity=identity,
            reference_aligned=ref_aligned,
            query_aligned=query_aligned,
        ),
        translation=TranslateResponse(
            sequence_id=query_record.sequence_id,
            frame=frame,
            protein=protein,
            length=length,
            has_stop_codon=has_stop_codon,
        ),
        annotated_variants=annotated,
    )