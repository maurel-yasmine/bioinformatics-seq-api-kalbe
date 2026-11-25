from typing import List, Literal
from pydantic import BaseModel


class AlignmentDetail(BaseModel):
    score: float
    identity: float
    reference_aligned: str
    query_aligned: str


class VariantItem(BaseModel):
    position: int
    ref: str
    alt: str
    type: Literal["SNV", "insertion", "deletion"]


class AlignResponse(BaseModel):
    sequence_id: str
    reference_id: str
    alignment: AlignmentDetail
    variants: List[VariantItem]


class TranslateResponse(BaseModel):
    sequence_id: str
    frame: int
    protein: str
    length: int
    has_stop_codon: bool


class AnnotatedVariantItem(BaseModel):
    position: int
    ref: str
    alt: str
    type: Literal["SNV", "insertion", "deletion"]
    effect: Literal[
        "synonymous",
        "missense",
        "nonsense",
        "frameshift",
        "in-frame insertion",
        "in-frame deletion",
        "undetermined",
    ]


class AnalyzeResponse(BaseModel):
    sequence_id: str
    reference_id: str
    alignment: AlignmentDetail
    translation: TranslateResponse
    annotated_variants: List[AnnotatedVariantItem]