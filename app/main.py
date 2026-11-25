from pathlib import Path
from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse, JSONResponse
from app.routers.bio import router as bio_router
from app.services.fasta_parser import FastaValidationError

app = FastAPI(
    title="DNA Sequence Analysis API",
    description="REST API for global DNA sequence alignment, variant detection, protein translation, and variant effect annotation.",
    version="1.0.0",
)

INDEX_HTML = Path(__file__).resolve().parent / "static" / "index.html"


@app.exception_handler(FastaValidationError)
async def fasta_validation_exception_handler(
    request: Request, exc: FastaValidationError
):
    return JSONResponse(status_code=exc.status_code, content={"detail": exc.message})


@app.get("/", response_class=HTMLResponse, include_in_schema=False)
async def serve_ui():
    return INDEX_HTML.read_text(encoding="utf-8")


app.include_router(bio_router)