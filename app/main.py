from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse
from app.routers.bio import router as bio_router
from app.services.fasta_parser import FastaValidationError

app = FastAPI(
    title="DNA Sequence Analysis API",
    description="REST API for global DNA sequence alignment, variant detection, and protein translation.",
    version="1.0.0",
)


@app.exception_handler(FastaValidationError)
async def fasta_validation_exception_handler(
    request: Request, exc: FastaValidationError
):
    return JSONResponse(status_code=exc.status_code, content={"detail": exc.message})


app.include_router(bio_router)