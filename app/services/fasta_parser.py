from dataclasses import dataclass

ALLOWED_BASES = {"A", "C", "G", "T", "N"}
MAX_SEQUENCE_LENGTH = 1000


class FastaValidationError(Exception):
    def __init__(self, message: str, status_code: int = 400):
        self.message = message
        self.status_code = status_code
        super().__init__(message)


@dataclass(frozen=True)
class FastaRecord:
    sequence_id: str
    sequence: str


def parse_single_fasta(raw_bytes: bytes) -> FastaRecord:
    if not raw_bytes:
        raise FastaValidationError("Uploaded file is empty.")

    try:
        text = raw_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise FastaValidationError("File must be valid UTF-8 text.") from exc

    lines = text.splitlines()
    non_empty_lines = [line for line in lines if line.strip() != ""]

    if not non_empty_lines:
        raise FastaValidationError("FASTA file contains no content.")

    if not non_empty_lines[0].startswith(">"):
        raise FastaValidationError("Invalid FASTA format: must start with '>' header.")

    header_count = sum(1 for line in non_empty_lines if line.startswith(">"))
    if header_count > 1:
        raise FastaValidationError("Multiple FASTA records are not allowed.")

    header_body = non_empty_lines[0][1:].strip()
    if not header_body:
        raise FastaValidationError("FASTA header cannot be empty.")

    sequence_id = header_body.split()[0]
    seq_lines = non_empty_lines[1:]

    if not seq_lines:
        raise FastaValidationError("FASTA sequence cannot be empty.")

    for line in seq_lines:
        if any(ch.isspace() for ch in line):
            raise FastaValidationError("Embedded whitespace in sequence is not allowed.")

    sequence = "".join(seq_lines).upper()

    if len(sequence) == 0:
        raise FastaValidationError("FASTA sequence cannot be empty.")

    if len(sequence) > MAX_SEQUENCE_LENGTH:
        raise FastaValidationError(
            f"Sequence length ({len(sequence)}) exceeds maximum limit of {MAX_SEQUENCE_LENGTH} bp."
        )

    invalid_chars = sorted(set(sequence) - ALLOWED_BASES)
    if invalid_chars:
        raise FastaValidationError(
            f"Invalid nucleotide character(s) found: {', '.join(invalid_chars)}. Only A, C, G, T, N are allowed."
        )

    return FastaRecord(sequence_id=sequence_id, sequence=sequence)