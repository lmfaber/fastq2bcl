import re
import logging

_logger = logging.getLogger(__name__)

SEQDESC_RE = re.compile(
    r"(?P<instrument>[A-Za-z0-9_-]+):"
    + r"(?P<run_number>[0-9]+):"
    + r"(?P<flowcell_id>[A-Za-z0-9-]+):"
    + r"(?P<lane>[0-9]+):"
    + r"(?P<tile>[0-9]+):?"
    + r"(?P<x_pos>[0-9]+)?:?"
    + r"(?P<y_pos>[0-9]+)?:?"
    + r"(?P<UMI>[A-Z-]+)?"
    + r"\s"
    + r"(?P<read>[0-9]+):"
    + r"(?P<is_filtered>[YN]+):"
    + r"(?P<control_number>[0-9]+):"
    + r"(?P<index>[0-9A-Z+]+)"
)

VALID_KEYS = [
    "instrument",
    "run_number",
    "flowcell_id",
    "lane",
    "tile",
    "x_pos",
    "y_pos",
    "UMI",
    "read",
    "is_filtered",
    "control_number",
    "index",
]


def parse_seqdesc_fast_bytes(txt):
    """
    Parse hot-path fields from an Illumina FASTQ description as bytes.

    Returns ``(x_pos, y_pos, index, umi)`` where ``umi`` is ``None`` when absent.
    """
    first, sep, rest = txt.partition(b" ")
    if not sep:
        first, sep, rest = txt.partition(b"\t")
    if not sep:
        raise ValueError(f"Sequence identifier not recognized: {txt.decode(errors='replace')}")

    fields = first.split(b":")
    if len(fields) == 7:
        umi = None
    elif len(fields) == 8:
        umi = fields[7]
    else:
        raise ValueError(f"Sequence identifier not recognized: {txt.decode(errors='replace')}")

    if not fields[5] or not fields[6]:
        raise ValueError(f"Requested Key x_pos or y_pos not Found in fastq description")

    read_info = rest.split(None, 1)[0].split(b":", 3)
    if len(read_info) != 4 or not read_info[3]:
        raise ValueError(f"Sequence identifier not recognized: {txt.decode(errors='replace')}")

    return fields[5], fields[6], read_info[3], umi


def parse_seqdesc_fields(txt):
    """
    Parse the SeqIO description field using named groups.
    """
    match = SEQDESC_RE.match(txt)
    if not match:
        raise ValueError(f"Sequence identifier not recognized: {txt}")

    return validate_fields(match.groupdict())


def validate_fields(fields):
    """
    Validate the fields extracted from SeqIO description
    """
    _logger.debug(f"Verifying keys ...")

    for key in VALID_KEYS:
        if not fields[key]:
            if key == "UMI" and fields[key] == None:
                _logger.debug(f"Found None value for optional key {key}. This is ok.")
            else:
                raise ValueError(f"Requested Key {key} not Found in fastq description")

    return fields
