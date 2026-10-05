"""Verified parser for the complete corrected Bordeaux totally real tables.

Standard library only.  The table range/completeness is a published external
input.  Source bytes, every GP row, coefficient order, class-group data, and
the inclusive manuscript cutoffs are checked here.  Full downloaded table
texts are bundled, so regeneration does not need network access.
"""
from pathlib import Path
from math import prod
import ast
import gzip
import hashlib
import json

TABLE_DIRECTORY = Path(__file__).resolve().parent / "data" / "number_field_tables"
MANIFEST_PATH = TABLE_DIRECTORY / "manifest.json"
ENUMERATION_BOUNDS = {3: 4218, 4: 68133, 5: 209508, 6: 2429064}
EXPECTED_COUNTS = {3: 143, 4: 552, 5: 37, 6: 40, 7: 4}
SEPTIC_BOUND = 28162991


def load_manifest():
    manifest = json.loads(MANIFEST_PATH.read_text())
    if manifest["schema"] != "bordeaux-field-table-input-v1":
        raise ValueError("unsupported source-table manifest schema")
    return manifest


def table_rows(degree, cutoff=None):
    """Return every verified source row through the inclusive cutoff.

    The coefficient vector returned in ``polynomial`` is constant first,
    as required by Sage's polynomial constructor and the existing JSON API.
    Table line numbers and the original descending vector remain explicit.
    Tables must be sorted, but equal discriminants are retained separately.
    """
    manifest = load_manifest()
    source = manifest["tables"][str(degree)]
    cutoff = source["manuscript_inclusive_cutoff"] if cutoff is None else cutoff
    if not isinstance(cutoff, int) or cutoff < 1:
        raise ValueError("cutoff must be a positive integer")
    if cutoff >= source["published_discriminant_bound_exclusive"]:
        raise ValueError("requested cutoff exceeds the complete source range")
    packed = (TABLE_DIRECTORY / source["file"]).read_bytes()
    if hashlib.sha256(packed).hexdigest() != source["repacked_gzip_sha256"]:
        raise ValueError("source gzip checksum mismatch")
    raw = gzip.decompress(packed)
    if len(raw) != source["decompressed_bytes"]:
        raise ValueError("source GP text byte count mismatch")
    if hashlib.sha256(raw).hexdigest() != source["source_text_sha256"]:
        raise ValueError("source GP text checksum mismatch")
    selected = []
    seen = {}
    duplicates = []
    previous = 0
    row_count = 0
    first_discriminant = None
    for line_number, text in enumerate(raw.decode("ascii").splitlines(), 1):
        if not text.strip():
            continue
        row = ast.literal_eval(text)
        if not isinstance(row, list) or len(row) != 4:
            raise ValueError("invalid GP row at line %s" % line_number)
        discriminant, descending, class_number, elementary_divisors = row
        if type(discriminant) is not int or discriminant < previous:
            raise ValueError("invalid or unsorted field discriminant")
        if not (0 < discriminant < source["published_discriminant_bound_exclusive"]):
            raise ValueError("field outside the published table range")
        if not (isinstance(descending, list) and len(descending) == degree + 1
                and descending[0] == 1
                and all(type(c) is int for c in descending)):
            raise ValueError("invalid monic defining polynomial")
        if type(class_number) is not int or class_number < 1:
            raise ValueError("invalid ordinary class number")
        if not (isinstance(elementary_divisors, list)
                and all(type(v) is int and v >= 2 for v in elementary_divisors)
                and prod(elementary_divisors) == class_number
                and all(a % b == 0 for a, b in zip(elementary_divisors,
                                                  elementary_divisors[1:]))):
            raise ValueError("invalid class-group elementary divisors")
        key = (discriminant, tuple(descending))
        if key in seen:
            duplicates.append({"first_line": seen[key],
                               "repeated_line": line_number,
                               "discriminant": discriminant})
            if discriminant <= cutoff:
                raise ValueError("duplicate field within the requested cutoff")
        seen[key] = line_number
        previous = discriminant
        row_count += 1
        if first_discriminant is None:
            first_discriminant = discriminant
        if discriminant <= cutoff:
            selected.append({
                "degree": degree,
                "discriminant": discriminant,
                "polynomial": list(reversed(descending)),
                "table_ordinary_class_number": class_number,
                "table_class_group_invariants": elementary_divisors,
                "table_file": source["file"],
                "table_line": line_number,
                "table_polynomial_highest_degree_first": descending,
            })
    if not (row_count == source["row_count"]
            and first_discriminant == source["first_discriminant"]
            and previous == source["last_discriminant"]
            and duplicates == source["exact_duplicate_source_rows"]):
        raise ValueError("source table extent mismatch")
    if cutoff == source["manuscript_inclusive_cutoff"]:
        if len(selected) != EXPECTED_COUNTS[degree]:
            raise ValueError("manuscript cutoff count mismatch")
    return selected
