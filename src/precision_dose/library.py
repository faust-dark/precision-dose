"""Drug library loading, validation and hashing.

The library is data, not code, and it is the part a clinical pharmacist owns.
Three rules are enforced mechanically so curation cannot rot:

1. every parameter used in a calculation names a citation that exists,
2. an entry whose status is not "resolved" cannot be dosed (no guessed numbers),
3. the entry's content hash is computed over its canonical form and printed in
   every output, so a regimen traces back to the exact data that produced it.

Format is TOML parsed by the standard library (tomllib). YAML was the original
plan; TOML was chosen because it keeps the zero-runtime-dependency promise
without hand-rolling a parser, which would be another unaudited component in a
tool whose whole point is auditability.
"""

import hashlib
import json
import tomllib
from pathlib import Path

from .errors import E_LIBRARY_INVALID, E_LIBRARY_UNRESOLVED, E_UNKNOWN_DRUG, EngineError

DEFAULT_DIR = Path(__file__).with_name("library")

REQUIRED_TOP = (
    "schema_version", "drug", "display_name", "drug_class", "entry_version",
    "status", "decision_basis", "model", "targets", "safety", "citations",
    "provenance",
)

VALID_STRUCTURES = ("one_compartment", "two_compartment", "michaelis_menten", "none")


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def entry_hash(entry):
    return "sha256:" + hashlib.sha256(canonical(entry).encode("utf-8")).hexdigest()


def _problems(entry):
    out = []
    for key in REQUIRED_TOP:
        if key not in entry:
            out.append(f"missing top-level key {key!r}")
    if out:
        return out
    if entry["schema_version"] != "1.0":
        out.append(f"unsupported schema_version {entry['schema_version']!r}")
    if entry["status"] not in ("resolved", "unresolved"):
        out.append("status must be 'resolved' or 'unresolved'")
    if entry["decision_basis"] not in ("pk_target", "label_criteria"):
        out.append("decision_basis must be 'pk_target' or 'label_criteria'")

    model = entry["model"]
    structure = model.get("structure")
    if structure not in VALID_STRUCTURES:
        out.append(f"model.structure {structure!r} is not one of {VALID_STRUCTURES}")
    if structure in ("one_compartment", "two_compartment"):
        for k in ("cl_ref", "v1_ref"):
            if not isinstance(model.get(k), (int, float)):
                out.append(f"model.{k} must be numeric")
        if structure == "two_compartment":
            for k in ("v2_ref", "q_ref"):
                if not isinstance(model.get(k), (int, float)):
                    out.append(f"model.{k} must be numeric for a two-compartment model")
        if model.get("cl_ref", 0) <= 0 or model.get("v1_ref", 0) <= 0:
            out.append("cl_ref and v1_ref must be positive")
    if structure == "michaelis_menten":
        for k in ("vmax_ref_mg_day", "km_mg_L", "vd_ref_l_per_kg"):
            if not isinstance(model.get(k), (int, float)):
                out.append(f"model.{k} must be numeric for Michaelis-Menten")

    citations = entry["citations"]
    if not isinstance(citations, dict) or not citations:
        out.append("citations must be a non-empty table")
    sources = model.get("sources") or {}
    for param, cite in sources.items():
        if cite not in citations:
            out.append(f"model.sources.{param} references unknown citation {cite!r}")

    iiv = model.get("iiv") or {}
    for k, v in iiv.items():
        if k.endswith("_cv_percent") and (not isinstance(v, (int, float)) or v < 0):
            out.append(f"model.iiv.{k} must be a non-negative number")

    for i, tgt in enumerate(entry["targets"]):
        for k in ("id", "metric", "low", "high", "source"):
            if k not in tgt:
                out.append(f"targets[{i}] missing {k!r}")
        if tgt.get("low") is not None and tgt.get("high") is not None:
            if tgt["low"] >= tgt["high"]:
                out.append(f"targets[{i}] low must be < high")
    for i, s in enumerate(entry["safety"]):
        for k in ("id", "metric", "threshold", "source"):
            if k not in s:
                out.append(f"safety[{i}] missing {k!r}")

    prov = entry["provenance"]
    for k in ("curator", "curated_on"):
        if not prov.get(k):
            out.append(f"provenance.{k} is required")

    # unresolved placeholder marker
    text = canonical(entry)
    if "PLACEHOLDER" in text:
        out.append("entry contains a PLACEHOLDER value (still being curated)")
    return out


def validate_entry(entry, allow_unresolved=True):
    probs = _problems(entry)
    if probs:
        raise EngineError(
            E_LIBRARY_INVALID,
            f"library entry for {entry.get('drug', '?')} is invalid: " + "; ".join(probs),
            details={"problems": probs},
        )
    if entry["status"] != "resolved" and not allow_unresolved:
        raise EngineError(
            E_LIBRARY_UNRESOLVED,
            f"library entry for {entry['drug']} is {entry['status']}: "
            + "; ".join(entry.get("provenance", {}).get("known_gaps", ["see provenance"])),
        )
    return entry


def load_entry(path):
    with open(path, "rb") as fh:
        entry = tomllib.load(fh)
    return entry


def load_library(directory=None):
    directory = Path(directory) if directory else DEFAULT_DIR
    drugs_dir = directory if directory.name == "drugs" else directory / "drugs"
    if not drugs_dir.exists():
        raise EngineError(
            E_LIBRARY_INVALID, f"library directory not found: {drugs_dir}"
        )
    library = {}
    problems = {}
    for path in sorted(drugs_dir.glob("*.toml")):
        entry = load_entry(path)
        try:
            validate_entry(entry)
        except EngineError as exc:
            problems[entry.get("drug", path.stem)] = exc.details.get("problems", [exc.message])
            entry["_invalid"] = problems[entry.get("drug", path.stem)]
        entry["_path"] = str(path)
        entry["content_hash"] = entry_hash({k: v for k, v in entry.items() if not k.startswith("_")})
        library[entry["drug"]] = entry
    interactions = {}
    ipath = drugs_dir.parent / "interactions.toml"
    if ipath.exists():
        with open(ipath, "rb") as fh:
            interactions = tomllib.load(fh)
    return library, interactions, problems


def get_drug(library, name):
    key = name.strip().lower()
    if key not in library:
        known = ", ".join(sorted(library))
        raise EngineError(E_UNKNOWN_DRUG, f"unknown drug {name!r}; library has: {known}")
    entry = library[key]
    if entry.get("_invalid"):
        raise EngineError(
            E_LIBRARY_INVALID,
            f"entry for {key} failed validation: " + "; ".join(entry["_invalid"]),
            details={"problems": entry["_invalid"]},
        )
    validate_entry(entry, allow_unresolved=False)
    return entry


def list_library(library):
    rows = []
    for name in sorted(library):
        e = library[name]
        rows.append(
            {
                "drug": name,
                "display_name": e.get("display_name"),
                "class": e.get("drug_class"),
                "status": e.get("status"),
                "entry_version": e.get("entry_version"),
                "decision_basis": e.get("decision_basis"),
                "structure": (e.get("model") or {}).get("structure"),
                "targets": [t.get("id") for t in e.get("targets", [])],
                "content_hash": e.get("content_hash"),
            }
        )
    return rows
