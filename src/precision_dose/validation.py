"""Validation harness: reproduce published regimens and guideline targets.

A case is a JSON document with a patient, a drug and a list of checks against
named output fields. Passing is a gate (non-zero exit), not a report. Where a
cited source does not supply the number a check would need, the check is
declared `source_verified: false` and the suite says so in its output rather
than quietly skipping it.
"""

import json
from pathlib import Path

from . import engine
from .errors import EngineError

DEFAULT_CASE_DIRS = (
    Path("validation/cases"),
    Path(__file__).resolve().parents[2] / "validation" / "cases",
)


def _find_cases(cases_dir):
    if cases_dir:
        p = Path(cases_dir)
        if not p.exists():
            raise EngineError("E_SCHEMA", f"case directory not found: {cases_dir}")
        return sorted(p.glob("*.json"))
    for d in DEFAULT_CASE_DIRS:
        if d.exists():
            return sorted(d.glob("*.json"))
    raise EngineError("E_SCHEMA", "no validation/cases directory found")


def _lookup(obj, path):
    cur = obj
    for part in path.split("."):
        if isinstance(cur, list):
            cur = cur[int(part)]
        elif isinstance(cur, dict):
            if part not in cur:
                return None
            cur = cur[part]
        else:
            return None
    return cur


def _check(value, check):
    op = check["op"]
    target = check.get("value")
    if op == "is_null":
        return value is None, f"{value!r} is null: {value is None}"
    if value is None:
        return False, "value is missing"
    if op == "eq":
        ok = value == target
    elif op == "ge":
        ok = value >= target
    elif op == "le":
        ok = value <= target
    elif op == "between":
        ok = target[0] <= value <= target[1]
    elif op == "approx":
        ok = abs(value - target) <= check.get("tolerance", 1e-6)
    elif op == "in":
        ok = value in target
    elif op == "max_abs_diff":
        ok = abs(value - target) <= check.get("tolerance", 1e-6)
    else:
        return False, f"unknown op {op!r}"
    return ok, f"{value!r} {op} {target!r}"


def run_cases(cases_dir=None, library_dir=None):
    rows = []
    for path in _find_cases(cases_dir):
        case = json.loads(path.read_text(encoding="utf-8"))
        row = {
            "case": case.get("id", path.stem),
            "description": case.get("description"),
            "citation": case.get("citation"),
            "source_verified": case.get("source_verified", True),
            "checks": [],
            "passed": False,
            "error": None,
        }
        try:
            result = engine.run_dose(
                case["patient"],
                case["drug"],
                target_id=case.get("target"),
                seed=case.get("seed", engine.DEFAULT_SEED),
                pta_n=case.get("pta_n", 2000),
                library_dir=library_dir,
            )
            row["status"] = result.get("status")
            for check in case.get("checks", []):
                value = _lookup(result, check["path"])
                ok, detail = _check(value, check)
                row["checks"].append(
                    {
                        "path": check["path"],
                        "ok": ok,
                        "detail": detail,
                        "note": check.get("note"),
                    }
                )
            row["passed"] = all(c["ok"] for c in row["checks"]) and row["status"] in (
                "ok",
                "target_unattainable",
            )
            if row["status"] == "target_unattainable" and case.get("expect_status") == "ok":
                row["passed"] = False
        except EngineError as exc:
            row["error"] = exc.as_dict()
            if exc.code == case.get("expect_error"):
                row["passed"] = True
        rows.append(row)
    return rows


def render(rows):
    lines = [
        f"{'case':<34} {'status':<8} {'checks':>7}  result",
        "-" * 78,
    ]
    for r in rows:
        n_ok = sum(1 for c in r["checks"] if c["ok"])
        n = len(r["checks"])
        mark = "PASS" if r["passed"] else "FAIL"
        verified = "" if r["source_verified"] else "  [UNVERIFIED SOURCE]"
        lines.append(
            f"{r['case'][:34]:<34} {mark:<8} {n_ok:>3}/{n:<3}  "
            + (r["error"]["code"] if r["error"] else r["status"] or "")
            + verified
        )
        for c in r["checks"]:
            if not c["ok"]:
                lines.append(f"    x {c['path']}: {c['detail']}")
    total = len(rows)
    passed = sum(1 for r in rows if r["passed"])
    lines.append("-" * 78)
    lines.append(f"{passed}/{total} cases passed")
    unverified = [r["case"] for r in rows if not r["source_verified"]]
    if unverified:
        lines.append(
            "cases whose primary source could not be verified in this build: "
            + ", ".join(unverified)
        )
    return "\n".join(lines)
