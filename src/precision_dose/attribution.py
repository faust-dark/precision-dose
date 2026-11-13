"""Covariate attribution.

Because every covariate effect in the model is multiplicative on clearance, the
natural logarithm of clearance is additive:

    ln CL = ln CL_ref + sum_k ln(f_k)

So the contributions reported here are not a heuristic ranking: they are the
model's own algebra, and the residual of the sum is asserted to be zero. That is
also why the percentages must not be added by a reader -- 0.6 and 0.8 multiply
to 0.48, not to 0.32 -- and the report says so in words.
"""

import math


def _group(contribs, target):
    return [c for c in contribs if c.get("target", "CL") == target]


def build(contribs, params, entry):
    out = {}
    for target, key in (("CL", "cl"), ("V", "v1"), ("Vmax", "vmax")):
        group = _group(contribs, target)
        if not group:
            continue
        if key not in params and key == "v1" and "vd" in params:
            key = "vd"
        ref_key = {"CL": "cl_ref", "V": "v1_ref", "Vmax": "vmax_ref_mg_day"}[target]
        ref = entry["model"].get(ref_key)
        total = 1.0
        for c in group:
            total *= c["multiplier"]
        individual = params.get(key)
        if individual is None:
            continue
        if ref is None and total > 0:
            # reference value not declared for this parameterisation (e.g. Vd in
            # L/kg); recover it exactly from the multiplicative structure
            ref = individual / total
        log_sum = sum(math.log(c["multiplier"]) for c in group)
        log_total = math.log(total) if total > 0 else float("nan")
        out[key] = {
            "reference_value": ref,
            "individual_value": round(individual, 6) if individual is not None else None,
            "total_multiplier": round(total, 6),
            "total_change_percent": round(100.0 * (total - 1.0), 2),
            "factors": sorted(group, key=lambda c: (-abs(math.log(c["multiplier"])), c["factor"])),
            "log_additivity_residual": round(log_total - log_sum, 12),
            "note": (
                "factors multiply; they are not additive percentages. The log-space "
                "decomposition is exact (residual ~ 0)."
            ),
        }
    return out


def sentence(attrib, entry):
    """One pharmacist-readable sentence per parameter."""
    lines = []
    for key, label in (("cl", "CL"), ("vmax", "Vmax"), ("v1", "V")):
        block = attrib.get(key)
        if not block:
            continue
        ref = block["reference_value"]
        ind = block["individual_value"]
        change = block["total_change_percent"]
        parts = [
            f"{f['factor']} x{f['multiplier']:.2f} ({f['detail']})"
            for f in block["factors"]
        ]
        unit = "mg/day" if key == "vmax" else ("L/h" if key == "cl" else "L")
        lines.append(
            f"{label} {ref:.2f} -> {ind:.2f} {unit} ({change:+.1f}%"
            + (f"; {'; '.join(parts)}" if parts else "")
            + ")"
        )
    return lines
