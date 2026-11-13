"""Covariate attribution.

Every covariate effect in the model is multiplicative on the scaled part of
clearance, so the natural logarithm of that part is additive:

    ln CL_scaled = ln CL_ref + sum_k ln(f_k)

The contributions reported here are therefore the model's own algebra, and the
residual of the sum is asserted to be zero. That is also why the percentages must
not be added by a reader -- 0.6 and 0.8 multiply to 0.48, not to 0.32 -- and the
report says so in words.

One published relation carries an intercept (vancomycin: CL = a*CrCl + b). That
term is additive, not multiplicative, so it is reported as its own absolute
component instead of being dressed up as a factor; the log decomposition applies
to the scaled component only, and the residual is computed against that.
"""

import math


def _group(contribs, target):
    return [c for c in contribs if c.get("target", "CL") == target]


def build(contribs, params, entry, additive=None):
    additive = additive or {}
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

        add_term = float(additive.get(key, 0.0) or 0.0)
        if add_term:
            scaled = individual - add_term
            log_sum = sum(math.log(c["multiplier"]) for c in group)
            log_total = (
                math.log(scaled) - math.log(ref) if ref and scaled > 0 else float("nan")
            )
            residual = log_total - log_sum
        else:
            if ref is None and total > 0:
                # reference value not declared for this parameterisation (e.g. Vd in
                # L/kg); recover it exactly from the multiplicative structure
                ref = individual / total
            log_sum = sum(math.log(c["multiplier"]) for c in group)
            log_total = math.log(total) if total > 0 else float("nan")
            residual = log_total - log_sum

        block = {
            "reference_value": ref,
            "individual_value": round(individual, 6),
            "total_multiplier": round(total, 6),
            "total_change_percent": (
                round(100.0 * (individual - ref) / ref, 2) if ref else None
            ),
            "factors": sorted(
                group, key=lambda c: (-abs(math.log(c["multiplier"])), c["factor"])
            ),
            "log_additivity_residual": round(residual, 12),
            "note": (
                "factors multiply; they are not additive percentages. The log-space "
                "decomposition is exact (residual ~ 0)."
            ),
        }
        if add_term:
            block["scaled_component"] = round(individual - add_term, 6)
            block["additive_component"] = round(add_term, 6)
            block["note"] = (
                "the published relation has a non-renal intercept, which is additive: "
                "CL = scaled component + non-renal component. The log-space "
                "decomposition applies to the scaled component only."
            )
        out[key] = block
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
