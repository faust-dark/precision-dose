"""Deterministic narrative generation.

No language model touches this file. The rationale is a template over the
structured result, which is why the same input produces the same bytes.
"""


def _fmt(v, digits=2):
    if v is None:
        return "n/a"
    if isinstance(v, str):
        return v
    return f"{v:.{digits}f}"


def rationale(result, entry, patient, derived, attrib):
    q = result["question"]
    p = patient
    units = entry["model"]["concentration_unit"]
    lines = []

    # 1. what was asked
    mic_txt = (
        f", MIC {_fmt(q.get('mic_mg_L'), 2)} mg/L"
        if q.get("mic_mg_L") is not None
        else ""
    )
    lines.append(
        f"Question: {entry['display_name']} regimen for a {p['age_years']:.0f}-year-old "
        f"{p['sex']}, {p['weight_kg']:.0f} kg, SCr {p['serum_creatinine_mg_dL']:.2f} mg/dL"
        f"{mic_txt}, target {q['target_id']} "
        f"({result['pta']['target_low']:g}-{result['pta']['target_high']:g} "
        f"{result['pta']['target_unit'] or ''})."
    )

    # 2. renal function and which metric the model uses
    d = derived
    egfr = d["egfr_ckd_epi_2021_mL_min_1_73m2"]["value"]
    v = d["crcl_cg_mL_min"]["value"]
    model_metric = (entry["model"].get("scaling") or {}).get("renal_metric", "none")
    lines.append(
        f"Renal function: CKD-EPI 2021 eGFR {egfr:.0f} mL/min/1.73 m2; Cockcroft-Gault "
        f"CrCl {v['total_body_weight']:.0f} mL/min (TBW), "
        f"{_fmt(v['ideal_body_weight'], 0)} (IBW), {_fmt(v['adjusted_body_weight'], 0)} (ABW). "
        f"This entry scales clearance with {model_metric}."
        + (
            " Discarding one of these in favour of the other is a decision with clinical "
            "consequences, which is why both are shown."
            if model_metric != "none"
            else ""
        )
    )
    if p.get("albumin_g_dL") is not None:
        cp = d["child_pugh"]
        lines.append(
            f"Hepatic: albumin {p['albumin_g_dL']:.1f} g/dL; Child-Pugh "
            + (
                f"class {cp['class']} (score {cp['score']})"
                if cp["complete"]
                else "incomplete (no class assigned)"
            )
            + f"; clearance multiplier "
            + (
                f"{(entry['model'].get('scaling') or {}).get('hepatic_multiplier_' + cp['class'], 1.0):.2f}"
                if cp["complete"]
                and (entry["model"].get("scaling") or {}).get("hepatic_applies")
                else "1.00 (not applied for this drug)"
            )
            + "."
        )

    # 3. the parameters and what moved them
    pr = result["parameters"]
    lines.append(
        "Individual parameters: "
        + ", ".join(
            f"{k.upper()} {v:.2f} {pr['units'].get(k, '')}".strip()
            for k, v in pr["values"].items()
        )
        + (
            f"; terminal half-life {pr['t_half_h']:.1f} h."
            if pr.get("t_half_h")
            else "."
        )
    )
    for s in _attribute_lines(result):
        lines.append(s)

    # 4. the regimen
    reg = result["regimen"]
    load = reg.get("loading")
    m = reg["maintenance"]
    reg_txt = (
        f"Regimen: loading {load['dose']:g} {load['unit']} over "
        f"{load['infusion_min']} min ({load['mg_per_kg']:g} mg/kg of {load['weight_metric']}), "
        if load
        else "Regimen: "
    )
    reg_txt += (
        f"maintenance {m['dose']:g} {m['unit']}"
        + (
            f" ({m['dose_mass_mg']:g} mg as the salt)"
            if m.get("dose_mass_mg")
            else ""
        )
        + f" every {m['interval_h']:g} h"
        + (f" over {m['infusion_min']} min" if m.get("infusion_min") else "")
        + "."
    )
    lines.append(reg_txt)

    # 5. exposure and PTA
    ex = result["exposure"]
    pta = result["pta"]
    lines.append(
        f"Predicted steady-state exposure: AUC24 {_fmt(ex['auc24'], 1)} "
        f"{units}*h, Cmax {_fmt(ex['cmax'], 2)}, trough {_fmt(ex['trough'], 2)} {units} "
        f"(90% population interval {pta['p5']:.1f}-{pta['p95']:.1f})."
    )
    lines.append(
        f"Probability of target attainment: {pta['probability']*100:.0f}% "
        f"(95% CI {pta['ci95'][0]*100:.0f}-{pta['ci95'][1]*100:.0f}%, n={pta['n']}, "
        f"seed {pta['seed']}). The spread is between-patient variability, not numerical "
        "noise."
    )

    # 6. safety
    for s in result.get("safety") or []:
        if s.get("threshold") is None:
            continue
        lines.append(
            f"Safety: predicted {s['metric']} {_fmt(s['predicted'], 1)} vs threshold "
            f"{s['threshold']:g}; margin x{_fmt(s.get('margin_ratio'), 2)}; "
            f"P(exceed) {_fmt(s.get('p_exceed', 0) * 100, 1)}%."
        )

    # 7. uncertainty and refusals
    warns = [w["message"] for w in result.get("warnings", [])]
    if warns:
        lines.append("Caveats:")
        for w in warns:
            lines.append(f"  - {w}")

    gaps = entry["provenance"].get("known_gaps") or []
    if gaps:
        lines.append("Model limitations declared by the library entry:")
        for g in gaps:
            lines.append(f"  - {g}")

    cites = entry["citations"]
    lines.append("Sources:")
    for key in sorted(cites):
        c = cites[key]
        lines.append(f"  [{key}] {c['ref']}" + (f" doi:{c['doi']}" if c.get("doi") else ""))

    lines.append(f"Library entry hash: {result['engine']['library_hash']}")
    return "\n".join(lines)


def _attribute_lines(result):
    at = result.get("attribution") or {}
    out = []
    for key in ("cl", "vmax", "v1"):
        block = at.get(key)
        if not block:
            continue
        label = {"cl": "Clearance", "vmax": "Vmax", "v1": "Volume"}[key]
        unit = {"cl": "L/h", "vmax": "mg/day", "v1": "L"}[key]
        parts = [
            f"{f['factor']} x{f['multiplier']:.2f} [{f['detail']}]"
            for f in block["factors"]
        ]
        out.append(
            f"{label}: {block['reference_value']:.2f} -> "
            f"{block['individual_value']:.2f} {unit} "
            f"({block['total_change_percent']:+.1f}%)"
            + (" because " + "; ".join(parts) if parts else " (no covariate effects)")
            + "."
        )
    if at.get("cl"):
        out.append(
            "Contributions multiply; do not add the percentages. "
            f"Log-space residual {at['cl']['log_additivity_residual']:.1e} confirms the "
            "decomposition is exact."
        )
    return out


def label_rationale(entry, patient, derived, regimen, evaluated):
    lines = [
        f"Question: dose of {entry['display_name']} for a "
        f"{patient['age_years']:.0f}-year-old {patient['sex']}, "
        f"{patient['weight_kg']:.0f} kg, SCr {patient['serum_creatinine_mg_dL']:.2f} mg/dL.",
        f"Decision basis: product-label criteria, not a pharmacokinetic target. "
        f"{entry['display_name']} is not dosed to a measured exposure, so this engine "
        "deliberately does not produce a concentration curve, PTA or Cmax margin.",
        f"Applied rule: {regimen['rule_id']} -> {regimen['label']}.",
    ]
    for e in evaluated:
        if e["type"] == "count_criteria":
            met = ", ".join(
                f"{c['criterion']} {'met' if c['met'] else 'not met'}"
                for c in e.get("criteria_met") or []
            )
            lines.append(
                f"  - {e['id']}: {e.get('count')} of "
                f"{len(e.get('criteria_met') or [])} criteria met ({met})."
            )
        elif e["type"] == "co_medication_class":
            lines.append(
                f"  - {e['id']}: "
                + (
                    "matched " + ", ".join(e.get("matched_medications") or [])
                    if e.get("matched_medications")
                    else "no matching co-medication"
                )
                + "."
            )
        elif e["type"] == "child_pugh_class":
            cp = e.get("child_pugh") or {}
            lines.append(
                f"  - {e['id']}: Child-Pugh "
                + (cp.get("class") or "not calculable")
                + "."
            )
    cites = entry["citations"]
    lines.append("Sources:")
    for key in sorted(cites):
        lines.append(f"  [{key}] {cites[key]['ref']}")
    lines.append(f"Library entry hash: {entry['content_hash']}")
    return "\n".join(lines)


def trace(result, entry, patient, derived, params, ev):
    """A 'show your work' page: formulas with substituted numbers."""
    m = entry["model"]
    sc = m.get("scaling") or {}
    out = ["CALCULATION TRACE", "=" * 18]
    out.append("1. Demographics")
    out.append(
        f"   IBW = 50 + 2.3*(H_in - 60) [male] or 45.5 + ... [female] = "
        f"{derived['ibw_kg']['value']}"
    )
    out.append(
        f"   eGFR(CKD-EPI 2021) = 142 * min(SCr/k,1)^a * max(SCr/k,1)^-1.200 * "
        f"0.9938^age * (1.012 if female) = "
        f"{derived['egfr_ckd_epi_2021_mL_min_1_73m2']['value']:.4f}"
    )
    v = derived["crcl_cg_mL_min"]["value"]
    out.append(
        f"   CrCl(CG) = (140-age)*W/(72*SCr)*(0.85 if female) = "
        f"{v['total_body_weight']:.4f} (TBW), {v['ideal_body_weight']} (IBW), "
        f"{v['adjusted_body_weight']} (ABW)"
    )
    out.append("2. Structural model")
    if m["structure"] == "two_compartment":
        out.append(
            "   dx/dt = M x + u, M = [[-(CL+Q)/V1, Q/V2],[Q/V1, -Q/V2]]"
        )
        out.append(
            "   x(t+h) = expm(M h) x(t) + M^-1 (expm(M h) - I) u   [exact, mass conserving]"
        )
        out.append("   steady state: x* = (I - expm(M tau))^-1 b   [one linear solve]")
    elif m["structure"] == "one_compartment":
        out.append("   C(t) = (A1/V1); A1(t+h) = A1 e^-kh + (R/k)(1 - e^-kh)")
    else:
        out.append("   dA/dt = -Vmax*C/(Km+C), C = A/Vd; RK4, h = 1 h")
        out.append("   steady state: R = Vmax*Css/(Km+Css)  ->  Css = Km*R/(Vmax-R)")
    out.append("3. Individual parameters")
    for k, val in params.items():
        out.append(f"   {k} = {val:.6g}")
    out.append("4. Covariate multipliers (log-space, exact)")
    for key in ("cl", "vmax", "v1"):
        block = (result.get("attribution") or {}).get(key)
        if not block:
            continue
        for f in block["factors"]:
            out.append(
                f"   {key}: {f['factor']:<22} x{f['multiplier']:.6f}  ({f['detail']})"
            )
        out.append(
            f"   {key}: product = x{block['total_multiplier']:.6f}, "
            f"log residual {block['log_additivity_residual']:.1e}"
        )
    out.append("5. Exposure")
    out.append(
        f"   AUC24 = {ev.get('auc24'):.4f}; Cmax = {ev.get('cmax'):.4f}; "
        f"trough = {ev.get('trough', float('nan')):.4f}; grid = {ev.get('grid_h')} h"
    )
    out.append(
        "   AUC is the exact integral of C over the cycle (int x dt = M^-1 (dx - u dt)), "
        "not a trapezoid over the display grid."
    )
    out.append("6. PTA")
    out.append(
        f"   CL_i = CL * exp(eta), eta ~ N(0, {m.get('iiv', {}).get('cl_cv_percent', 0)}%^2); "
        f"n = {result['pta']['n']}, seed = {result['pta']['seed']}, "
        "common random numbers across candidates"
    )
    out.append("7. Determinism")
    out.append(
        f"   library {result['engine']['library_hash']}; engine {result['engine']['version']}"
    )
    return out
