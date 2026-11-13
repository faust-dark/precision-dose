"""Engine: one patient, one drug, one answer.

Pipeline (each stage records what it used, and nothing is imputed silently):

  patient -> derived -> library entry -> individual parameters -> candidates
  -> PTA (common random numbers) -> selection -> exposure, safety, attribution
  -> rationale + output document
"""

import math

from . import attribution, clin, ddi, library, pgx, pk, report
from .errors import (
    E_LIBRARY_NO_MODEL,
    E_MISSING_FIELD,
    E_NOT_CONVERGED,
    E_TARGET_UNATTAINABLE,
    EngineError,
    W_INDUCTION_NOT_STEADY,
    W_MIC_ASSUMED,
    W_MODEL_NO_WEIGHT_V,
    W_PTA_WIDE_CI,
    W_SALT_FACTOR,
    W_TARGET_UNATTAINABLE,
)
from .matmath import DeterministicRNG, norm_ppf

DEFAULT_SEED = 20260101
DEFAULT_PTA_N = 5000


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------


def round_to(value, step):
    if not step:
        return round(value, 6)
    return round(round(value / step) * step, 6)


def _weight(patient, derived, metric, warnings):
    if metric in (None, "total"):
        return patient["weight_kg"], "total body weight"
    key = {"ibw": "ibw_kg", "abw": "abw_kg"}.get(metric)
    if key is None:
        raise EngineError(E_MISSING_FIELD, f"unknown weight metric {metric!r}")
    v = derived[key]["value"]
    if v is None:
        warnings.append(
            {
                "code": "W_MODEL_NO_WEIGHT_V",
                "message": f"{metric.upper()} requested by the drug entry but height is "
                "missing; total body weight was used instead.",
            }
        )
        return patient["weight_kg"], f"total body weight (fallback for {metric})"
    return v, metric.upper()


def _phase_multiplier(entry, patient):
    scaling = entry["model"].get("scaling") or {}
    phase = scaling.get("phase_days")
    if not phase:
        return 1.0, None
    if "days_post_transplant" not in entry["model"].get("requires", []):
        return 1.0, None
    days = patient.get("days_post_transplant")
    if days is None:
        raise EngineError(
            E_MISSING_FIELD,
            f"{entry['drug']} needs days_post_transplant: this model's clearance "
            "changes over the first post-transplant months",
            "days_post_transplant",
        )
    bounds = phase["bounds"]
    mults = phase["multipliers"]
    for i, m in enumerate(mults):
        if bounds[i] <= days < bounds[i + 1]:
            return float(m), f"day {days:.0f} post-transplant"
    return float(mults[-1]), f"day {days:.0f} post-transplant"


def _pgx_multiplier(spec, patient, warnings):
    """One gene's effect. `spec` is a [model.pgx] table."""
    gene = spec.get("gene")
    if not gene:
        return 1.0, None
    resolved = pgx.resolve(patient, gene)
    if resolved is None:
        warnings.append(
            {
                "code": "W_GENOTYPE_INCOMPLETE",
                "message": f"{gene} status is unknown; the population-average multiplier "
                "of 1.0 was used. Dosing on an unknown genotype is dosing the average.",
            }
        )
        return 1.0, None
    if spec.get("mode") == "diplotype_multipliers":
        table = spec.get("multipliers") or {}
        dip = resolved.get("diplotype")
        if dip is None:
            warnings.append(
                {
                    "code": "W_GENOTYPE_INCOMPLETE",
                    "message": f"{gene}: a phenotype was supplied without a diplotype, but "
                    "this model is parameterised on diplotype multipliers; population "
                    "average used.",
                }
            )
            return 1.0, None
        key = "/".join(sorted(dip.split("/")))
        table_norm = {"/".join(sorted(k.split("/"))): v for k, v in table.items()}
        if key not in table_norm:
            raise EngineError(
                E_MISSING_FIELD,
                f"{gene} diplotype {dip!r} has no multiplier in this entry "
                f"(known: {sorted(table)}); refusing to guess",
                "genotypes",
            )
        return float(table_norm[key]), resolved
    return float(resolved["relative_activity"]), resolved


def individualize(entry, patient, derived, interactions, warnings):
    """Return (params, contributions). Contributions are exact log factors."""
    model = entry["model"]
    scaling = model.get("scaling") or {}
    contribs = []

    def add(factor, detail, multiplier, target="CL"):
        if abs(multiplier - 1.0) < 1e-12:
            return
        contribs.append(
            {
                "factor": factor,
                "detail": detail,
                "multiplier": round(multiplier, 6),
                "delta_percent": round(100.0 * (multiplier - 1.0), 4),
                "target": target,
            }
        )

    weight_cl, wlabel_cl = _weight(
        patient, derived, scaling.get("cl_weight_metric", "total"), warnings
    )
    weight_v, wlabel_v = _weight(
        patient, derived, scaling.get("v_weight_metric", "total"), warnings
    )

    f_weight_cl = (weight_cl / 70.0) ** scaling.get("cl_weight_exponent", 0.0)
    f_weight_v = (weight_v / 70.0) ** scaling.get("v_weight_exponent", 0.0)
    add("body_weight", f"{wlabel_cl} {weight_cl:.1f} kg vs 70 kg", f_weight_cl)
    add(
        "body_weight",
        f"{wlabel_v} {weight_v:.1f} kg vs 70 kg",
        f_weight_v,
        target="V",
    )

    renal = 1.0
    rf_val = None
    if scaling.get("cl_renal_exponent"):
        metric = scaling.get("renal_metric", "crcl_cg_total")
        rf_val = clin.renal_value(patient, derived, metric)
        ref = scaling.get("renal_reference", 100.0)
        renal = (max(rf_val, 0.0) / ref) ** scaling["cl_renal_exponent"]
        add(f"renal_{metric}", f"{metric} {rf_val:.1f} vs {ref:.0f}", renal)

    f_age_cl = 1.0
    if scaling.get("cl_age_exponent"):
        ref = scaling.get("cl_age_reference", 70.0)
        f_age_cl = (patient["age_years"] / ref) ** scaling["cl_age_exponent"]
        add("age", f"age {patient['age_years']:.0f} y vs {ref:.0f} y", f_age_cl)

    f_age_v = 1.0
    if scaling.get("v1_age_exponent"):
        ref = scaling.get("v1_age_reference", 70.0)
        f_age_v = (patient["age_years"] / ref) ** scaling["v1_age_exponent"]
        add("age", f"age {patient['age_years']:.0f} y vs {ref:.0f} y", f_age_v, target="V")

    # Continuous covariates a published model may carry beyond weight/renal/age.
    # Each is a ratio to the model's own reference value; absent data refuses
    # rather than silently defaulting the factor to 1.
    f_bsa = f_alb = f_hct = f_scr = 1.0
    if scaling.get("cl_bsa_exponent"):
        bsa = derived["bsa_m2"]["value"]
        if bsa is None:
            raise EngineError(
                E_MISSING_FIELD,
                f"{entry['drug']}: this model scales clearance by body surface area, "
                "which needs height_cm",
                "height_cm",
            )
        ref = scaling.get("bsa_reference", 1.93)
        f_bsa = (bsa / ref) ** scaling["cl_bsa_exponent"]
        add("body_surface_area", f"BSA {bsa:.2f} vs {ref} m2", f_bsa)
    if scaling.get("cl_albumin_exponent"):
        alb = patient.get("albumin_g_dL")
        if alb is None:
            raise EngineError(
                E_MISSING_FIELD,
                f"{entry['drug']}: this model scales clearance by albumin",
                "albumin_g_dL",
            )
        ref = scaling.get("albumin_reference_g_L", 42.0)
        f_alb = ((alb * 10.0) / ref) ** scaling["cl_albumin_exponent"]
        add("albumin", f"albumin {alb*10.0:.0f} vs {ref:.0f} g/L", f_alb)
    if scaling.get("cl_hematocrit_exponent"):
        hct = patient.get("hematocrit_percent")
        if hct is None:
            raise EngineError(
                E_MISSING_FIELD,
                f"{entry['drug']}: this model scales clearance by haematocrit",
                "hematocrit_percent",
            )
        ref = scaling.get("hematocrit_reference", 0.34)
        f_hct = ((hct / 100.0) / ref) ** scaling["cl_hematocrit_exponent"]
        add("haematocrit", f"haematocrit {hct:.0f}% vs {ref*100:.0f}%", f_hct)
    if scaling.get("cl_creatinine_exponent"):
        ref = scaling.get("creatinine_reference_umol_L", 135.0)
        scr_umol = patient["serum_creatinine_mg_dL"] * 88.4
        f_scr = (scr_umol / ref) ** scaling["cl_creatinine_exponent"]
        add("serum_creatinine", f"SCr {scr_umol:.0f} vs {ref:.0f} umol/L", f_scr)

    f_hep = 1.0
    if scaling.get("hepatic_applies"):
        cp = derived["child_pugh"]
        if not cp["complete"]:
            warnings.append(
                {
                    "code": "W_HEPATIC_DATA_MISSING",
                    "message": "hepatic clearance scaling is enabled for this drug but "
                    "Child-Pugh is incomplete; no hepatic effect was applied.",
                }
            )
        else:
            mult = scaling.get(f"hepatic_multiplier_{cp['class']}", 1.0)
            f_hep = float(mult)
            add(
                "hepatic_child_pugh",
                f"Child-Pugh {cp['class']} (score {cp['score']})",
                f_hep,
            )

    f_phase, phase_detail = _phase_multiplier(entry, patient)
    add("post_transplant_phase", phase_detail or "", f_phase)

    pgx_specs = entry["model"].get("pgx") or []
    if isinstance(pgx_specs, dict):
        pgx_specs = [pgx_specs]
    f_pgx = 1.0
    for spec in pgx_specs:
        m, resolved = _pgx_multiplier(spec, patient, warnings)
        f_pgx *= m
        if resolved:
            add(
                "pharmacogenomics",
                f"{resolved['gene']} "
                f"{resolved.get('diplotype') or resolved['phenotype']} "
                f"({resolved['phenotype']})",
                m,
                target="Vmax" if spec.get("target") == "vmax" else "CL",
            )

    f_ddi, ddi_contrib, _ = ddi.evaluate(
        patient, model.get("ddi") or {}, interactions, warnings
    )
    for c in ddi_contrib:
        if c.get("applied"):
            add(
                "ddi",
                f"{c['perpetrator']} {c['mechanism']} on {c['pathway']}"
                + (
                    f" ([I]u {c['iu_umol_L']:.1f} umol/L)"
                    if c.get("iu_umol_L") is not None
                    else " (measured clearance ratio)"
                ),
                c["cl_multiplier"],
            )

    structure = model["structure"]
    params = {}
    if structure in ("one_compartment", "two_compartment"):
        cl_renal_part = (
            model["cl_ref"]
            * f_weight_cl
            * renal
            * f_age_cl
            * f_hep
            * f_phase
            * f_pgx
            * f_bsa
            * f_alb
            * f_hct
            * f_scr
        )
        cl_nonrenal = float(scaling.get("cl_nonrenal_l_h", 0.0) or 0.0)
        for flag, value in (scaling.get("cl_nonrenal_by_comorbidity") or {}).items():
            if flag in (patient.get("comorbidities") or []):
                cl_nonrenal = float(value)
                add("comorbidity", f"{flag} present", 1.0)  # recorded, not a multiplier
                break
        # the published renal relation may carry an intercept (Matzke 1984):
        # CL = a*CrCl + b, with b an absolute non-renal term, not a multiplier
        params["cl"] = (cl_renal_part + cl_nonrenal) * f_ddi
        params["_cl_renal_component"] = cl_renal_part * f_ddi
        params["_cl_nonrenal_component"] = cl_nonrenal * f_ddi
        params["v1"] = model["v1_ref"] * f_weight_v * f_age_v * f_hep
        if structure == "two_compartment":
            params["v2"] = model["v2_ref"] * f_weight_v
            params["q"] = model["q_ref"]
    elif structure == "michaelis_menten":
        w = weight_cl
        exp = scaling.get("vmax_weight_exponent", 0.6)
        params["vmax"] = model["vmax_ref_mg_day"] * (w / 70.0) ** exp * f_pgx
        params["km"] = model["km_mg_L"]
        params["vd"] = model["vd_ref_l_per_kg"] * weight_v
        add("body_weight", f"{wlabel_cl} {w:.1f} kg vs 70 kg", (w / 70.0) ** exp, target="Vmax")
    else:
        raise EngineError(
            E_LIBRARY_NO_MODEL,
            f"{entry['drug']} has structure={structure!r}: it is dosed by label criteria, "
            "not by a pharmacokinetic target",
        )

    if structure in ("one_compartment", "two_compartment") and scaling.get(
        "fraction_renal_clearance", 0.0
    ) == 0.0:
        if patient.get("dialysis"):
            warnings.append(
                {
                    "code": "W_DIALYSIS_ACTIVE",
                    "message": "dialysis is recorded but this drug's entry declares no "
                    "intradialytic clearance; dialysis has NOT been modelled.",
                }
            )
    return params, contribs


def build_model(entry, params):
    model = entry["model"]
    structure = model["structure"]
    f = model.get("f", 1.0)
    cf = float(model.get("concentration_factor", 1.0))
    if structure == "two_compartment":
        return pk.LinearPK(
            params["cl"], params["v1"], params.get("v2"), params.get("q"), f=f,
            conc_factor=cf,
        )
    if structure == "one_compartment":
        return pk.LinearPK(
            params["cl"],
            params["v1"],
            ka=model.get("ka_per_h"),
            f=f,
            lag_h=model.get("lag_h", 0.0),
            conc_factor=cf,
        )
    if structure == "michaelis_menten":
        return pk.MichaelisMentenPK(
            params["vmax"], params["km"], params["vd"], f=f, conc_factor=cf
        )
    raise EngineError(E_LIBRARY_NO_MODEL, f"structure {structure!r} is not doseable")


# --------------------------------------------------------------------------
# targets
# --------------------------------------------------------------------------


def select_target(entry, requested=None, patient=None):
    targets = entry["targets"]
    if not targets:
        raise EngineError(E_TARGET_UNATTAINABLE, f"{entry['drug']} declares no target")
    if requested:
        for t in targets:
            if t["id"] == requested:
                return t
        raise EngineError(
            E_TARGET_UNATTAINABLE,
            f"target {requested!r} not declared for {entry['drug']}; "
            f"valid: {[t['id'] for t in targets]}",
        )
    day = (patient or {}).get("days_post_transplant")
    day = 0.0 if day is None else day
    usable = [
        t
        for t in targets
        if t.get("applies_from_day") is None
        or (t["applies_from_day"] <= day < t.get("applies_to_day", 1e9))
    ]
    for t in usable:
        if t.get("default"):
            return t
    return usable[0] if usable else targets[0]


def candidates(entry):
    c = entry.get("constraints") or {}
    step = c.get("dose_step") or 1.0
    max_dose = c.get("max_single_dose") or 0.0
    max_daily = c.get("max_daily_dose") or math.inf
    out = []
    n = int(round(max_dose / step))
    for i in range(1, n + 1):
        dose = round(step * i, 6)
        for tau in c.get("allowed_intervals", [12]):
            if dose * 24.0 / tau <= max_daily + 1e-9:
                out.append((dose, float(tau)))
    return out


# --------------------------------------------------------------------------
# exposure evaluation
# --------------------------------------------------------------------------


def _point_times(target, interval):
    metric = target["metric"]
    if metric == "trough":
        return [interval]
    if metric == "conc_at":
        return [target.get("at_h", interval)]
    return []


def evaluate_regimen(entry, model, regimen, target, mic, dialysis=None, hd_cl_h=0.0):
    """Full evaluation of one regimen: curve, metrics, target value."""
    structure = entry["model"]["structure"]
    if structure == "michaelis_menten":
        ev = model.evaluate(regimen, grid_h=1.0)
        ev["css_avg"] = model.steady_state_concentration(
            sum(a for _, a in regimen.doses) * 24.0 / regimen.cycle_h
        )
        if ev["css_avg"] is None:
            ev["css_avg"] = math.inf
    else:
        ev = model.evaluate(regimen, dialysis=dialysis, hd_cl_h=hd_cl_h)
        points = _point_times(target, regimen.cycle_h)
        if points:
            samples = dict(model.sample(regimen, points, dialysis, hd_cl_h))
            for t, c in samples.items():
                ev.setdefault("sampled", {})[round(t, 6)] = c
    ev["target_value"] = target_value(entry, ev, target, mic)
    return ev


def target_value(entry, ev, target, mic):
    metric = target["metric"]
    if metric == "auc24":
        return ev.get("auc24")
    if metric == "auc24_over_mic":
        return ev["auc24"] / mic if mic else None
    if metric == "trough":
        if ev.get("sampled"):
            return max(ev["sampled"].values())
        return ev.get("trough")
    if metric == "cmax":
        return ev.get("cmax")
    if metric == "css_avg":
        return ev.get("css_avg")
    if metric == "css_free_avg":
        return ev.get("css_free_avg")
    if metric == "conc_at":
        if ev.get("sampled"):
            return max(ev["sampled"].values())
        return ev.get("trough")
    return None


# --------------------------------------------------------------------------
# PTA
# --------------------------------------------------------------------------


def _draw(entry, params, rng, iiv, corr):
    e_cl = rng.normal()
    e_v = rng.normal() * (1.0 - corr * corr) ** 0.5 + corr * e_cl
    out = dict(params)
    if "cl" in params:
        out["cl"] = params["cl"] * math.exp(iiv.get("cl_cv_percent", 0.0) / 100.0 * e_cl)
    if "v1" in params:
        out["v1"] = params["v1"] * math.exp(iiv.get("v1_cv_percent", 0.0) / 100.0 * e_v)
    if "v2" in params:
        out["v2"] = params["v2"] * math.exp(
            iiv.get("v2_cv_percent", 0.0) / 100.0 * rng.normal()
        )
    if "q" in params:
        out["q"] = params["q"] * math.exp(
            iiv.get("q_cv_percent", 0.0) / 100.0 * rng.normal()
        )
    if "vmax" in params:
        out["vmax"] = params["vmax"] * math.exp(
            iiv.get("vmax_cv_percent", 0.0) / 100.0 * e_cl
        )
    if "km" in params:
        out["km"] = params["km"] * math.exp(
            iiv.get("km_cv_percent", 0.0) / 100.0 * e_v
        )
    if "vd" in params:
        out["vd"] = params["vd"] * math.exp(
            iiv.get("vd_cv_percent", 0.0) / 100.0 * rng.normal()
        )
    return out


def wilson_interval(k, n, z=1.959963984540054):
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def run_pta(
    entry,
    params,
    target,
    mic,
    seed=DEFAULT_SEED,
    n=DEFAULT_PTA_N,
    dialysis=None,
    hd_cl_h=0.0,
    iiv=None,
):
    """Monte Carlo PTA with common random numbers across all candidates."""
    model_spec = entry["model"]
    iiv = iiv if iiv is not None else (model_spec.get("iiv") or {})
    corr = float(iiv.get("correlation_cl_v1", 0.0) or 0.0)
    rng = DeterministicRNG(seed)
    population = [_draw(entry, params, rng, iiv, corr) for _ in range(n)]
    structure = model_spec["structure"]
    metric = target["metric"]
    c = entry.get("constraints") or {}
    infusion = c.get("default_infusion_min", 0) / 60.0
    route = model_spec.get("route_default", "iv")
    salt = float(model_spec.get("salt_factor", 1.0))
    low, high = float(target["low"]), float(target["high"])

    # cache of unit-dose metric values per (interval) across the population
    unit_cache = {}

    def unit_values(interval):
        key = float(interval)
        if key in unit_cache:
            return unit_cache[key]
        vals = []
        if structure == "michaelis_menten":
            unit_cache[key] = None
            return None
        points = _point_times(target, key)
        for p in population:
            if metric in ("auc24", "auc24_over_mic"):
                vals.append(1.0 / p["cl"])  # AUC per unit daily dose
                continue
            m = build_model(entry, p)
            reg = pk.Regimen([(0.0, salt)], key, infusion, route)
            if points:
                s = m.sample(reg, points, dialysis, hd_cl_h)
                vals.append(max(cv for _, cv in s))
            else:
                vals.append(1.0)
        unit_cache[key] = vals
        return vals

    results = []
    for dose, tau in candidates(entry):
        daily = dose * salt * 24.0 / tau
        no_ss = 0
        if structure == "michaelis_menten":
            vals = []
            for p in population:
                model = build_model(entry, p)
                css = model.steady_state_concentration(daily)
                if css is None:
                    no_ss += 1
                    vals.append(math.inf)
                else:
                    vals.append(css)
        else:
            unit = unit_values(tau)
            if metric in ("auc24", "auc24_over_mic"):
                vals = [daily * u for u in unit]
            else:
                vals = [dose * u for u in unit]
        if metric == "auc24_over_mic" and mic:
            vals = [v / mic for v in vals]
        hits = sum(1 for v in vals if low - 1e-12 <= v <= high + 1e-12)
        sorted_vals = sorted(vals)
        median = sorted_vals[len(sorted_vals) // 2]
        lo90 = sorted_vals[int(0.05 * len(sorted_vals))]
        hi90 = sorted_vals[min(len(sorted_vals) - 1, int(0.95 * len(sorted_vals)))]
        results.append(
            {
                "dose": dose,
                "interval_h": tau,
                "daily_dose": round(daily, 6),
                "no_steady_state_fraction": (no_ss / n) if structure == "michaelis_menten" else 0.0,
                "pta": hits / n,
                "pta_ci95": list(wilson_interval(hits, n)),
                "median_target_value": median,
                "p5": lo90,
                "p95": hi90,
                "n": n,
            }
        )
    return results, {
        "n": n,
        "seed": seed,
        "rng": "mt19937+AS241",
        "common_random_numbers": True,
    }


def safety_probability(
    entry, params, safety, dose, interval, mic, seed, n, dialysis, hd_cl_h, iiv=None
):
    """P(metric exceeds threshold) for the chosen regimen, same population."""
    iiv = iiv if iiv is not None else (entry["model"].get("iiv") or {})
    corr = float(iiv.get("correlation_cl_v1", 0.0) or 0.0)
    rng = DeterministicRNG(seed)
    population = [_draw(entry, params, rng, iiv, corr) for _ in range(n)]
    structure = entry["model"]["structure"]
    metric = safety["metric"]
    c = entry.get("constraints") or {}
    infusion = c.get("default_infusion_min", 0) / 60.0
    route = entry["model"].get("route_default", "iv")
    salt = float(entry["model"].get("salt_factor", 1.0))
    daily = dose * salt * 24.0 / interval
    exceed = 0
    for p in population:
        if structure == "michaelis_menten":
            model = build_model(entry, p)
            v = model.steady_state_concentration(daily)
            v = math.inf if v is None else v
        elif metric == "auc24":
            v = daily / p["cl"]
        else:
            model = build_model(entry, p)
            reg = pk.Regimen([(0.0, dose * salt)], interval, infusion, route)
            pts = _point_times({"metric": metric}, interval) or [interval]
            v = max(cv for _, cv in model.sample(reg, pts, dialysis, hd_cl_h))
        if metric == "auc24_over_mic" and mic:
            v = v / mic
        if v > safety["threshold"]:
            exceed += 1
    return exceed / n


# --------------------------------------------------------------------------
# main entry point
# --------------------------------------------------------------------------


def run_dose(
    patient_doc,
    drug,
    target_id=None,
    seed=DEFAULT_SEED,
    pta_n=DEFAULT_PTA_N,
    library_dir=None,
    curve_points=200,
    renal_metric=None,
    weight_metric=None,
    param_override=None,
    iiv_override=None,
    posterior=None,
):
    patient = clin.validate_patient(patient_doc)
    derived, warnings = clin.derive(patient)
    lib, interactions, problems = library.load_library(library_dir)
    entry = library.get_drug(lib, drug)

    if entry["decision_basis"] == "label_criteria":
        return _run_label(entry, patient, derived, warnings, interactions, lib)

    _require_fields(entry, patient, derived)
    target = select_target(entry, target_id, patient)
    mic = patient.get("mic_mg_L")
    if target["metric"] == "auc24_over_mic" and not mic:
        mic = 1.0
        warnings.append(
            {
                "code": "W_MIC_ASSUMED",
                "message": "no MIC supplied: MIC 1.0 mg/L assumed, as in the 2020 "
                "vancomycin guideline. Supply mic_mg_L to remove this assumption.",
            }
        )

    params, contribs = individualize(entry, patient, derived, interactions, warnings)

    iiv_used = dict(entry["model"].get("iiv") or {})
    if param_override:
        for k, v in param_override.items():
            if k in params:
                params[k] = float(v)
        # The attribution must multiply out to the parameter actually used, so a
        # Bayesian shift enters as its own exact factor rather than being hidden.
        for tgt, key, ref_key in (
            ("CL", "cl", "cl_ref"),
            ("V", "v1", "v1_ref"),
            ("Vmax", "vmax", "vmax_ref_mg_day"),
        ):
            if key not in param_override or key not in params:
                continue
            prod = 1.0
            for c in contribs:
                if c.get("target", "CL") == tgt:
                    prod *= c["multiplier"]
            ref = entry["model"].get(ref_key)
            if not ref or prod <= 0:
                continue
            extra = params[key] / (ref * prod)
            if abs(extra - 1.0) > 1e-9:
                contribs.append(
                    {
                        "factor": "tdm_posterior",
                        "detail": "MAP estimate from the measured levels",
                        "multiplier": round(extra, 6),
                        "delta_percent": round(100.0 * (extra - 1.0), 4),
                        "target": tgt,
                    }
                )
    if iiv_override:
        for k, v in iiv_override.items():
            base = float(iiv_used.get(k) or 0.0)
            iiv_used[k] = (base * base + float(v) * float(v)) ** 0.5

    dialysis = patient.get("dialysis")
    hd_cl_h = 0.0
    dial = entry["model"].get("dialysis") or {}
    if dialysis:
        if dialysis["modality"] != "hd":
            raise EngineError(
                E_LIBRARY_NO_MODEL,
                f"dialysis modality {dialysis['modality']!r} is not supported for "
                f"{entry['drug']} in v1",
            )
        if dial.get("hd_additional_cl_l_h"):
            hd_cl_h = float(dial["hd_additional_cl_l_h"])
            warnings.append(
                {
                    "code": "W_DIALYSIS_ACTIVE",
                    "message": f"intermittent haemodialysis modelled: +{hd_cl_h:.1f} L/h "
                    f"intradialytic clearance on {dialysis['sessions_per_week']} sessions "
                    f"per week of {dialysis['hours_per_session']:.0f} h.",
                }
            )
        else:
            frac = float((entry["model"].get("scaling") or {}).get(
                "fraction_renal_clearance", 0.0))
            if frac >= 0.5:
                raise EngineError(
                    E_LIBRARY_NO_MODEL,
                    f"{entry['drug']} is {frac*100:.0f}% renally cleared and its entry "
                    "declares no intradialytic clearance: dosing a dialysis patient "
                    "without it would under-dose by a large margin, so the engine "
                    "refuses instead of guessing.",
                )
            warnings.append(
                {
                    "code": "W_DIALYSIS_ACTIVE",
                    "message": "dialysis is recorded but this drug's entry declares no "
                    "intradialytic clearance; dialysis has NOT been modelled.",
                }
            )

    pta_rows, pta_meta = run_pta(
        entry, params, target, mic, seed=seed, n=pta_n,
        dialysis=dialysis, hd_cl_h=hd_cl_h, iiv=iiv_used,
    )
    chosen = _select(entry, pta_rows, target, warnings)

    interval = chosen["interval_h"]
    infusion = (entry.get("constraints") or {}).get("default_infusion_min", 0) / 60.0
    route = entry["model"].get("route_default", "iv")
    salt = float(entry["model"].get("salt_factor", 1.0))
    regimen = pk.Regimen([(0.0, chosen["dose"] * salt)], interval, infusion, route)
    if salt != 1.0:
        warnings.append(
            {
                "code": W_SALT_FACTOR,
                "message": f"dose is prescribed as the salt and converted to base at "
                f"x{salt:g} (e.g. phenytoin sodium 100 mg = {100*salt:.0f} mg acid) before "
                "the model is applied.",
            }
        )
    model = build_model(entry, params)
    ev = evaluate_regimen(entry, model, regimen, target, mic, dialysis, hd_cl_h)

    if dialysis and entry["model"]["structure"] != "michaelis_menten":
        weekly = pk.Regimen(
            [
                (
                    float(dialysis["first_session_h"])
                    + i * float(dialysis["day_interval_h"])
                    + float(dialysis["hours_per_session"]),
                    chosen["dose"] * salt,
                )
                for i in range(int(dialysis["sessions_per_week"]))
            ],
            168.0,
            infusion,
            route,
        )
        ev_week = model.evaluate(weekly, dialysis=dialysis, hd_cl_h=hd_cl_h)
        ev["dialysis_week"] = _summarise_dialysis(ev_week, dialysis)

    safety_rows = []
    for s in entry["safety"]:
        if s["metric"] == "none":
            safety_rows.append(
                {
                    "id": s["id"],
                    "metric": "none",
                    "threshold": None,
                    "predicted": None,
                    "p_exceed": None,
                    "note": s.get("note"),
                }
            )
            continue
        p_ex = safety_probability(
            entry, params, s, chosen["dose"], interval, mic, seed, min(pta_n, 2000),
            dialysis, hd_cl_h, iiv=iiv_used,
        )
        predicted = (
            ev["target_value"] if s["metric"] == target["metric"] else ev.get(s["metric"])
        )
        if s["metric"] == "auc24":
            predicted = ev["auc24"]
        elif s["metric"] == "trough":
            predicted = ev.get("trough")
        elif s["metric"] == "conc_at":
            pts = model.sample(regimen, [s["at_h"]], dialysis, hd_cl_h)
            predicted = pts[0][1]
        elif s["metric"] == "css_avg":
            predicted = ev.get("css_avg")
        safety_rows.append(
            {
                "id": s["id"],
                "metric": s["metric"],
                "threshold": s["threshold"],
                "predicted": predicted,
                "margin_ratio": (
                    s["threshold"] / predicted if predicted and predicted > 0 else None
                ),
                "p_exceed": p_ex,
                "direction": s.get("direction"),
                "source": s.get("source"),
            }
        )

    deduped = []
    seen = set()
    for w in warnings:
        key = (w.get("code"), w.get("message"))
        if key in seen:
            continue
        seen.add(key)
        deduped.append(w)
    warnings[:] = deduped

    loading = _loading(entry, patient, derived, warnings)
    attrib = attribution.build(
        contribs, params, entry, additive={"cl": params.get("_cl_nonrenal_component", 0.0)}
    )

    result = {
        "schema_version": "1.0",
        "status": "ok",
        "engine": {
            "name": "precision-dose",
            "version": _version(),
            "seed": seed,
            "pta_n": pta_n,
            "rng": "mt19937+AS241",
            "library_hash": entry["content_hash"],
            "entry_version": entry["entry_version"],
        },
        "question": {
            "drug": entry["drug"],
            "display_name": entry["display_name"],
            "target_id": target["id"],
            "target_metric": target["metric"],
            "indication": patient.get("indication"),
            "mic_mg_L": mic,
            "decision_basis": "pk_target",
        },
        "patient": {
            "patient_id": patient.get("patient_id"),
            "input": {
                "age_years": patient["age_years"],
                "sex": patient["sex"],
                "weight_kg": patient["weight_kg"],
                "height_cm": patient.get("height_cm"),
                "serum_creatinine_mg_dL": patient["serum_creatinine_mg_dL"],
                "albumin_g_dL": patient.get("albumin_g_dL"),
                "days_post_transplant": patient.get("days_post_transplant"),
                "genotypes": patient.get("genotypes"),
                "co_medications": [m["name"] for m in patient.get("co_medications", [])],
                "dialysis": bool(dialysis),
            },
            "derived": {k: v for k, v in derived.items()},
        },
        "library_entry": {
            "drug": entry["drug"],
            "entry_version": entry["entry_version"],
            "content_hash": entry["content_hash"],
            "citations": entry["citations"],
            "known_gaps": entry["provenance"].get("known_gaps", []),
        },
        "parameters": {
            "structure": entry["model"]["structure"],
            "values": {k: round(v, 6) for k, v in params.items() if not k.startswith("_")},
            "units": {k: _param_unit(k, entry) for k in params if not k.startswith("_")},
            "individualized": bool(contribs),
            "t_half_h": _t_half(model, regimen, entry),
            "t_half_basis": (
                "effective over the dosing interval"
                if entry["model"]["structure"] == "two_compartment"
                else "terminal"
            ),
            "posterior": posterior,
        },
        "attribution": attrib,
        "regimen": {
            "loading": loading,
            "maintenance": {
                "dose": chosen["dose"],
                "unit": entry["model"]["amount_unit"],
                "dose_mass_mg": (
                    round(chosen["dose"] * entry["model"]["mass_per_unit_mg"], 1)
                    if entry["model"].get("mass_per_unit_mg")
                    else None
                ),
                "interval_h": interval,
                "infusion_min": int(round(infusion * 60)) or None,
                "route": route,
            },
            "decision_basis": "pk_target",
            "candidates_considered": len(pta_rows),
        },
        "exposure": {
            "auc24": round(ev["auc24"], 4) if ev.get("auc24") is not None else None,
            "auc24_over_mic": (
                round(ev["auc24"] / mic, 4) if mic and ev.get("auc24") else None
            ),
            "cmax": round(ev["cmax"], 4),
            "cmin": round(ev["cmin"], 4),
            "trough": round(ev.get("trough", ev["cmin"]), 4),
            "tmax_h": round(ev["tmax_h"], 3),
            "css_avg": round(ev["css_avg"], 4) if ev.get("css_avg") is not None else None,
            "unit": entry["model"]["concentration_unit"],
            "matrix": entry["model"]["concentration_matrix"],
            "dialysis_week": ev.get("dialysis_week"),
        },
        "pta": {
            "target_low": target["low"],
            "target_high": target["high"],
            "target_unit": target.get("unit"),
            "probability": round(chosen["pta"], 6),
            "ci95": [round(x, 6) for x in chosen["pta_ci95"]],
            "n": pta_meta["n"],
            "seed": seed,
            "median_target_value": round(chosen["median_target_value"], 4),
            "p5": round(chosen["p5"], 4),
            "p95": round(chosen["p95"], 4),
            "table": _pta_table(pta_rows, target),
        },
        "safety": safety_rows,
        "uncertainty": {
            "iiv": entry["model"].get("iiv"),
            "target_pi90": [round(chosen["p5"], 4), round(chosen["p95"], 4)],
            "pta_ci_width": round(chosen["pta_ci95"][1] - chosen["pta_ci95"][0], 6),
        },
        "warnings": warnings,
        "curve": _curve(ev, curve_points, entry),
        "rationale": "",
        "disclaimer": (
            "Computational decision support for a pharmacist. Not a prescribing "
            "instruction. Every recommendation requires professional review against "
            "the patient's full clinical picture."
        ),
    }
    max_daily = max(r["daily_dose"] for r in pta_rows) if pta_rows else None
    top = [r for r in pta_rows if r["daily_dose"] == max_daily]
    unreachable = bool(top) and max(r["median_target_value"] for r in top) < target["low"]
    if unreachable:
        result["status"] = "target_unattainable"
        result["warnings"].append(
            {
                "code": W_TARGET_UNATTAINABLE,
                "message": f"even the highest permitted dose (median {top[0]['median_target_value']:.1f}) "
                f"does not reach the target low bound {target['low']:g}; the target is out of "
                "reach within this entry's dosing constraints.",
            }
        )
    if chosen.get("no_steady_state_fraction", 0) > 0.05:
        result["warnings"].append(
            {
                "code": "W_NONLINEAR_NO_STEADY_STATE",
                "message": f"{chosen['no_steady_state_fraction']*100:.0f}% of simulated "
                "patients have no steady state at this dose (input exceeds Vmax): "
                "capacity-limited drugs can accumulate without bound, which is exactly "
                "why phenytoin needs measured levels rather than a calculated dose.",
            }
        )
    if chosen["pta"] < 0.8:
        result["warnings"].append(
            {
                "code": "W_LOW_PTA",
                "message": f"no fixed regimen attains {target['id']} in 80% of patients "
                f"(best {chosen['pta']*100:.0f}%): the target window is narrow relative to "
                "between-patient variability. This is the quantitative case for "
                "therapeutic drug monitoring rather than a reason to widen the target.",
            }
        )
    result["rationale"] = report.rationale(result, entry, patient, derived, attrib)
    result["trace"] = report.trace(result, entry, patient, derived, params, ev)
    return result


def _t_half(model, regimen, entry):
    """Effective half-life for multi-compartment fits, terminal otherwise."""
    if entry["model"]["structure"] == "two_compartment" and hasattr(
        model, "effective_half_life"
    ):
        val = model.effective_half_life(regimen)
        return round(val, 4) if val else None
    if hasattr(model, "half_life"):
        return round(model.half_life(), 4)
    return None


def _require_fields(entry, patient, derived):
    """Refuse early, with a named field, when the model needs data that is absent."""
    for field in entry["model"].get("requires", []):
        if field == "bsa_m2":
            if derived["bsa_m2"]["value"] is None:
                raise EngineError(
                    E_MISSING_FIELD,
                    f"{entry['drug']} needs body surface area, which needs height_cm",
                    "height_cm",
                )
            continue
        if patient.get(field) is None:
            raise EngineError(
                E_MISSING_FIELD,
                f"{entry['drug']}: this pharmacokinetic model requires {field}; it has no "
                "population default, so the engine will not substitute one",
                field,
            )


def _select(entry, rows, target, warnings):
    feasible = [r for r in rows if r.get("safety_ok") is not False] if rows else []
    if not feasible:
        feasible = rows
    # Tie-break, in order: highest PTA; highest PST; lowest total daily dose; then
    # the interval closest to the one the guideline actually uses (entries declare
    # preferred_interval_h), then the smaller dose. Regimens with the same daily
    # dose attain the target equally often, so this only decides which of several
    # equivalent exposures is the more conventional one to write on a chart.
    preferred = (entry.get("constraints") or {}).get("preferred_interval_h")
    feasible.sort(
        key=lambda r: (
            -round(r["pta"], 9),
            -round(r.get("pst", r["pta"]), 9),
            r["daily_dose"],
            abs(r["interval_h"] - preferred) if preferred else r["interval_h"],
            r["dose"],
        )
    )
    best = feasible[0]
    if best["pta_ci95"][1] - best["pta_ci95"][0] > 0.10:
        warnings.append(
            {
                "code": W_PTA_WIDE_CI,
                "message": "the PTA 95% interval is wider than 0.10: increase --pta-n for a "
                "sharper estimate.",
            }
        )
    return best


def _pta_table(rows, target):
    out = []
    for r in sorted(rows, key=lambda r: (-r["pta"], r["dose"], r["interval_h"]))[:12]:
        out.append(
            {
                "dose": r["dose"],
                "interval_h": r["interval_h"],
                "pta": round(r["pta"], 4),
                "median_target": round(r["median_target_value"], 3),
            }
        )
    return out


def _curve(ev, points, entry):
    curve = ev.get("curve") or []
    if not curve:
        return {"time_h": [], "concentration": [], "unit": entry["model"]["concentration_unit"]}
    stride = max(1, len(curve) // max(points, 2))
    sampled = curve[::stride]
    if sampled[-1] != curve[-1]:
        sampled.append(curve[-1])
    return {
        "time_h": [round(t, 3) for t, _ in sampled],
        "concentration": [round(c, 4) for _, c in sampled],
        "unit": entry["model"]["concentration_unit"],
        "matrix": entry["model"]["concentration_matrix"],
    }


def _loading(entry, patient, derived, warnings):
    c = entry.get("constraints") or {}
    per_kg = c.get("loading_per_kg") or c.get("loading_mg_per_kg") or 0.0
    if per_kg <= 0:
        return None
    w, label = _weight(patient, derived, c.get("loading_weight_metric", "total"), warnings)
    raw = per_kg * w
    cap = c.get("loading_max") or raw
    dose = round_to(min(raw, cap), c.get("dose_step") or 250.0)
    return {
        "dose": dose,
        "unit": entry["model"]["amount_unit"],
        "per_kg": per_kg,
        "per_kg_unit": entry["model"]["amount_unit"],
        "weight_kg": round(w, 2),
        "weight_metric": label,
        "infusion_min": c.get("loading_infusion_min"),
        "capped": raw > cap,
        "note": "loading dose sets time to target; it does not change steady-state exposure",
    }


def _summarise_dialysis(ev, dialysis):
    curve = ev["curve"]
    out = {"sessions": []}
    dur = float(dialysis["hours_per_session"])
    for i in range(int(dialysis["sessions_per_week"])):
        start = float(dialysis["first_session_h"]) + i * float(dialysis["day_interval_h"])
        pre = _at(curve, start)
        post = _at(curve, start + dur)
        out["sessions"].append(
            {
                "session": i + 1,
                "start_h": round(start, 2),
                "pre_hd": round(pre, 4),
                "post_hd": round(post, 4),
                "removal_percent": round(100.0 * (pre - post) / pre, 2) if pre else None,
                "rebound_4h": round(_at(curve, start + dur + 4.0) - post, 4),
            }
        )
    out["auc24_mean"] = round(ev["auc24"], 4)
    out["note"] = "post-HD rebound is a redistribution phenomenon and is reproduced by the two-compartment model, not assumed."
    return out


def _at(curve, t):
    best = curve[0]
    for pt in curve:
        if pt[0] <= t + 1e-9:
            best = pt
        else:
            break
    return best[1]


def _param_unit(name, entry):
    m = entry["model"]
    if name in ("cl",):
        return "L/h"
    if name in ("v1", "v2", "vd"):
        return "L"
    if name == "q":
        return "L/h"
    if name == "vmax":
        return "mg/day"
    if name == "km":
        return m["concentration_unit"]
    return ""


def _version():
    from . import __version__

    return __version__


# --------------------------------------------------------------------------
# label-criteria pathway
# --------------------------------------------------------------------------


def _run_label(entry, patient, derived, warnings, interactions, lib):
    fired = []
    evaluated = []
    for rule in entry.get("label_rules", []):
        kind = rule["type"]
        if kind == "count_criteria":
            met = []
            for crit in rule.get("criteria", []):
                v = patient.get(crit["field"])
                op = crit["op"]
                ok = (
                    v is not None
                    and (
                        (op == "ge" and v >= crit["value"])
                        or (op == "le" and v <= crit["value"])
                        or (op == "gt" and v > crit["value"])
                        or (op == "lt" and v < crit["value"])
                    )
                )
                met.append(
                    {
                        "criterion": f"{crit['field']} {op} {crit['value']}",
                        "value": v,
                        "met": bool(ok),
                    }
                )
            count = sum(1 for m in met if m["met"])
            hit = count >= rule["threshold"]
            evaluated.append({**rule, "criteria_met": met, "count": count, "fired": hit})
            if hit:
                fired.append(rule)
        elif kind == "co_medication_class":
            names = [
                m["name"]
                for m in patient.get("co_medications", [])
                if rule["class"]
                in (interactions.get("perpetrators", {}).get(m["name"], {}).get("classes") or [])
            ]
            hit = bool(names)
            evaluated.append({**rule, "matched_medications": names, "fired": hit})
            if hit:
                fired.append(rule)
        elif kind == "child_pugh_class":
            cp = derived["child_pugh"]
            hit = cp.get("class") == rule["class"]
            evaluated.append({**rule, "child_pugh": cp, "fired": hit})
            if hit:
                fired.append(rule)
        elif kind == "default":
            evaluated.append({**rule, "fired": not fired})
            if not fired:
                fired.append(rule)

    severity = {"not_recommended": 3, "avoid": 2, "reduce": 1, "full_dose": 0}
    fired.sort(key=lambda r: (-severity.get(r.get("action"), 0), r.get("dose") or 0))
    chosen = fired[0] if fired else None
    if chosen is None:
        raise EngineError(E_TARGET_UNATTAINABLE, "no label rule applied to this patient")

    warnings.append(
        {
            "code": "W_LABEL_CRITERIA",
            "message": f"{entry['drug']}: the dose is determined by label criteria "
            f"(rule {chosen['id']}), not by a pharmacokinetic target. No concentration "
            "prediction is produced.",
        }
    )
    regimen = {
        "dose": chosen.get("dose"),
        "unit": entry["model"]["amount_unit"],
        "interval_h": chosen.get("interval_h"),
        "route": chosen.get("route"),
        "label": chosen.get("label"),
        "rule_id": chosen["id"],
        "action": chosen.get("action"),
    }
    return {
        "schema_version": "1.0",
        "status": "ok" if chosen.get("action") != "not_recommended" else "not_recommended",
        "engine": {
            "name": "precision-dose",
            "version": _version(),
            "library_hash": entry["content_hash"],
            "entry_version": entry["entry_version"],
        },
        "question": {
            "drug": entry["drug"],
            "display_name": entry["display_name"],
            "decision_basis": "label_criteria",
            "target_id": None,
            "indication": patient.get("indication"),
        },
        "patient": {
            "patient_id": patient.get("patient_id"),
            "derived": derived,
        },
        "library_entry": {
            "drug": entry["drug"],
            "entry_version": entry["entry_version"],
            "content_hash": entry["content_hash"],
            "citations": entry["citations"],
            "known_gaps": entry["provenance"].get("known_gaps", []),
        },
        "regimen": {**regimen, "decision_basis": "label_criteria"},
        "label_evaluation": [
            {
                "rule": e["id"],
                "type": e["type"],
                "fired": e["fired"],
                "criteria_met": e.get("criteria_met"),
                "count": e.get("count"),
                "matched_medications": e.get("matched_medications"),
            }
            for e in evaluated
        ],
        "exposure": None,
        "pta": None,
        "safety": [
            {"id": s["id"], "metric": s["metric"], "note": s.get("note")}
            for s in entry["safety"]
        ],
        "warnings": warnings,
        "rationale": report.label_rationale(entry, patient, derived, regimen, evaluated),
        "disclaimer": (
            "Computational decision support for a pharmacist. Not a prescribing "
            "instruction; verify against the current product label."
        ),
    }
