"""Patient validation, body size, renal function, hepatic status.

Clinical distinctions that are deliberately *not* glossed over:

- eGFR (CKD-EPI 2021) is reported for staging and CKD context. Drug labels and
  most dosing guidelines were built on Cockcroft-Gault creatinine clearance, so
  CrCl is computed too, with three weight metrics. Which one drives clearance
  is declared by the drug entry, never guessed here.
- Child-Pugh classifies hepatic impairment on a 5-15 point scale and is used
  only where a drug entry supplies a per-drug multiplier. MELD is a mortality
  score, not a clearance surrogate, and is not computed.
- Absolute body weight, ideal body weight (Devine) and adjusted body weight are
  all returned; obesity is flagged because the choice changes CrCl by a large
  margin and the guidelines disagree.
"""

import math

from .errors import (
    E_BAD_UNIT,
    E_BAD_VALUE,
    E_MISSING_FIELD,
    E_SCHEMA,
    EngineError,
    W_CG_EGFR_DISCORDANCE,
    W_OBESITY_WEIGHT_AMBIGUITY,
)

SCHEMA_VERSION = "1.0"

REQUIRED = ("age_years", "sex", "weight_kg", "serum_creatinine_mg_dL")


def _num(raw, field, required=True, lo=None, hi=None):
    if field not in raw or raw[field] is None:
        if required:
            raise EngineError(E_MISSING_FIELD, f"{field} is required", field)
        return None
    v = raw[field]
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise EngineError(E_BAD_VALUE, f"{field} must be a number", field)
    v = float(v)
    if math.isnan(v) or math.isinf(v):
        raise EngineError(E_BAD_VALUE, f"{field} is not finite", field)
    if lo is not None and v < lo:
        raise EngineError(E_BAD_VALUE, f"{field} must be >= {lo}", field)
    if hi is not None and v > hi:
        raise EngineError(E_BAD_VALUE, f"{field} must be <= {hi}", field)
    return v


def validate_patient(raw):
    """Normalise and check a patient document. Raises EngineError on refusal."""
    if not isinstance(raw, dict):
        raise EngineError(E_SCHEMA, "patient must be a JSON object")
    version = raw.get("schema_version", SCHEMA_VERSION)
    if version != SCHEMA_VERSION:
        raise EngineError(
            E_SCHEMA, f"unsupported patient schema_version {version!r}", "schema_version"
        )

    p = {
        "schema_version": SCHEMA_VERSION,
        "patient_id": raw.get("patient_id"),
        "age_years": _num(raw, "age_years", lo=0, hi=120),
        "sex": raw.get("sex"),
        "weight_kg": _num(raw, "weight_kg", lo=2, hi=400),
        "height_cm": _num(raw, "height_cm", required=False, lo=30, hi=250),
    }
    if p["sex"] not in ("male", "female"):
        raise EngineError(E_BAD_VALUE, "sex must be 'male' or 'female'", "sex")
    if p["age_years"] < 18:
        raise EngineError(
            "E_AGE_OUT_OF_SCOPE",
            "precision-dose v1 is validated for adults (>= 18 years); paediatric "
            "allometry and maturation are not implemented",
            "age_years",
        )

    scr = raw.get("serum_creatinine_mg_dL")
    if scr is None and raw.get("serum_creatinine_umol_L") is not None:
        scr = _num(raw, "serum_creatinine_umol_L", lo=10, hi=2000) / 88.4
    p["serum_creatinine_mg_dL"] = _num(
        {"serum_creatinine_mg_dL": scr}, "serum_creatinine_mg_dL", lo=0.1, hi=20.0
    )
    p["creatinine_assay"] = raw.get("creatinine_assay", "unspecified")
    if p["creatinine_assay"] not in ("enzymatic", "jaffe", "unspecified"):
        raise EngineError(E_BAD_UNIT, "unknown creatinine_assay", "creatinine_assay")
    p["creatinine_stable"] = bool(raw.get("creatinine_stable", True))

    p["albumin_g_dL"] = _num(raw, "albumin_g_dL", required=False, lo=0.5, hi=6.5)
    p["hematocrit_percent"] = _num(
        raw, "hematocrit_percent", required=False, lo=5, hi=70
    )
    p["days_post_transplant"] = _num(
        raw, "days_post_transplant", required=False, lo=0.0, hi=20000.0
    )
    p["mic_mg_L"] = _num(raw, "mic_mg_L", required=False, lo=0.0, hi=256.0)
    p["indication"] = raw.get("indication")
    p["pregnancy"] = bool(raw.get("pregnancy", False))
    p["smoking_status"] = raw.get("smoking_status", "unknown")

    # genotypes / phenotypes
    genotypes = raw.get("genotypes") or {}
    phenotypes = raw.get("phenotypes") or {}
    if not isinstance(genotypes, dict) or not isinstance(phenotypes, dict):
        raise EngineError(E_SCHEMA, "genotypes/phenotypes must be objects")
    for gene, dip in genotypes.items():
        if not isinstance(dip, str) or "/" not in dip:
            raise EngineError(
                E_BAD_VALUE, f"genotype for {gene} must look like '*1/*3'", "genotypes"
            )
    p["genotypes"] = {str(k): str(v) for k, v in genotypes.items()}
    p["phenotypes"] = {str(k): str(v) for k, v in phenotypes.items()}

    # dialysis
    dialysis = raw.get("renal", {}).get("dialysis") if isinstance(raw.get("renal"), dict) else None
    if dialysis is not None:
        if not isinstance(dialysis, dict):
            raise EngineError(E_SCHEMA, "renal.dialysis must be null or an object")
        dialysis = {
            "modality": dialysis.get("modality", "hd"),
            "sessions_per_week": int(dialysis.get("sessions_per_week", 3)),
            "hours_per_session": float(dialysis.get("hours_per_session", 4.0)),
            "first_session_h": float(dialysis.get("first_session_h", 24.0)),
            "day_interval_h": float(
                dialysis.get("day_interval_h", 48.0)  # Mon/Wed/Fri pattern
            ),
            "blood_flow_mL_min": float(dialysis.get("blood_flow_mL_min", 300.0)),
            "dialysate_flow_mL_min": float(
                dialysis.get("dialysate_flow_mL_min", 500.0)
            ),
            "dialyzer_koA_mL_min": dialysis.get("dialyzer_koA_mL_min"),
        }
        if dialysis["modality"] not in ("hd", "pd", "crrt"):
            raise EngineError(E_BAD_VALUE, "unknown dialysis modality", "renal.dialysis")
    p["dialysis"] = dialysis
    p["renal_status"] = (raw.get("renal") or {}).get("status", "unspecified")

    # hepatic components (all optional; class is only computed when complete)
    hep = raw.get("hepatic") or {}
    comp = hep.get("child_pugh_components") if isinstance(hep, dict) else None
    if comp is not None and not isinstance(comp, dict):
        raise EngineError(E_SCHEMA, "hepatic.child_pugh_components must be an object")
    p["child_pugh_components"] = comp

    cos = raw.get("co_medications") or []
    if not isinstance(cos, list):
        raise EngineError(E_SCHEMA, "co_medications must be a list")
    cleaned = []
    for i, med in enumerate(cos):
        if not isinstance(med, dict) or "name" not in med:
            raise EngineError(E_BAD_VALUE, f"co_medications[{i}] needs a name")
        cleaned.append(
            {
                "name": str(med["name"]).strip().lower(),
                "dose_mg": _num(med, "dose_mg", required=False, lo=0.0),
                "unbound_conc_umol_L": _num(
                    med, "unbound_conc_umol_L", required=False, lo=0.0
                ),
                "steady_state": bool(med.get("steady_state", True)),
                "effects": med.get("effects") or {},
            }
        )
    p["co_medications"] = cleaned

    # Comorbidities are a small controlled vocabulary, not free text: an entry may
    # change a parameter on one of them (digoxin's non-renal clearance is 20 mL/min
    # in heart failure and 41 mL/min without it).
    comorb = raw.get("comorbidities") or []
    if not isinstance(comorb, list):
        raise EngineError(E_SCHEMA, "comorbidities must be a list of strings")
    p["comorbidities"] = [str(c).strip().lower() for c in comorb]
    return p


# --------------------------------------------------------------------------
# body size
# --------------------------------------------------------------------------


def ibw_kg(sex, height_cm):
    """Devine. Returns None when height is unavailable (never imputed)."""
    if height_cm is None:
        return None
    inches = height_cm / 2.54
    base = 50.0 if sex == "male" else 45.5
    return round(base + 2.3 * (inches - 60.0), 4)


def abw_kg(tbw, ibw):
    if ibw is None:
        return None
    return round(ibw + 0.4 * (max(tbw, 0.0) - ibw), 4) if tbw > ibw else round(tbw, 4)


def bsa_m2(height_cm, weight_kg):
    if height_cm is None:
        return None
    return round(math.sqrt(height_cm * weight_kg / 3600.0), 4)


def bmi(weight_kg, height_cm):
    if height_cm is None:
        return None
    return round(weight_kg / (height_cm / 100.0) ** 2, 2)


# --------------------------------------------------------------------------
# renal function
# --------------------------------------------------------------------------


def egfr_ckd_epi_2021(scr_mg_dL, age_years, sex):
    """Inker LA et al., NEJM 2021;385:1737-49. Race-free. mL/min/1.73 m2."""
    if sex == "female":
        k, alpha, female = 0.7, -0.241, 1.012
    else:
        k, alpha, female = 0.9, -0.302, 1.0
    ratio = scr_mg_dL / k
    return (
        142.0
        * min(ratio, 1.0) ** alpha
        * max(ratio, 1.0) ** -1.200
        * 0.9938 ** age_years
        * female
    )


def crcl_cockcroft_gault(scr_mg_dL, age_years, weight_kg, sex):
    """Cockcroft DW, Gault MH. Nephron 1976;16:31-41. mL/min."""
    raw = (140.0 - age_years) * weight_kg / (72.0 * scr_mg_dL)
    return raw * (0.85 if sex == "female" else 1.0)


# --------------------------------------------------------------------------
# hepatic status
# --------------------------------------------------------------------------


def child_pugh(components, albumin_fallback=None):
    """Return {score, class, points, complete}. Missing components -> incomplete."""
    if not components:
        return {"score": None, "class": None, "points": {}, "complete": False}
    points = {}

    bili = components.get("bilirubin_mg_dL")
    if bili is not None:
        points["bilirubin"] = 1 if bili < 2.0 else 2 if bili <= 3.0 else 3
    inr = components.get("inr")
    if inr is not None:
        points["inr"] = 1 if inr < 1.7 else 2 if inr <= 2.3 else 3
    alb = components.get("albumin_g_dL", albumin_fallback)
    if alb is not None:
        points["albumin"] = 1 if alb > 3.5 else 2 if alb >= 2.8 else 3
    asc = components.get("ascites")
    if asc is not None:
        points["ascites"] = {"none": 1, "mild": 2, "moderate_severe": 3}.get(asc)
    ence = components.get("encephalopathy")
    if ence is not None:
        points["encephalopathy"] = {"none": 1, "grade_1_2": 2, "grade_3_4": 3}.get(ence)

    if len(points) < 5 or any(v is None for v in points.values()):
        return {
            "score": sum(v for v in points.values() if v is not None) if points else None,
            "class": None,
            "points": points,
            "complete": False,
        }
    score = sum(points.values())
    cls = "A" if score <= 6 else "B" if score <= 9 else "C"
    return {"score": score, "class": cls, "points": points, "complete": True}


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------


def derive(patient):
    """All derived quantities, each labelled with the formula that produced it."""
    sex = patient["sex"]
    age = patient["age_years"]
    tbw = patient["weight_kg"]
    height = patient.get("height_cm")
    scr = patient["serum_creatinine_mg_dL"]

    ibw = ibw_kg(sex, height)
    abw = abw_kg(tbw, ibw)
    derived = {
        "bmi_kg_m2": {"value": bmi(tbw, height), "formula": "F_BMI"},
        "bsa_m2": {"value": bsa_m2(height, tbw), "formula": "F_BSA_MOSTELLER"},
        "ibw_kg": {"value": ibw, "formula": "F_IBW_DEVINE"},
        "abw_kg": {"value": abw, "formula": "F_ABW_04"},
        "egfr_ckd_epi_2021_mL_min_1_73m2": {
            "value": round(egfr_ckd_epi_2021(scr, age, sex), 4),
            "formula": "F_EGFR_CKDEPI2021",
            "source": "Inker LA et al. NEJM 2021;385:1737-49",
        },
        "crcl_cg_mL_min": {
            "value": {
                "total_body_weight": round(
                    crcl_cockcroft_gault(scr, age, tbw, sex), 4
                ),
                "ideal_body_weight": (
                    round(crcl_cockcroft_gault(scr, age, ibw, sex), 4)
                    if ibw is not None
                    else None
                ),
                "adjusted_body_weight": (
                    round(crcl_cockcroft_gault(scr, age, abw, sex), 4)
                    if abw is not None
                    else None
                ),
            },
            "formula": "F_CRCL_CG",
            "source": "Cockcroft DW, Gault MH. Nephron 1976;16:31-41",
        },
        "child_pugh": child_pugh(
            patient.get("child_pugh_components"), patient.get("albumin_g_dL")
        ),
    }

    warnings = []
    if not patient["creatinine_stable"]:
        warnings.append(
            {
                "code": "W_SCr_NOT_STEADY",
                "message": "serum creatinine flagged as not at steady state (AKI or "
                "rapidly changing renal function): CrCl and eGFR are unreliable; "
                "prefer TDM-guided dosing.",
            }
        )

    crcl_tbw = derived["crcl_cg_mL_min"]["value"]["total_body_weight"]
    egfr = derived["egfr_ckd_epi_2021_mL_min_1_73m2"]["value"]
    if crcl_tbw and egfr:
        disc = abs(crcl_tbw - egfr) / max(egfr, 1e-9)
        derived["crcl_vs_egfr_relative_difference"] = {
            "value": round(disc, 4),
            "formula": "F_DISCORDANCE",
        }
        if disc > 0.30:
            warnings.append(
                {
                    "code": "W_CG_EGFR_DISCORDANCE",
                    "message": f"Cockcroft-Gault CrCl ({crcl_tbw:.0f} mL/min) and "
                    f"CKD-EPI 2021 eGFR ({egfr:.0f}) differ by {disc*100:.0f}%. "
                    "This is common in obesity, low muscle mass and extremes of age; "
                    "the drug entry declares which metric its clearance model used.",
                }
            )
    bmi_v = derived["bmi_kg_m2"]["value"]
    if bmi_v is not None and bmi_v >= 30.0:
        warnings.append(
            {
                "code": "W_OBESITY_WEIGHT_AMBIGUITY",
                "message": f"BMI {bmi_v:.1f}: CrCl depends on the weight metric "
                f"(TBW {crcl_tbw:.0f}, IBW "
                f"{derived['crcl_cg_mL_min']['value']['ideal_body_weight']}, ABW "
                f"{derived['crcl_cg_mL_min']['value']['adjusted_body_weight']} "
                "mL/min). TBW overestimates CrCl in obesity; ABW is the usual choice "
                "for maintenance dosing, TBW for loading (volume).",
            }
        )
    if height is None:
        warnings.append(
            {
                "code": "W_HEPATIC_DATA_MISSING",
                "message": "height is missing: IBW/ABW are unavailable, so "
                "weight-metric-dependent calculations cannot be performed.",
            }
        )
    cp = derived["child_pugh"]
    if not cp["complete"]:
        warnings.append(
            {
                "code": "W_HEPATIC_DATA_MISSING",
                "message": "Child-Pugh is incomplete: hepatic clearance scaling is not "
                "applied. Supply bilirubin, INR, albumin, ascites and encephalopathy "
                "to enable it.",
            }
        )
    return derived, warnings


def renal_value(patient, derived, metric):
    """Resolve a declared renal metric id to a number, or raise."""
    v = derived["crcl_cg_mL_min"]["value"]
    mapping = {
        "crcl_cg_total": v["total_body_weight"],
        "crcl_cg_ibw": v["ideal_body_weight"],
        "crcl_cg_abw": v["adjusted_body_weight"],
        "egfr_ckd_epi_2021": derived["egfr_ckd_epi_2021_mL_min_1_73m2"]["value"],
    }
    if metric not in mapping:
        raise EngineError(E_BAD_VALUE, f"unknown renal metric {metric!r}", "renal")
    if mapping[metric] is None:
        raise EngineError(
            E_MISSING_FIELD,
            f"renal metric {metric} needs a height (for IBW/ABW); provide height_cm",
            "height_cm",
        )
    return float(mapping[metric])
