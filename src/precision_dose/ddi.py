"""Drug-drug interaction algebra.

Two mechanisms are supported, and they are not interchangeable:

- competitive / noncompetitive: instantaneous reversible inhibition described by
  Ki and the unbound systemic concentration [I]u:
      CL_int,app / CL_int = 1 / (1 + [I]u / Ki)
- mechanism-based inactivation: time-dependent loss of enzyme activity:
      activity = kdeg / (kdeg + kinact * [I]u / (KI + [I]u))
- induction: Emax model on intrinsic clearance, with a time to steady state.

All effects are applied through the fraction of clearance that the pathway
actually accounts for (fm). An inhibitor of an enzyme the drug does not use
changes nothing, and the attribution output says so.
"""

from .errors import E_BAD_VALUE, EngineError


def unbound_concentration(med, perpetrator_pk):
    """Estimate [I]u (umol/L) at steady state for a co-medication.

    Preference order: caller-supplied [I]u, then a one-compartment Cmax,ss
    estimate from the curated perpetrator PK. If neither is available we return
    None and the caller must warn -- never assume zero and never assume a value.
    """
    if med.get("unbound_conc_umol_L") is not None:
        return float(med["unbound_conc_umol_L"]), "caller-supplied unbound concentration"
    if not perpetrator_pk:
        return None, "no perpetrator PK in the interaction library"
    dose = med.get("dose_mg")
    if dose is None or dose <= 0:
        return None, "no dose supplied for the co-medication"
    v_l = perpetrator_pk["v_l_per_kg"] * perpetrator_pk.get("weight_kg", 70.0)
    cl = perpetrator_pk["cl_l_h"]
    f = perpetrator_pk.get("f", 1.0)
    fu = perpetrator_pk.get("fu", 1.0)
    mw = perpetrator_pk.get("mw_g_per_mol")
    if not mw:
        return None, "perpetrator MW missing, cannot convert mg/L to umol/L"
    import math

    k = cl / v_l
    tau = 24.0
    cmax_ss = (f * dose / v_l) / (1.0 - math.exp(-k * tau))  # mg/L, 1-cmt approx
    iu = cmax_ss * fu / mw * 1000.0  # mg/L -> umol/L
    return iu, "1-cmt Cmax,ss estimate from curated perpetrator PK"


def pathway_multiplier(effect, iu):
    """Multiplier on intrinsic clearance for one pathway."""
    mech = effect.get("mechanism")
    if mech in ("competitive", "noncompetitive"):
        ki = effect.get("ki_umol_L")
        if ki is None or ki <= 0:
            raise EngineError(
                E_BAD_VALUE, "competitive inhibition needs a positive ki_umol_L", "effects"
            )
        return 1.0 / (1.0 + iu / ki)
    if mech == "mechanism_based":
        ki = effect.get("ki_umol_L")
        kinact = effect.get("kinact_per_h")
        kdeg = effect.get("kdeg_per_h")
        if not ki or not kinact or not kdeg:
            raise EngineError(
                E_BAD_VALUE,
                "mechanism-based inhibition needs ki_umol_L, kinact_per_h and kdeg_per_h",
                "effects",
            )
        occupancy = iu / (ki + iu)
        return kdeg / (kdeg + kinact * occupancy)
    if mech == "induction":
        emax = effect.get("emax")
        ec50 = effect.get("ec50_umol_L")
        if emax is None or not ec50:
            raise EngineError(
                E_BAD_VALUE, "induction needs emax and ec50_umol_L", "effects"
            )
        return 1.0 + emax * iu / (ec50 + iu)
    if mech == "empirical_fold":
        fold = effect.get("cl_fold")
        if fold is None:
            raise EngineError(E_BAD_VALUE, "empirical_fold needs cl_fold", "effects")
        return float(fold)
    raise EngineError(E_BAD_VALUE, f"unknown DDI mechanism {mech!r}", "effects")


def evaluate(patient, fm_by_pathway, interactions, warnings):
    """Return (cl_multiplier, contributions, notes).

    contributions: one entry per (perpetrator, pathway) with the log-effect,
    which is what the attribution engine later reports verbatim.
    """
    contributions = []
    multiplier = 1.0
    for med in patient.get("co_medications", []):
        name = med["name"]
        perp = (interactions.get("perpetrators") or {}).get(name)
        declared = med.get("effects") or {}
        effects = dict(perp.get("effects", {})) if perp else {}
        if declared:
            # caller-supplied effects override library effects for this patient
            effects.update({k: v for k, v in declared.items()})
        if not effects:
            warnings.append(
                {
                    "code": "W_DDI_PRESENT",
                    "message": f"{name} is on the medication list but has no curated "
                    "interaction entry; no DDI effect was applied.",
                }
            )
            continue
        pk_entry = (interactions.get("perpetrator_pk") or {}).get(name)
        iu, basis = unbound_concentration(med, pk_entry)
        if iu is None:
            warnings.append(
                {
                    "code": "W_DDI_PRESENT",
                    "message": f"{name}: unbound concentration unavailable ({basis}); "
                    "the interaction was NOT applied.",
                }
            )
            continue
        for pathway, effect in effects.items():
            fm = fm_by_pathway.get(pathway)
            if fm is None or fm <= 0:
                contributions.append(
                    {
                        "perpetrator": name,
                        "pathway": pathway,
                        "applied": False,
                        "reason": "the drug's clearance does not depend on this pathway",
                        "iu_umol_L": round(iu, 4),
                    }
                )
                continue
            path_mult = pathway_multiplier(effect, iu)
            combined = 1.0 + fm * (path_mult - 1.0)
            multiplier *= combined
            contributions.append(
                {
                    "perpetrator": name,
                    "pathway": pathway,
                    "mechanism": effect.get("mechanism"),
                    "applied": True,
                    "iu_umol_L": round(iu, 4),
                    "ki_umol_L": effect.get("ki_umol_L"),
                    "pathway_multiplier": round(path_mult, 6),
                    "fm": fm,
                    "cl_multiplier": round(combined, 6),
                    "delta_cl_percent": round(100.0 * (combined - 1.0), 4),
                    "iu_basis": basis,
                    "source": (perp or {}).get("source") or effect.get("source"),
                }
            )
        if effects:
            warnings.append(
                {
                    "code": "W_DDI_PRESENT",
                    "message": f"{name}: interaction applied at [I]u {iu:.2f} umol/L "
                    f"({basis}).",
                }
            )
    return multiplier, contributions, warnings
