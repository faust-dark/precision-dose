"""Bayesian TDM: MAP estimation of individual PK parameters.

The objective is the standard empirical Bayes one:

    OFV(eta) = sum_j (Cobs_j - Cpred_j)^2 / var_j  +  sum_i eta_i^2 / omega_i^2

with a combined additive/proportional assay-plus-model error and the library's
between-subject variability as the prior. The optimiser is deterministic (no
random restarts) and the caller is forced to look at convergence: an
unconverged fit raises rather than returning a plausible-looking number.

Two things this deliberately reports that a point estimate hides:

- the prior-to-posterior information ratio, which answers "is this answer
  actually informed by the levels, or is it the prior wearing a lab coat?",
- the posterior standard errors, which are propagated into the re-dosing PTA.
"""

import math

from . import engine, pk
from .errors import E_NOT_CONVERGED, E_SCHEMA, EngineError, W_LEVEL_NEAR_DOSE, W_WEAK_INFORMATION
from .matmath import nelder_mead, solve


def _parse(dosing_doc, entry):
    if not isinstance(dosing_doc, dict):
        raise EngineError(E_SCHEMA, "the TDM document must be a JSON object")
    doses_raw = dosing_doc.get("doses")
    levels_raw = dosing_doc.get("levels")
    if not doses_raw or not levels_raw:
        raise EngineError(E_SCHEMA, "the TDM document needs 'doses' and 'levels'")
    default_inf = (entry.get("constraints") or {}).get("default_infusion_min", 0)
    doses = []
    for i, d in enumerate(doses_raw):
        amount = d.get("dose", d.get("dose_mg"))
        if amount is None or "start_h" not in d:
            raise EngineError(E_SCHEMA, f"doses[{i}] needs dose and start_h")
        doses.append(
            (
                float(d["start_h"]),
                float(amount),
                float(d.get("infusion_min", default_inf)) / 60.0,
                d.get("route", entry["model"].get("route_default", "iv")),
            )
        )
    levels = []
    for i, l in enumerate(levels_raw):
        if "time_h" not in l or "value" not in l:
            raise EngineError(E_SCHEMA, f"levels[{i}] needs time_h and value")
        levels.append(
            {
                "time_h": float(l["time_h"]),
                "value": float(l["value"]),
                "matrix": l.get("matrix", entry["model"].get("concentration_matrix")),
                "assay": l.get("assay"),
            }
        )
    levels.sort(key=lambda x: x["time_h"])
    return doses, levels


def map_estimate(entry, params, dosing_doc):
    doses, levels = _parse(dosing_doc, entry)
    iiv = entry["model"].get("iiv") or {}
    residual = entry["model"].get("residual") or {}
    cv = float(residual.get("cv_percent", 0.0)) / 100.0
    sd = float(residual.get("sd", 0.0))

    free = ["cl", "v1"] + (["v2"] if "v2" in params else [])
    prior_cv = {
        "cl": float(iiv.get("cl_cv_percent", 0.0)) / 100.0,
        "v1": float(iiv.get("v1_cv_percent", 0.0)) / 100.0,
        "v2": float(iiv.get("v2_cv_percent", 0.0)) / 100.0,
    }
    times = [l["time_h"] for l in levels]
    obs = [l["value"] for l in levels]

    warnings = []
    min_time = ((entry.get("monitoring") or {}).get("level_min_time_h")) or 0.0
    for l in levels:
        for start, _, dur, _ in doses:
            elapsed = l["time_h"] - (start + dur)
            if 0.0 <= elapsed < float(min_time):
                warnings.append(
                    {
                        "code": "W_LEVEL_TOO_EARLY",
                        "message": f"a level {elapsed:.1f} h after the end of an infusion is "
                        f"inside this drug's distribution phase ({min_time:g} h): the "
                        "concentration over-reads the steady state and the fit will "
                        "under-estimate clearance.",
                    }
                )
    for l in levels:
        for start, _, dur, _ in doses:
            if 0.0 <= l["time_h"] - (start + dur) < 0.5:
                warnings.append(
                    {
                        "code": W_LEVEL_NEAR_DOSE,
                        "message": f"a level at {l['time_h']:.2f} h is within 30 min of the "
                        "end of an infusion: distribution is incomplete and the residual "
                        "will be inflated.",
                    }
                )

    def predict(eta):
        p = dict(params)
        for i, name in enumerate(free):
            p[name] = params[name] * math.exp(eta[i])
        model = engine.build_model(entry, p)
        preds = [c for _, c in model.simulate(doses, times)]
        return p, preds

    def ofv(eta):
        _, preds = predict(eta)
        total = 0.0
        for o, pr in zip(obs, preds):
            var = (cv * max(pr, 1e-9)) ** 2 + sd ** 2
            total += (o - pr) ** 2 / var
        for i, name in enumerate(free):
            if prior_cv[name] > 0:
                total += (eta[i] / prior_cv[name]) ** 2
        return total

    x0 = [0.0] * len(free)
    x_best, f_best, info = nelder_mead(ofv, x0, step=0.35, max_iter=4000)
    if not info["converged"]:
        raise EngineError(
            E_NOT_CONVERGED,
            "the MAP fit did not converge; refusing to report parameters. Check level "
            "timing and whether the doses in the history match what was given.",
            details=info,
        )
    post_params, preds = predict(x_best)

    # numerical Hessian -> covariance (OFV = -2 log L, so cov = 2 H^-1)
    k = len(free)
    h = 1e-3
    hess = [[0.0] * k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            def at(di, dj):
                e = list(x_best)
                e[i] += di * h
                e[j] += dj * h
                return ofv(e)

            hess[i][j] = (at(1, 1) - at(1, -1) - at(-1, 1) + at(-1, -1)) / (4 * h * h)
    try:
        ident = [[1.0 if a == b else 0.0 for b in range(k)] for a in range(k)]
        inv_cols = [solve(hess, ident[r]) for r in range(k)]
        cov = [[2.0 * inv_cols[c][r] for c in range(k)] for r in range(k)]
        var = [max(cov[i][i], 0.0) for i in range(k)]
    except ZeroDivisionError:
        cov = None
        var = [float("nan")] * k

    se = [math.sqrt(v) if v == v else float("nan") for v in var]
    info_ratio = {}
    for i, name in enumerate(free):
        prior_var = prior_cv[name] ** 2
        if prior_var > 0 and var[i] == var[i]:
            info_ratio[name] = round(1.0 - var[i] / prior_var, 4)
        else:
            info_ratio[name] = None
    if any(v is not None and v < 0.2 for v in info_ratio.values()):
        warnings.append(
            {
                "code": W_WEAK_INFORMATION,
                "message": "at least one parameter is still dominated by the prior "
                f"(information ratio {info_ratio}): this is Bayesian shrinkage, not "
                "measurement. More or better-timed levels are needed.",
            }
        )

    residuals = []
    for l, o, pr in zip(levels, obs, preds):
        sigma = math.sqrt((cv * max(pr, 1e-9)) ** 2 + sd ** 2)
        residuals.append(
            {
                "time_h": l["time_h"],
                "observed": o,
                "predicted": round(pr, 4),
                "weighted_residual": round((o - pr) / sigma, 4) if sigma else None,
                "percent_error": round(100.0 * (pr - o) / o, 2) if o else None,
            }
        )

    return {
        "converged": True,
        "optimizer": info,
        "ofv": round(f_best, 6),
        "free_parameters": free,
        "eta": {name: round(x_best[i], 6) for i, name in enumerate(free)},
        "prior": {name: round(params[name], 6) for name in free},
        "posterior": {name: round(post_params[name], 6) for name in free},
        "posterior_cv_percent": {
            name: (round(100.0 * se[i] / post_params[name], 4) if se[i] == se[i] else None)
            for i, name in enumerate(free)
        },
        "information_ratio": info_ratio,
        "levels": residuals,
        "n_levels": len(levels),
        "covariance": [[round(v, 10) for v in row] for row in cov] if cov else None,
        "warnings": warnings,
        "note": (
            "information_ratio = 1 - posterior_variance/prior_variance. Below 0.2 the "
            "level has told us little; the number is mostly the population prior. "
            "Posterior CVs are propagated into the re-dosing PTA."
        ),
    }
