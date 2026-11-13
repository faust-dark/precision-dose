"""Bayesian estimation: recovery, refusal and honesty about weak information."""

import math
import unittest

from precision_dose import bayes, engine, library, pk

LIB, INTER, _ = library.load_library()


def patient():
    return {
        "schema_version": "1.0",
        "age_years": 60,
        "sex": "male",
        "weight_kg": 75.0,
        "height_cm": 175.0,
        "serum_creatinine_mg_dL": 1.2,
        "albumin_g_dL": 3.8,
        "indication": "mrsa_bacteremia",
        "mic_mg_L": 1.0,
    }


def prior_params(entry, p_doc):
    from precision_dose import clin

    p = clin.validate_patient(p_doc)
    derived, warns = clin.derive(p)
    params, _ = engine.individualize(entry, p, derived, INTER, warns)
    return params


class TestMAPRecovery(unittest.TestCase):
    def test_recovers_clearance_within_tolerance(self):
        entry = LIB["vancomycin"]
        truth = {"cl": 2.4, "v1": 40.0, "v2": 45.0, "q": 7.0}
        model = pk.LinearPK(truth["cl"], truth["v1"], truth["v2"], truth["q"])
        dosing = {
            "doses": [
                {"dose": 1500, "start_h": 0.0, "infusion_min": 120, "route": "iv"},
                {"dose": 1000, "start_h": 12.0, "infusion_min": 120, "route": "iv"},
                {"dose": 1000, "start_h": 24.0, "infusion_min": 120, "route": "iv"},
            ],
            "levels": [],
        }
        times = [14.0, 24.0]
        preds = model.simulate(
            [(d["start_h"], d["dose"], 2.0, "iv") for d in dosing["doses"]], times
        )
        dosing["levels"] = [
            {"time_h": t, "value": round(c, 3)} for t, c in preds
        ]
        prior = prior_params(entry, patient())
        fit = bayes.map_estimate(entry, prior, dosing)
        self.assertTrue(fit["converged"])
        recovered = fit["posterior"]["cl"]
        # empirical Bayes shrinks toward the prior, so the correct assertion is
        # "the level moved the estimate toward the truth", not "it hit the truth"
        self.assertLess(
            abs(recovered - truth["cl"]), abs(prior["cl"] - truth["cl"]), fit
        )
        self.assertLess(abs(recovered - truth["cl"]) / truth["cl"], 0.25, fit)

    def test_single_level_reports_weak_information(self):
        entry = LIB["vancomycin"]
        prior = prior_params(entry, patient())
        dosing = {
            "doses": [{"dose": 1500, "start_h": 0.0, "infusion_min": 120, "route": "iv"}],
            "levels": [{"time_h": 14.0, "value": 18.0}],
        }
        fit = bayes.map_estimate(entry, prior, dosing)
        ir = fit["information_ratio"]["cl"]
        self.assertIsNotNone(ir)
        self.assertLess(ir, 0.9)

    def test_absurd_levels_still_produce_a_finite_fit_or_refuse(self):
        entry = LIB["vancomycin"]
        prior = prior_params(entry, patient())
        dosing = {
            "doses": [{"dose": 1000, "start_h": 0.0, "infusion_min": 60, "route": "iv"}],
            "levels": [{"time_h": 8.0, "value": 400.0}],
        }
        try:
            fit = bayes.map_estimate(entry, prior, dosing)
        except Exception as exc:  # a refusal is acceptable, a crash is not
            self.assertIn("not converge", str(exc).lower())
            return
        self.assertTrue(math.isfinite(fit["posterior"]["cl"]))

    def test_levels_without_doses_refused(self):
        entry = LIB["vancomycin"]
        prior = prior_params(entry, patient())
        with self.assertRaises(Exception):
            bayes.map_estimate(entry, prior, {"doses": [], "levels": [{"time_h": 2, "value": 10}]})


if __name__ == "__main__":
    unittest.main()
