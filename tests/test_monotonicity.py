"""Monotonicity and physical-invariance properties.

These are the properties a mechanistic model must not violate. They are checked
through the public engine, not the internals, because the engine is what a
pharmacist runs.
"""

import unittest

from precision_dose import engine, library

LIB, INTER, _ = library.load_library()


def patient(**over):
    base = {
        "schema_version": "1.0",
        "age_years": 60,
        "sex": "male",
        "weight_kg": 75.0,
        "height_cm": 175.0,
        "serum_creatinine_mg_dL": 1.2,
        "albumin_g_dL": 3.8,
        "days_post_transplant": 100,
        "genotypes": {"CYP3A5": "*3/*3"},
    }
    base.update(over)
    return base


class TestMonotonicity(unittest.TestCase):
    def test_auc_increases_with_dose(self):
        entry = LIB["vancomycin"]
        from precision_dose import clin

        p = clin.validate_patient(patient())
        derived, warns = clin.derive(p)
        params, _ = engine.individualize(entry, p, derived, INTER, warns)
        model = engine.build_model(entry, params)
        target = entry["targets"][0]
        prev = 0.0
        for dose in (500.0, 1000.0, 2000.0):
            ev = engine.evaluate_regimen(
                entry, model, __import__("precision_dose.pk", fromlist=["x"]).Regimen(
                    [(0.0, dose)], 12.0, 2.0, "iv"
                ),
                target, 1.0,
            )
            self.assertGreater(ev["auc24"], prev)
            prev = ev["auc24"]

    def test_clearance_increases_with_renal_function(self):
        entry = LIB["vancomycin"]
        from precision_dose import clin

        prev = 0.0
        for scr in (3.0, 2.0, 1.2, 0.8):
            p = clin.validate_patient(patient(serum_creatinine_mg_dL=scr))
            derived, warns = clin.derive(p)
            params, _ = engine.individualize(entry, p, derived, INTER, warns)
            self.assertGreater(params["cl"], prev)
            prev = params["cl"]

    def test_stronger_inhibition_lowers_clearance(self):
        entry = LIB["tacrolimus"]
        from precision_dose import clin

        cl = []
        for meds in ([], [{"name": "fluconazole", "dose_mg": 400}],
                     [{"name": "ketoconazole", "dose_mg": 400}]):
            p = clin.validate_patient(patient(co_medications=meds))
            derived, warns = clin.derive(p)
            params, _ = engine.individualize(entry, p, derived, INTER, warns)
            cl.append(params["cl"])
        self.assertGreater(cl[0], cl[1])
        self.assertGreater(cl[1], cl[2])

    def test_volume_increases_with_weight(self):
        entry = LIB["lithium"]
        from precision_dose import clin

        v = []
        for w in (55.0, 70.0, 95.0):
            p = clin.validate_patient(patient(weight_kg=w, serum_creatinine_mg_dL=0.9))
            derived, warns = clin.derive(p)
            params, _ = engine.individualize(entry, p, derived, INTER, warns)
            v.append(params["v1"])
        self.assertEqual(v, sorted(v))


class TestAttributionExactness(unittest.TestCase):
    def test_log_additivity_residual_is_zero(self):
        for drug in ("vancomycin", "tacrolimus", "phenytoin", "lithium"):
            res = engine.run_dose(
                patient(co_medications=[{"name": "fluconazole", "dose_mg": 400}]),
                drug, pta_n=200,
            )
            for key, block in (res.get("attribution") or {}).items():
                self.assertLess(abs(block["log_additivity_residual"]), 1e-9, (drug, key))

    def test_attribution_reproduces_clearance(self):
        res = engine.run_dose(patient(), "vancomycin", pta_n=200)
        block = res["attribution"]["cl"]
        self.assertAlmostEqual(
            block["reference_value"] * block["total_multiplier"],
            block["individual_value"],
            delta=1e-6 * block["individual_value"],
        )


if __name__ == "__main__":
    unittest.main()
