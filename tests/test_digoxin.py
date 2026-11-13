"""Digoxin-specific behaviour: the heart-failure switch, the P-gp interaction, and
the sampling window a level has to respect."""

import unittest

from precision_dose import bayes, engine, library

LIB, INTER, _ = library.load_library()

BASE = {
    "schema_version": "1.0",
    "age_years": 78,
    "sex": "female",
    "weight_kg": 60.0,
    "height_cm": 158.0,
    "serum_creatinine_mg_dL": 1.4,
    "albumin_g_dL": 3.6,
    "indication": "hfref_with_atrial_fibrillation",
}


def patient(**over):
    doc = dict(BASE)
    doc.update(over)
    return doc


class TestDigoxin(unittest.TestCase):
    def test_heart_failure_lowers_clearance(self):
        without = engine.run_dose(patient(), "digoxin", pta_n=300)
        with_hf = engine.run_dose(patient(comorbidities=["heart_failure"]), "digoxin", pta_n=300)
        cl_no = without["parameters"]["values"]["cl"]
        cl_hf = with_hf["parameters"]["values"]["cl"]
        # 41 vs 20 mL/min of non-renal clearance
        self.assertAlmostEqual(cl_no - cl_hf, (41 - 20) * 60 / 1000, 2)
        self.assertGreater(with_hf["regimen"]["maintenance"]["dose"],
                           0)  # still doseable, but a step lower or equal
        self.assertLessEqual(
            with_hf["regimen"]["maintenance"]["dose"], without["regimen"]["maintenance"]["dose"]
        )

    def test_amiodarone_halves_clearance(self):
        plain = engine.run_dose(patient(comorbidities=["heart_failure"]), "digoxin", pta_n=300)
        combo = engine.run_dose(
            patient(
                comorbidities=["heart_failure"],
                co_medications=[{"name": "amiodarone", "dose_mg": 200}],
            ),
            "digoxin",
            pta_n=300,
        )
        self.assertAlmostEqual(
            combo["parameters"]["values"]["cl"] / plain["parameters"]["values"]["cl"], 0.5, 2
        )
        # the additive component is scaled by the interaction too, and says so
        f = combo["attribution"]["cl"]["factors"]
        self.assertTrue(any(x["factor"] == "ddi" for x in f))

    def test_amiodarone_does_not_touch_vancomycin(self):
        # the P-gp fold must be gated on the substrate declaring that pathway
        res = engine.run_dose(
            patient(
                serum_creatinine_mg_dL=1.2,
                comorbidities=[],
                co_medications=[{"name": "amiodarone", "dose_mg": 200}],
            ),
            "vancomycin",
            pta_n=200,
        )
        self.assertFalse(
            any(f["factor"] == "ddi" for f in res["attribution"]["cl"]["factors"])
        )

    def test_level_inside_distribution_phase_warns(self):
        entry = LIB["digoxin"]
        params = engine.individualize(
            entry, __import__("precision_dose.clin", fromlist=["x"]).validate_patient(patient()),
            __import__("precision_dose.clin", fromlist=["x"]).derive(
                __import__("precision_dose.clin", fromlist=["x"]).validate_patient(patient())
            )[0],
            INTER, [],
        )[0]
        fit = bayes.map_estimate(
            entry,
            params,
            {
                "doses": [{"dose": 125, "start_h": 0, "infusion_min": 0, "route": "oral"}],
                "levels": [{"time_h": 3.0, "value": 1.4}],
            },
        )
        codes = {w["code"] for w in fit["warnings"]}
        self.assertIn("W_LEVEL_TOO_EARLY", codes)

    def test_default_target_is_the_heart_failure_range(self):
        res = engine.run_dose(patient(comorbidities=["heart_failure"]), "digoxin", pta_n=200)
        self.assertEqual(res["question"]["target_id"], "trough_hf")
        self.assertEqual((res["pta"]["target_low"], res["pta"]["target_high"]), (0.5, 0.9))

    def test_dose_is_a_tablet_multiple(self):
        res = engine.run_dose(patient(comorbidities=["heart_failure"]), "digoxin", pta_n=200)
        dose = res["regimen"]["maintenance"]["dose"]
        self.assertAlmostEqual(dose % 62.5, 0.0, 6)


if __name__ == "__main__":
    unittest.main()
