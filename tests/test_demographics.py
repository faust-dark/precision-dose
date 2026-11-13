"""Body size, renal function and hepatic classification.

The reference values in these tests are hand-computed from the published
equations (the arithmetic is written out in the comments) rather than copied
from another implementation, so a systematic error in this module cannot hide
behind a matching bug elsewhere.
"""

import unittest

from precision_dose import clin
from precision_dose.errors import EngineError


def patient(**over):
    base = {
        "schema_version": "1.0",
        "age_years": 60,
        "sex": "male",
        "weight_kg": 70.0,
        "height_cm": 175.0,
        "serum_creatinine_mg_dL": 1.0,
    }
    base.update(over)
    return clin.validate_patient(base)


class TestCKDEPI2021(unittest.TestCase):
    def test_male_reference(self):
        # ratio = 1.0/0.9 = 1.1111; max(ratio,1)^-1.2 = 0.88125
        # 142 * 0.88125 = 125.14; 0.9938^60 = 0.68876; product = 86.19
        got = clin.egfr_ckd_epi_2021(1.0, 60, "male")
        self.assertAlmostEqual(got, 86.19, delta=0.05)

    def test_female_reference(self):
        # ratio = 1.42857; ^-1.2 = 0.65183; 142*0.65183 = 92.56;
        # * 0.68876 * 1.012 = 64.53
        got = clin.egfr_ckd_epi_2021(1.0, 60, "female")
        self.assertAlmostEqual(got, 64.53, delta=0.05)

    def test_no_race_term_and_monotone_in_scr(self):
        vals = [clin.egfr_ckd_epi_2021(s, 55, "male") for s in (0.7, 1.0, 1.5, 2.5)]
        self.assertEqual(vals, sorted(vals, reverse=True))

    def test_young_beats_old(self):
        self.assertGreater(
            clin.egfr_ckd_epi_2021(1.0, 30, "male"), clin.egfr_ckd_epi_2021(1.0, 80, "male")
        )


class TestCockcroftGault(unittest.TestCase):
    def test_male(self):
        # (140-60)*70/(72*1.0) = 77.78
        self.assertAlmostEqual(clin.crcl_cockcroft_gault(1.0, 60, 70, "male"), 77.78, 2)

    def test_female_factor(self):
        self.assertAlmostEqual(
            clin.crcl_cockcroft_gault(1.0, 60, 70, "female"), 77.78 * 0.85, 2
        )


class TestBodySize(unittest.TestCase):
    def test_devine(self):
        # 175 cm = 68.9 in -> 50 + 2.3*8.9 = 70.47 kg
        self.assertAlmostEqual(clin.ibw_kg("male", 175.0), 70.47, 1)
        # 163 cm = 64.17 in -> 45.5 + 2.3*4.17 = 55.09 kg
        self.assertAlmostEqual(clin.ibw_kg("female", 163.0), 55.09, 1)

    def test_abw(self):
        self.assertAlmostEqual(clin.abw_kg(90.0, 70.0), 78.0, 2)
        # below IBW the actual weight is used, not an inflated value
        self.assertAlmostEqual(clin.abw_kg(60.0, 70.0), 60.0, 2)

    def test_missing_height_is_none_not_guess(self):
        self.assertIsNone(clin.ibw_kg("male", None))
        self.assertIsNone(clin.abw_kg(80.0, None))


class TestChildPugh(unittest.TestCase):
    def test_class_a(self):
        cp = clin.child_pugh(
            {
                "bilirubin_mg_dL": 1.0,
                "inr": 1.1,
                "albumin_g_dL": 4.0,
                "ascites": "none",
                "encephalopathy": "none",
            }
        )
        self.assertTrue(cp["complete"])
        self.assertEqual((cp["score"], cp["class"]), (5, "A"))

    def test_class_c(self):
        cp = clin.child_pugh(
            {
                "bilirubin_mg_dL": 4.0,
                "inr": 2.8,
                "albumin_g_dL": 2.2,
                "ascites": "moderate_severe",
                "encephalopathy": "grade_3_4",
            }
        )
        self.assertEqual((cp["score"], cp["class"]), (15, "C"))

    def test_incomplete_never_assigns_a_class(self):
        cp = clin.child_pugh({"bilirubin_mg_dL": 1.0, "inr": 1.1})
        self.assertFalse(cp["complete"])
        self.assertIsNone(cp["class"])


class TestPatientValidation(unittest.TestCase):
    def test_paediatric_refused(self):
        with self.assertRaises(EngineError) as ctx:
            patient(age_years=15)
        self.assertEqual(ctx.exception.code, "E_AGE_OUT_OF_SCOPE")

    def test_missing_creatinine_refused(self):
        doc = {
            "schema_version": "1.0",
            "age_years": 50,
            "sex": "male",
            "weight_kg": 70,
            "height_cm": 175,
        }
        with self.assertRaises(EngineError) as ctx:
            clin.validate_patient(doc)
        self.assertEqual(ctx.exception.code, "E_MISSING_FIELD")

    def test_micromolar_creatinine_converted(self):
        p = clin.validate_patient(
            {
                "schema_version": "1.0",
                "age_years": 50,
                "sex": "male",
                "weight_kg": 70,
                "height_cm": 175,
                "serum_creatinine_umol_L": 88.4,
            }
        )
        self.assertAlmostEqual(p["serum_creatinine_mg_dL"], 1.0, 6)

    def test_obesity_warning_and_metric_spread(self):
        p = patient(weight_kg=120.0, height_cm=170.0)
        derived, warnings = clin.derive(p)
        codes = {w["code"] for w in warnings}
        self.assertIn("W_OBESITY_WEIGHT_AMBIGUITY", codes)
        v = derived["crcl_cg_mL_min"]["value"]
        self.assertGreater(v["total_body_weight"], v["ideal_body_weight"])

    def test_renal_metric_requiring_height_raises(self):
        p = patient(height_cm=None, weight_kg=80.0)
        derived, _ = clin.derive(p)
        with self.assertRaises(EngineError):
            clin.renal_value(p, derived, "crcl_cg_ibw")


if __name__ == "__main__":
    unittest.main()
