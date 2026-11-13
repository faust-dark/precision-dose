"""Edge cases: dialysis, severe hepatic impairment, obesity, age extremes,
missing data, and the refusal paths that must not silently impute."""

import copy
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from precision_dose import engine, library
from precision_dose.errors import EngineError

LIB, INTER, _ = library.load_library()
ROOT = Path(__file__).resolve().parents[1]


def base(**over):
    doc = {
        "schema_version": "1.0",
        "age_years": 60,
        "sex": "male",
        "weight_kg": 75.0,
        "height_cm": 175.0,
        "serum_creatinine_mg_dL": 1.2,
        "albumin_g_dL": 3.8,
        "hematocrit_percent": 33.0,
        "days_post_transplant": 100,
        "genotypes": {"CYP3A5": "*3/*3"},
        "indication": "test",
    }
    doc.update(over)
    return doc


class TestDialysis(unittest.TestCase):
    def test_haemodialysis_modelled_with_rebound(self):
        doc = base(
            serum_creatinine_mg_dL=6.5,
            weight_kg=70.0,
            renal={
                "dialysis": {
                    "modality": "hd",
                    "sessions_per_week": 3,
                    "hours_per_session": 4.0,
                    "first_session_h": 24.0,
                    "day_interval_h": 48.0,
                    "blood_flow_mL_min": 300,
                    "dialysate_flow_mL_min": 500,
                }
            },
        )
        res = engine.run_dose(doc, "vancomycin", pta_n=200)
        week = res["exposure"].get("dialysis_week")
        self.assertIsNotNone(week, "intradialytic clearance was not applied")
        self.assertEqual(len(week["sessions"]), 3)
        first = week["sessions"][0]
        self.assertLess(first["post_hd"], first["pre_hd"])
        # redistribution after the session raises the concentration again
        self.assertGreater(first["rebound_4h"], 0.0)

    def test_renally_cleared_drug_without_intradialytic_data_is_refused(self):
        # lithium is 100% renally cleared; ignoring dialysis would under-dose by
        # a wide margin, so a warning is not good enough here
        doc = base(serum_creatinine_mg_dL=6.0, renal={"dialysis": {"modality": "hd"}})
        with self.assertRaises(EngineError) as ctx:
            engine.run_dose(doc, "lithium", pta_n=200)
        self.assertIn("intradialytic", str(ctx.exception))


class TestHepatic(unittest.TestCase):
    def test_child_pugh_c_lowers_tacrolimus_clearance(self):
        def cl_for(child_pugh):
            doc = base(
                hepatic={"child_pugh_components": child_pugh},
                serum_creatinine_mg_dL=1.1,
            )
            res = engine.run_dose(doc, "tacrolimus", pta_n=200)
            return res["parameters"]["values"]["cl"]

        a = cl_for(
            {
                "bilirubin_mg_dL": 0.8,
                "inr": 1.0,
                "albumin_g_dL": 4.2,
                "ascites": "none",
                "encephalopathy": "none",
            }
        )
        c = cl_for(
            {
                "bilirubin_mg_dL": 4.0,
                "inr": 2.6,
                "albumin_g_dL": 2.2,
                "ascites": "moderate_severe",
                "encephalopathy": "grade_3_4",
            }
        )
        self.assertAlmostEqual(c / a, 0.5, 2)

    def test_incomplete_child_pugh_warns_and_does_not_scale(self):
        doc = base(hepatic={"child_pugh_components": {"bilirubin_mg_dL": 2.0}})
        res = engine.run_dose(doc, "tacrolimus", pta_n=200)
        codes = {w["code"] for w in res["warnings"]}
        self.assertIn("W_HEPATIC_DATA_MISSING", codes)


class TestObesityAndAge(unittest.TestCase):
    def test_obesity_flagged_and_weight_metrics_differ(self):
        doc = base(weight_kg=140.0, height_cm=170.0, serum_creatinine_mg_dL=1.4)
        res = engine.run_dose(doc, "vancomycin", pta_n=200)
        codes = {w["code"] for w in res["warnings"]}
        self.assertIn("W_OBESITY_WEIGHT_AMBIGUITY", codes)
        crcl = res["patient"]["derived"]["crcl_cg_mL_min"]["value"]
        self.assertGreater(crcl["total_body_weight"], crcl["adjusted_body_weight"])

    def test_extreme_ages_run(self):
        for age in (19, 95):
            res = engine.run_dose(
                base(age_years=age, serum_creatinine_mg_dL=1.1), "vancomycin", pta_n=200
            )
            self.assertEqual(res["status"], "ok")

    def test_missing_height_falls_back_with_warning(self):
        doc = base(height_cm=None, weight_kg=70.0, serum_creatinine_mg_dL=0.9)
        res = engine.run_dose(doc, "lithium", pta_n=200)
        codes = {w["code"] for w in res["warnings"]}
        self.assertIn("W_MODEL_NO_WEIGHT_V", codes)


class TestRefusals(unittest.TestCase):
    def test_unknown_drug(self):
        with self.assertRaises(EngineError) as ctx:
            engine.run_dose(base(), "unobtainium", pta_n=100)
        self.assertEqual(ctx.exception.code, "E_UNKNOWN_DRUG")

    def test_unknown_target(self):
        with self.assertRaises(EngineError) as ctx:
            engine.run_dose(base(), "vancomycin", target_id="made_up", pta_n=100)
        self.assertIn("not declared", str(ctx.exception))

    def test_unresolved_entry_refused(self):
        tmp = Path(tempfile.mkdtemp())
        try:
            src = ROOT / "src" / "precision_dose" / "library"
            shutil.copytree(src, tmp / "library")
            path = tmp / "library" / "drugs" / "vancomycin.toml"
            text = path.read_text()
            path.write_text(text.replace('status = "resolved"', 'status = "unresolved"'))
            with self.assertRaises(EngineError) as ctx:
                engine.run_dose(base(), "vancomycin", pta_n=100, library_dir=tmp / "library")
            self.assertEqual(ctx.exception.code, "E_LIBRARY_UNRESOLVED")
        finally:
            shutil.rmtree(tmp)

    def test_paediatric_refused(self):
        with self.assertRaises(EngineError) as ctx:
            engine.run_dose(base(age_years=8), "vancomycin", pta_n=100)
        self.assertEqual(ctx.exception.code, "E_AGE_OUT_OF_SCOPE")

    def test_missing_model_covariate_refused(self):
        # Andrews 2019 carries albumin and haematocrit as covariates and publishes
        # no population default for them, so their absence is a refusal, not a 1.0
        for field in ("hematocrit_percent", "albumin_g_dL"):
            doc = base()
            doc.pop(field)
            with self.assertRaises(EngineError) as ctx:
                engine.run_dose(doc, "tacrolimus", pta_n=100)
            self.assertEqual(ctx.exception.code, "E_MISSING_FIELD")
            self.assertEqual(ctx.exception.field, field)

    def test_post_transplant_day_no_longer_required(self):
        # Andrews 2019 found time post-transplant not significant; EBE from the
        # previous (Passey) model is gone, so the input is optional now
        doc = base()
        doc.pop("days_post_transplant", None)
        res = engine.run_dose(doc, "tacrolimus", pta_n=100)
        self.assertEqual(res["status"], "ok")


class TestLabelCriteria(unittest.TestCase):
    def test_two_of_three_criteria_reduce_dose(self):
        doc = base(age_years=82, weight_kg=58.0, serum_creatinine_mg_dL=1.6)
        res = engine.run_dose(doc, "apixaban", pta_n=100)
        self.assertEqual(res["regimen"]["dose"], 2.5)
        self.assertEqual(res["regimen"]["decision_basis"], "label_criteria")

    def test_full_dose_when_no_criteria_met(self):
        doc = base(age_years=60, weight_kg=80.0, serum_creatinine_mg_dL=0.9)
        res = engine.run_dose(doc, "apixaban", pta_n=100)
        self.assertEqual(res["regimen"]["dose"], 5.0)

    def test_strong_inhibitor_reduces_dose(self):
        doc = base(
            age_years=60, weight_kg=80.0, serum_creatinine_mg_dL=0.9,
            co_medications=[{"name": "ketoconazole", "dose_mg": 400}],
        )
        res = engine.run_dose(doc, "apixaban", pta_n=100)
        self.assertEqual(res["regimen"]["dose"], 2.5)

    def test_child_pugh_c_not_recommended(self):
        doc = base(
            hepatic={
                "child_pugh_components": {
                    "bilirubin_mg_dL": 4.0,
                    "inr": 2.6,
                    "albumin_g_dL": 2.2,
                    "ascites": "moderate_severe",
                    "encephalopathy": "grade_3_4",
                }
            }
        )
        res = engine.run_dose(doc, "apixaban", pta_n=100)
        self.assertEqual(res["status"], "not_recommended")


if __name__ == "__main__":
    unittest.main()
