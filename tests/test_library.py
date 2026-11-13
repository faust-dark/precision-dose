"""Library integrity: the curation contract, enforced.

If a curator adds a parameter without a citation, or leaves an entry unresolved,
this suite fails before any patient is dosed.
"""

import json
import shutil
import tempfile
import unittest
from pathlib import Path

from precision_dose import library
from precision_dose.errors import EngineError

ROOT = Path(__file__).resolve().parents[1]


class TestBundledLibrary(unittest.TestCase):
    def setUp(self):
        self.lib, self.inter, self.problems = library.load_library()

    def test_no_invalid_entries(self):
        self.assertEqual(self.problems, {})

    def test_all_five_drugs_present(self):
        self.assertEqual(
            sorted(self.lib),
            ["apixaban", "lithium", "phenytoin", "tacrolimus", "vancomycin"],
        )

    def test_every_entry_has_citations_and_provenance(self):
        for name, entry in self.lib.items():
            self.assertTrue(entry["citations"], name)
            self.assertTrue(entry["provenance"]["curator"], name)
            self.assertTrue(entry["provenance"]["curated_on"], name)

    def test_every_parameter_source_exists(self):
        for name, entry in self.lib.items():
            sources = entry["model"].get("sources") or {}
            for param, cite in sources.items():
                self.assertIn(cite, entry["citations"], f"{name}:{param}")

    def test_targets_and_safety_have_sources(self):
        for name, entry in self.lib.items():
            for t in entry["targets"]:
                self.assertIn(t["source"], entry["citations"], f"{name}:{t['id']}")
            for s in entry["safety"]:
                self.assertIn(s["source"], entry["citations"], f"{name}:{s['id']}")

    def test_hashes_are_stable(self):
        again, _, _ = library.load_library()
        for name in self.lib:
            self.assertEqual(self.lib[name]["content_hash"], again[name]["content_hash"])

    def test_interaction_perpetrators_have_sources(self):
        for name, perp in self.inter["perpetrators"].items():
            self.assertIn(perp["source"], self.inter["citations"], name)

    def test_tacrolimus_requires_its_model_covariates(self):
        requires = self.lib["tacrolimus"]["model"]["requires"]
        self.assertIn("bsa_m2", requires)
        self.assertIn("albumin_g_dL", requires)
        self.assertIn("hematocrit_percent", requires)

    def test_every_parameter_names_one_model(self):
        # the tacrolimus entry was rebuilt on Andrews 2019 precisely so that no
        # parameter is borrowed from a second cohort
        src = self.lib["tacrolimus"]["model"]["sources"]
        for param in ("cl_ref", "v1_ref", "v2_ref", "q_ref", "ka_per_h", "iiv", "residual"):
            self.assertEqual(src[param], "S1", param)


class TestValidationContract(unittest.TestCase):
    def test_missing_source_is_rejected(self):
        tmp = Path(tempfile.mkdtemp())
        try:
            src = ROOT / "src" / "precision_dose" / "library"
            shutil.copytree(src, tmp / "library")
            path = tmp / "library" / "drugs" / "vancomycin.toml"
            text = path.read_text().replace('cl_ref = "S1"', 'cl_ref = "S9"')
            path.write_text(text)
            lib, _, problems = library.load_library(tmp / "library")
            self.assertIn("vancomycin", problems)
            self.assertTrue(any("S9" in p for p in problems["vancomycin"]))
        finally:
            shutil.rmtree(tmp)

    def test_placeholder_value_is_rejected(self):
        tmp = Path(tempfile.mkdtemp())
        try:
            src = ROOT / "src" / "precision_dose" / "library"
            shutil.copytree(src, tmp / "library")
            path = tmp / "library" / "drugs" / "lithium.toml"
            path.write_text(path.read_text().replace("km_mg_L", "x", 0))
            text = path.read_text().replace(
                'curated_on = "2026-11-13"', 'curated_on = "PLACEHOLDER"'
            )
            path.write_text(text)
            _, _, problems = library.load_library(tmp / "library")
            self.assertIn("lithium", problems)
        finally:
            shutil.rmtree(tmp)


if __name__ == "__main__":
    unittest.main()
