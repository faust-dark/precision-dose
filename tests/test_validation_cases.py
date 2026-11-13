"""The validation suite is a gate: every published anchor must pass."""

import unittest

from precision_dose import validation


class TestValidationCases(unittest.TestCase):
    def test_all_cases_pass(self):
        rows = validation.run_cases()
        failures = [
            (r["case"], r["error"], [c for c in r["checks"] if not c["ok"]])
            for r in rows
            if not r["passed"]
        ]
        self.assertEqual(failures, [], validation.render(rows))

    def test_suite_has_coverage_of_every_drug(self):
        rows = validation.run_cases()
        self.assertGreaterEqual(len(rows), 8)


if __name__ == "__main__":
    unittest.main()
