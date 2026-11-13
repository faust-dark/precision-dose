"""Determinism: the same input must produce the same bytes.

Checked four ways: repeated calls in one process, a fresh process, a different
working directory, and different PYTHONHASHSEED values. If a dict iteration
order or a hash seed ever leaks into the output, this fails.
"""

import hashlib
import json
import os
import subprocess
import sys
import unittest
from pathlib import Path

from precision_dose import engine

ROOT = Path(__file__).resolve().parents[1]
PATIENT = ROOT / "demo" / "case_vanc.json"


def digest(doc):
    return hashlib.sha256(
        json.dumps(doc, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def run(patient_path, cwd=None):
    code = (
        "import json,sys;from precision_dose import engine;"
        f"p=json.load(open({str(patient_path)!r}));"
        "r=engine.run_dose(p,'vancomycin',pta_n=300);"
        "print(json.dumps(r,sort_keys=True,separators=(',',':')))"
    )
    out = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        cwd=cwd or str(ROOT),
        env={**os.environ},
        check=True,
    )
    return hashlib.sha256(out.stdout.strip().encode()).hexdigest()


class TestDeterminism(unittest.TestCase):
    def test_repeated_calls_identical(self):
        p = json.loads(PATIENT.read_text())
        a = engine.run_dose(p, "vancomycin", pta_n=300)
        b = engine.run_dose(p, "vancomycin", pta_n=300)
        self.assertEqual(digest(a), digest(b))

    def test_seed_changes_result(self):
        p = json.loads(PATIENT.read_text())
        a = engine.run_dose(p, "vancomycin", pta_n=300, seed=1)
        b = engine.run_dose(p, "vancomycin", pta_n=300, seed=2)
        self.assertNotEqual(
            a["pta"]["median_target_value"], b["pta"]["median_target_value"]
        )

    def test_across_processes_and_hash_seeds(self):
        hashes = set()
        for seed in ("0", "1", "random"):
            os.environ["PYTHONHASHSEED"] = seed
            hashes.add(run(PATIENT))
        os.environ.pop("PYTHONHASHSEED", None)
        self.assertEqual(len(hashes), 1, hashes)

    def test_independent_of_working_directory(self):
        a = run(PATIENT, cwd=str(ROOT))
        b = run(PATIENT, cwd="/tmp")
        self.assertEqual(a, b)

    def test_input_key_order_irrelevant(self):
        p = json.loads(PATIENT.read_text())
        shuffled = dict(reversed(list(p.items())))
        a = engine.run_dose(p, "vancomycin", pta_n=300)
        b = engine.run_dose(shuffled, "vancomycin", pta_n=300)
        self.assertEqual(digest(a), digest(b))


if __name__ == "__main__":
    unittest.main()
