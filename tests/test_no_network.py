"""No network. Not a promise in a README -- a test that fails if the engine so
much as opens a socket."""

import ast
import importlib
import pkgutil
import socket
import sys
import unittest

from precision_dose import engine

PATIENT = {
    "schema_version": "1.0",
    "age_years": 68,
    "sex": "male",
    "weight_kg": 78.0,
    "height_cm": 172.0,
    "serum_creatinine_mg_dL": 1.9,
    "albumin_g_dL": 2.8,
    "indication": "mrsa_bacteremia",
    "mic_mg_L": 1.0,
}

# the complete set of modules the engine is allowed to import
ALLOWED = {
    "precision_dose", "json", "math", "hashlib", "tomllib", "argparse", "http",
    "threading", "pathlib", "sys", "random", "collections", "ast",
}


class TestNoNetwork(unittest.TestCase):
    def test_engine_runs_with_sockets_disabled(self):
        original = socket.socket
        original_create = socket.create_connection

        class Blocked(socket.socket):
            def __init__(self, *a, **k):
                raise AssertionError("the engine attempted to open a socket")

        socket.socket = Blocked
        socket.create_connection = lambda *a, **k: (_ for _ in ()).throw(
            AssertionError("the engine attempted a network connection")
        )
        try:
            res = engine.run_dose(PATIENT, "vancomycin", pta_n=300)
        finally:
            socket.socket = original
            socket.create_connection = original_create
        self.assertEqual(res["status"], "ok")
        self.assertTrue(res["rationale"])

    def test_no_third_party_imports_in_the_engine(self):
        """Parsed with ast, not grepped: an import inside a function counts too."""
        import precision_dose

        offenders = []
        for mod in pkgutil.walk_packages(
            precision_dose.__path__, prefix="precision_dose."
        ):
            module = importlib.import_module(mod.name)
            path = getattr(module, "__file__", None)
            if not path or not path.endswith(".py"):
                continue
            with open(path, "r", encoding="utf-8") as fh:
                tree = ast.parse(fh.read())
            for node in ast.walk(tree):
                roots = []
                if isinstance(node, ast.Import):
                    roots = [a.name.split(".")[0] for a in node.names]
                elif isinstance(node, ast.ImportFrom) and node.level == 0:
                    roots = [(node.module or "").split(".")[0]]
                for root in roots:
                    if root and root not in ALLOWED:
                        offenders.append(f"{mod.name} imports {root}")
        self.assertEqual(offenders, [])


if __name__ == "__main__":
    unittest.main()
