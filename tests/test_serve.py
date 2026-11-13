"""The HTTP API is a transport, not a second implementation.

The test starts the real server on an ephemeral port and asserts that the JSON
it returns is byte-identical to the in-process engine result for the same input.
"""

import json
import threading
import unittest
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

from precision_dose import engine, serve

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


class TestServe(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.httpd = ThreadingHTTPServer(("127.0.0.1", 0), serve._make_handler(None))
        cls.port = cls.httpd.server_address[1]
        cls.thread = threading.Thread(target=cls.httpd.serve_forever, daemon=True)
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.httpd.shutdown()
        cls.httpd.server_close()

    def post(self, path, payload):
        req = urllib.request.Request(
            f"http://127.0.0.1:{self.port}{path}",
            data=json.dumps(payload).encode(),
            headers={"Content-Type": "application/json"},
        )
        with urllib.request.urlopen(req) as r:
            return json.load(r)

    def test_healthz(self):
        with urllib.request.urlopen(f"http://127.0.0.1:{self.port}/healthz") as r:
            self.assertEqual(json.load(r)["status"], "ok")

    def test_dose_matches_cli_engine(self):
        body = {"patient": PATIENT, "drug": "vancomycin",
                "options": {"pta_n": 400, "seed": 7}}
        via_http = self.post("/v1/dose", body)
        direct = engine.run_dose(PATIENT, "vancomycin", pta_n=400, seed=7)
        self.assertEqual(
            json.dumps(via_http, sort_keys=True), json.dumps(direct, sort_keys=True)
        )

    def test_unknown_drug_is_404(self):
        with self.assertRaises(urllib.error.HTTPError) as ctx:
            self.post("/v1/dose", {"patient": PATIENT, "drug": "unobtainium"})
        self.assertEqual(ctx.exception.code, 404)

    def test_bad_body_is_422(self):
        with self.assertRaises(urllib.error.HTTPError) as ctx:
            self.post("/v1/dose", {"patient": PATIENT})
        self.assertEqual(ctx.exception.code, 422)

    def test_library_endpoint(self):
        with urllib.request.urlopen(f"http://127.0.0.1:{self.port}/v1/library/vancomycin") as r:
            entry = json.load(r)
        self.assertEqual(entry["drug"], "vancomycin")
        self.assertIn("content_hash", entry)


if __name__ == "__main__":
    unittest.main()
