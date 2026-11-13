"""Local HTTP API. Same engine, same determinism, stateless.

Deliberately built on the standard library's ThreadingHTTPServer: no framework,
no dependency, nothing between the request and the calculation.
"""

import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from . import __version__, bayes, clin, engine, library
from .errors import EngineError

_LOCK = threading.Lock()


def _make_handler(library_dir):
    class Handler(BaseHTTPRequestHandler):
        server_version = f"precision-dose/{__version__}"

        def log_message(self, fmt, *args):  # keep stdout clean
            pass

        def _send(self, status, payload):
            body = json.dumps(payload, indent=2, sort_keys=True).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _error(self, status, exc_or_code, message=None, field=None):
            if isinstance(exc_or_code, EngineError):
                payload = {"error": exc_or_code.as_dict()}
            else:
                payload = {"error": {"code": exc_or_code, "message": message}}
                if field:
                    payload["error"]["field"] = field
            self._send(status, payload)

        def _read_body(self):
            length = int(self.headers.get("Content-Length") or 0)
            raw = self.rfile.read(length) if length else b"{}"
            try:
                return json.loads(raw.decode("utf-8"))
            except (UnicodeDecodeError, json.JSONDecodeError) as exc:
                self._error(400, "E_SCHEMA", f"malformed JSON body: {exc}")
                return None

        def do_GET(self):
            path = self.path.split("?")[0].rstrip("/")
            if path in ("/healthz", ""):
                self._send(200, {"status": "ok", "version": __version__})
                return
            if path == "/v1/library":
                lib, _, problems = library.load_library(library_dir)
                self._send(200, {"drugs": library.list_library(lib), "problems": problems})
                return
            if path.startswith("/v1/library/"):
                name = path.rsplit("/", 1)[-1]
                try:
                    lib, _, _ = library.load_library(library_dir)
                    entry = library.get_drug(lib, name)
                except EngineError as exc:
                    self._error(404, exc)
                    return
                clean = {k: v for k, v in entry.items() if not k.startswith("_")}
                self._send(200, clean)
                return
            self._error(404, "E_SCHEMA", "unknown route")

        def do_POST(self):
            path = self.path.split("?")[0].rstrip("/")
            body = self._read_body()
            if body is None:
                return
            try:
                with _LOCK:
                    if path == "/v1/dose":
                        patient = body.get("patient")
                        drug = body.get("drug")
                        if not patient or not drug:
                            self._error(
                                422, "E_MISSING_FIELD", "body needs 'patient' and 'drug'"
                            )
                            return
                        opts = body.get("options") or {}
                        result = engine.run_dose(
                            patient,
                            drug,
                            target_id=body.get("target"),
                            seed=int(opts.get("seed", engine.DEFAULT_SEED)),
                            pta_n=int(opts.get("pta_n", engine.DEFAULT_PTA_N)),
                            library_dir=library_dir,
                            curve_points=int(opts.get("curve_points", 200)),
                        )
                        self._send(200, result)
                        return
                    if path == "/v1/tdm":
                        patient = body.get("patient")
                        drug = body.get("drug")
                        dosing = body.get("dosing") or {
                            "doses": body.get("doses"),
                            "levels": body.get("levels"),
                        }
                        if not patient or not drug or not dosing:
                            self._error(
                                422,
                                "E_MISSING_FIELD",
                                "body needs 'patient', 'drug' and 'dosing'",
                            )
                            return
                        lib, interactions, _ = library.load_library(library_dir)
                        entry = library.get_drug(lib, drug)
                        p = clin.validate_patient(patient)
                        derived, warnings = clin.derive(p)
                        params, _ = engine.individualize(
                            entry, p, derived, interactions, warnings
                        )
                        fit = bayes.map_estimate(entry, params, dosing)
                        iiv_override = {
                            f"{k}_cv_percent": v
                            for k, v in fit["posterior_cv_percent"].items()
                            if v
                        }
                        result = engine.run_dose(
                            patient,
                            drug,
                            seed=int((body.get("options") or {}).get("seed", engine.DEFAULT_SEED)),
                            pta_n=int((body.get("options") or {}).get("pta_n", engine.DEFAULT_PTA_N)),
                            library_dir=library_dir,
                            param_override=fit["posterior"],
                            iiv_override=iiv_override,
                            posterior=fit,
                        )
                        self._send(200, {"tdm": fit, "regimen": result})
                        return
                self._error(404, "E_SCHEMA", "unknown route")
            except EngineError as exc:
                status = 422 if exc.code.startswith(("E_SCHEMA", "E_MISSING", "E_BAD")) else 400
                if exc.code == "E_UNKNOWN_DRUG":
                    status = 404
                self._error(status, exc)

    return Handler


def serve(host="127.0.0.1", port=8765, library_dir=None):
    httpd = ThreadingHTTPServer((host, port), _make_handler(library_dir))
    print(f"precision-dose {__version__} listening on http://{host}:{port}")
    print("POST /v1/dose  POST /v1/tdm  GET /v1/library/{drug}  GET /healthz")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")
    finally:
        httpd.server_close()
