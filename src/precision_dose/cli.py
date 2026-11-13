"""Command line interface.

    precision-dose dose     --patient p.json --drug vancomycin [--target ...]
    precision-dose tdm      --patient p.json --drug vancomycin --levels l.json
    precision-dose library  list|show|validate
    precision-dose validate [--cases DIR]
    precision-dose serve    [--host 127.0.0.1] [--port 8765]

Exit codes: 0 ok, 2 invalid input, 3 refusal / insufficient data, 4 target
unattainable (best candidate returned, never hidden).
"""

import argparse
import json
import sys
from pathlib import Path

from . import __version__, bayes, engine, library, validation
from .errors import EngineError

EXIT_INVALID = 2
EXIT_REFUSED = 3
EXIT_UNATTAINABLE = 4


def _read_json(path):
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except FileNotFoundError:
        raise EngineError("E_SCHEMA", f"file not found: {path}")
    except json.JSONDecodeError as exc:
        raise EngineError("E_SCHEMA", f"{path} is not valid JSON: {exc}")


def _dump(obj, path):
    text = json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False)
    if path:
        Path(path).write_text(text + "\n", encoding="utf-8")
        print(f"wrote {path}", file=sys.stderr)
    else:
        print(text)


def _exit_code(exc):
    code = getattr(exc, "code", "")
    if code in ("E_LIBRARY_UNRESOLVED", "E_NOT_CONVERGED", "E_UNSUPPORTED", "E_LIBRARY_NO_MODEL"):
        return EXIT_REFUSED
    return EXIT_INVALID


def cmd_dose(args):
    patient_doc = _read_json(args.patient)
    result = engine.run_dose(
        patient_doc,
        args.drug,
        target_id=args.target,
        seed=args.seed,
        pta_n=args.pta_n,
        library_dir=args.library,
        curve_points=args.curve_points,
    )
    if args.out:
        _dump(result, args.out)
        if not args.quiet:
            print(render_text(result))
    elif args.format == "json":
        _dump(result, None)
    else:
        print(render_text(result))
        if args.trace:
            print()
            print("\n".join(result.get("trace") or []))
    return EXIT_UNATTAINABLE if result.get("status") == "target_unattainable" else 0


def cmd_tdm(args):
    patient_doc = _read_json(args.patient)
    levels_doc = _read_json(args.levels)
    lib, interactions, _ = library.load_library(args.library)
    entry = library.get_drug(lib, args.drug)

    from . import clin

    patient = clin.validate_patient(patient_doc)
    derived, warnings = clin.derive(patient)
    params, contribs = engine.individualize(entry, patient, derived, interactions, warnings)
    fit = bayes.map_estimate(entry, params, levels_doc)

    override = dict(fit["posterior"])
    iiv_override = {
        f"{k}_cv_percent": v for k, v in fit["posterior_cv_percent"].items() if v
    }
    redose = engine.run_dose(
        patient_doc,
        args.drug,
        target_id=args.target,
        seed=args.seed,
        pta_n=args.pta_n,
        library_dir=args.library,
        param_override=override,
        iiv_override=iiv_override,
        posterior=fit,
    )
    out = {"tdm": fit, "regimen": redose}
    if args.out:
        _dump(out, args.out)
    elif args.format == "json":
        _dump(out, None)
    else:
        print(render_tdm_text(out))
    return 0


def cmd_library(args):
    lib, _, problems = library.load_library(args.library)
    if args.action == "list":
        rows = library.list_library(lib)
        if args.format == "json":
            _dump(rows, None)
        else:
            for r in rows:
                print(
                    f"{r['drug']:<12} {r['status']:<10} {str(r['structure']):<18} "
                    f"{r['decision_basis']:<14} targets={','.join(r['targets'])}"
                )
        return 0
    if args.action == "show":
        entry = lib.get(args.drug)
        if entry is None:
            raise EngineError("E_UNKNOWN_DRUG", f"unknown drug {args.drug!r}")
        clean = {k: v for k, v in entry.items() if not k.startswith("_")}
        _dump(clean, None)
        return 0
    # validate
    ok = True
    for name, probs in sorted(problems.items()):
        ok = False
        print(f"FAIL {name}")
        for p in probs:
            print(f"     - {p}")
    for name in sorted(lib):
        if name in problems:
            continue
        e = lib[name]
        print(
            f"OK   {name:<12} status={e['status']:<10} "
            f"hash={e['content_hash'][:19]}... targets={len(e['targets'])}"
        )
    return 0 if ok else EXIT_INVALID


def cmd_validate(args):
    rows = validation.run_cases(args.cases, library_dir=args.library)
    failed = [r for r in rows if not r["passed"]]
    if args.format == "json":
        _dump(rows, None)
    else:
        print(validation.render(rows))
    return EXIT_REFUSED if failed else 0


def cmd_serve(args):
    from .serve import serve

    serve(host=args.host, port=args.port, library_dir=args.library)
    return 0


# --------------------------------------------------------------------------
# text rendering
# --------------------------------------------------------------------------


def render_text(result):
    if result.get("exposure") is None:
        return result["rationale"]
    reg = result["regimen"]
    m = reg["maintenance"]
    lines = [
        result["rationale"],
        "",
        "TOP CANDIDATES (common random numbers, same seed)",
        f"{'dose':>10} {'interval':>9} {'PTA':>6} {'median target':>14}",
    ]
    for row in result["pta"]["table"][:8]:
        lines.append(
            f"{row['dose']:>10g} {row['interval_h']:>9g} {row['pta']*100:>5.1f}% "
            f"{row['median_target']:>14.3f}"
        )
    lines.append("")
    lines.append(
        f"SELECTED: {m['dose']:g} {m['unit']} every {m['interval_h']:g} h"
        + (f" over {m['infusion_min']} min" if m.get("infusion_min") else "")
        + f"  (PTA {result['pta']['probability']*100:.1f}%, CI95 "
        f"{result['pta']['ci95'][0]*100:.0f}-{result['pta']['ci95'][1]*100:.0f}%)"
    )
    return "\n".join(lines)


def render_tdm_text(out):
    fit = out["tdm"]
    reg = out["regimen"]
    lines = ["BAYESIAN TDM", "=" * 12]
    lines.append(f"OFV {fit['ofv']}; levels {fit['n_levels']}; converged {fit['converged']}")
    for k in fit["free_parameters"]:
        lines.append(
            f"  {k.upper():<4} prior {fit['prior'][k]:.3f} -> posterior "
            f"{fit['posterior'][k]:.3f} (CV {fit.get('posterior_cv_percent', {}).get(k)}%"
            f", information ratio {fit['information_ratio'].get(k)})"
        )
    lines.append("")
    lines.append(f"{'time':>8} {'observed':>9} {'predicted':>10} {'wres':>8} {'err %':>8}")
    for l in fit["levels"]:
        lines.append(
            f"{l['time_h']:>8.2f} {l['observed']:>9.3f} {l['predicted']:>10.3f} "
            f"{l['weighted_residual']:>8.2f} {l['percent_error']:>8.1f}"
        )
    for w in fit.get("warnings", []):
        lines.append(f"  ! {w['message']}")
    lines.append("")
    lines.append("RE-DOSING WITH POSTERIOR PARAMETERS")
    lines.append(reg["rationale"].split("Caveats:")[0].rstrip())
    return "\n".join(lines)


# --------------------------------------------------------------------------
# parser
# --------------------------------------------------------------------------


def build_parser():
    p = argparse.ArgumentParser(
        prog="precision-dose",
        description="Deterministic model-informed precision dosing for one patient.",
    )
    p.add_argument("--version", action="version", version=f"precision-dose {__version__}")
    p.add_argument("--library", default=None, help="library directory (default: bundled)")
    sub = p.add_subparsers(dest="command", required=True)

    d = sub.add_parser("dose", help="compute a regimen")
    d.add_argument("--patient", required=True)
    d.add_argument("--drug", required=True)
    d.add_argument("--target", default=None)
    d.add_argument("--seed", type=int, default=engine.DEFAULT_SEED)
    d.add_argument("--pta-n", type=int, default=engine.DEFAULT_PTA_N)
    d.add_argument("--curve-points", type=int, default=200)
    d.add_argument("--out", default=None)
    d.add_argument("--format", choices=("text", "json"), default="text")
    d.add_argument("--trace", action="store_true")
    d.add_argument("--quiet", action="store_true")
    d.set_defaults(func=cmd_dose)

    t = sub.add_parser("tdm", help="Bayesian dose revision from measured levels")
    t.add_argument("--patient", required=True)
    t.add_argument("--drug", required=True)
    t.add_argument("--levels", required=True)
    t.add_argument("--target", default=None)
    t.add_argument("--seed", type=int, default=engine.DEFAULT_SEED)
    t.add_argument("--pta-n", type=int, default=engine.DEFAULT_PTA_N)
    t.add_argument("--out", default=None)
    t.add_argument("--format", choices=("text", "json"), default="text")
    t.set_defaults(func=cmd_tdm)

    l = sub.add_parser("library", help="inspect and validate the drug library")
    l.add_argument("action", choices=("list", "show", "validate"))
    l.add_argument("drug", nargs="?", default=None)
    l.add_argument("--format", choices=("text", "json"), default="text")
    l.set_defaults(func=cmd_library)

    v = sub.add_parser("validate", help="run the published-regimen validation suite")
    v.add_argument("--cases", default=None)
    v.add_argument("--format", choices=("text", "json"), default="text")
    v.set_defaults(func=cmd_validate)

    s = sub.add_parser("serve", help="local HTTP API")
    s.add_argument("--host", default="127.0.0.1")
    s.add_argument("--port", type=int, default=8765)
    s.set_defaults(func=cmd_serve)
    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except EngineError as exc:
        payload = {"error": exc.as_dict()}
        print(json.dumps(payload, indent=2, sort_keys=True), file=sys.stderr)
        return _exit_code(exc)


if __name__ == "__main__":
    sys.exit(main())
