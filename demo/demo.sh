#!/usr/bin/env bash
# End-to-end demonstration: forward dosing, Bayesian revision, label criteria,
# validation suite. Run from the repository root after `pip install -e .`.
set -euo pipefail
cd "$(dirname "$0")/.."

# use the local virtualenv if it exists, otherwise rely on PATH
if [ -x ".venv/bin/precision-dose" ]; then
  PD=".venv/bin/precision-dose"
else
  PD="precision-dose"
fi

echo "=== 1. vancomycin, AUC-guided, renal impairment ==="
"$PD" dose --patient demo/case_vanc.json --drug vancomycin --out demo/out_vanc.json

echo
echo "=== 2. tacrolimus, CYP3A5 non-expresser on fluconazole ==="
"$PD" dose --patient demo/case_tacro.json --drug tacrolimus --pta-n 3000 --out demo/out_tacro.json

echo
echo "=== 3. apixaban, label criteria (two of three) ==="
"$PD" dose --patient demo/case_apixaban.json --drug apixaban --out demo/out_apixaban.json

echo
echo "=== 4. phenytoin, capacity-limited kinetics, CYP2C9 *1/*2 ==="
"$PD" dose --patient demo/case_phenytoin.json --drug phenytoin --out demo/out_phenytoin.json

echo
echo "=== 5. lithium, 12-hour level ==="
"$PD" dose --patient demo/case_lithium.json --drug lithium --out demo/out_lithium.json

echo
echo "=== 6. digoxin, heart failure with renal impairment ==="
"$PD" dose --patient demo/case_digoxin.json --drug digoxin --out demo/out_digoxin.json

echo
echo "=== 7. Bayesian TDM: levels that came back higher than predicted ==="
"$PD" tdm --patient demo/case_vanc.json --drug vancomycin --levels demo/levels_vanc.json

echo
echo "=== 8. validation suite ==="
"$PD" validate
