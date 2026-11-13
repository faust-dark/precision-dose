# precision-dose

Deterministic, mechanistic, auditable **model-informed precision dosing (MIPD)** for one
patient at a time. Python, no runtime dependencies, no network, no LLM in the decision
path.

Most dosing is still one-size-fits-all, while a patient's pathology changes drug
exposure: renal impairment, hepatic impairment, age and body size, CYP450 genotype,
protein-binding changes and polypharmacy. Model-informed precision dosing exists in
transplant, oncology and infectious disease, but the tooling is closed and hard to
audit. `precision-dose` is a small engine that a pharmacist can read, run and defend.

```
precision-dose dose --patient patient.json --drug vancomycin --target auc24 --out regimen.json
```

returns the regimen, the predicted concentration-time curve, the **probability of target
attainment (PTA)**, the safety margin, and a pharmacist-readable rationale — with the
population pharmacokinetics citation and the content hash of the exact library entry
that produced it.

## What this is, and what it is not

It is a **clinical decision support (CDS)** calculation: a mechanistic pharmacokinetic
model, individualised with the patient's covariates, evaluated against a guideline
target, with its uncertainty exposed. It is not a prescribing system, not an "AI
therapeutics" demo, and not a platform. One patient, one question, one answer, in
seconds, reproducibly.

Three commitments make it auditable rather than merely accurate-looking:

1. **Every number is traceable.** Each parameter in the drug library names a citation;
   each derived quantity in the output names the formula that produced it; each result
   carries the library entry hash.
2. **Missing data is refused, never imputed.** No silent defaults, no guessed allele
   activity scores, no interpolation of a Child-Pugh class.
3. **The model's own algebra is exposed.** Covariate effects multiply, so their
   log-contributions add exactly; the attribution output is a decomposition, not a
   heuristic ranking, and the residual of the decomposition is printed.

## A real case, run by this code

A 68-year-old man, 78 kg, 172 cm, serum creatinine 1.9 mg/dL, MRSA bacteraemia
(MIC 1.0 mg/L), albumin 2.8 g/dL. The question: what vancomycin regimen achieves
AUC24/MIC 400-600?

```
$ precision-dose dose --patient demo/case_vanc.json --drug vancomycin --out regimen.json

Question: Vancomycin regimen for a 68-year-old male, 78 kg, SCr 1.90 mg/dL, MIC 1.00 mg/L, target auc24_mic (400-600 h).
Renal function: CKD-EPI 2021 eGFR 38 mL/min/1.73 m2; Cockcroft-Gault CrCl 41 mL/min (TBW), 36 (IBW), 38 (ABW). This entry scales clearance with crcl_cg_total. Discarding one of these in favour of the other is a decision with clinical consequences, which is why both are shown.
Hepatic: albumin 2.8 g/dL; Child-Pugh class A (score 6); clearance multiplier 1.00 (not applied for this drug).
Individual parameters: CL 1.81 L/h, V1 36.86 L, V2 44.20 L, Q 6.95 L/h; terminal half-life 33.7 h.
Clearance: 4.40 -> 1.81 L/h (-59.0%) because renal_crcl_cg_total x0.41 [crcl_cg_total 41.1 vs 100].
Volume: 37.94 -> 36.86 L (-2.9%) because age x0.97 [age 68 y vs 70 y].
Contributions multiply; do not add the percentages. Log-space residual 0.0e+00 confirms the decomposition is exact.
Regimen: loading 2000 mg over 120 min (25 mg/kg of total body weight), maintenance 1750 mg every 48 h over 120 min.
Predicted steady-state exposure: AUC24 484.4 mg/L*h, Cmax 49.14, trough 11.10 mg/L (90% population interval 273.0-872.3).
Probability of target attainment: 43% (95% CI 42-44%, n=5000, seed 20260101). The spread is between-patient variability, not numerical noise.
Safety: predicted auc24 484.4 vs threshold 600; margin x1.24; P(exceed) 26.7%.
Caveats:
  - no fixed regimen attains auc24_mic in 80% of patients (best 43%): the target window is narrow relative to between-patient variability. This is the quantitative case for therapeutic drug monitoring rather than a reason to widen the target.

TOP CANDIDATES (common random numbers, same seed)
      dose  interval    PTA  median target
      1750        48  43.0%        483.551
       500        12  40.3%        552.630
      1000        24  40.3%        552.630
```

Three things in that output are the point of the tool:

- **43% PTA is the honest answer.** A fixed dose cannot put 80% of a population inside a
  400-600 window when clearance varies by 36% CV. The tool says so and names the
  consequence: TDM, not a wider target. `P(exceed) 26.7%` is the other half of that
  sentence.
- **CrCl 41 (TBW), 36 (IBW), 38 (ABW)** are all printed because the guidelines disagree
  about which weight to use, and the library declares which one its model was estimated
  against. CKD-EPI 2021 eGFR is shown alongside, because it is the right number for CKD
  staging and the wrong one for most drug labels.
- **`1750 mg q48h` beats `1000 mg q24h` on PTA** (43.0% vs 40.3%) because it centres the
  median AUC closer to the middle of the window. The full candidate table is in the JSON,
  so the choice can be argued with rather than merely trusted.

### Adding measured levels: Bayesian dose optimisation

Three doses of 1000 mg q12h later, a peak of 42.0 mg/L at 26 h and a trough of 28.5 mg/L
at 36 h. The population model predicted 39.7 and 27.3; the patient is clearing the drug
more slowly than the covariates suggested.

```
$ precision-dose tdm --patient demo/case_vanc.json --drug vancomycin --levels demo/levels_vanc.json

BAYESIAN TDM
============
OFV 1.739796; levels 2; converged True
  CL   prior 1.806 -> posterior 1.382 (CV 17.8947%, information ratio 0.5223)
  V1   prior 36.856 -> posterior 33.395 (CV 0.5305%, information ratio 0.2836)
  V2   prior 44.200 -> posterior 27.989 (CV 1.4259%, information ratio 0.5144)

    time  observed  predicted     wres    err %
   26.00    42.000     39.718     0.51     -5.4
   36.00    28.500     27.349     0.26     -4.0
...
Clearance: 4.40 -> 1.38 L/h (-68.6%) because renal_crcl_cg_total x0.41 [crcl_cg_total 41.1 vs 100]; tdm_posterior x0.77 [MAP estimate from the measured levels].
Regimen: loading 2000 mg over 120 min (25 mg/kg of total body weight), maintenance 1250 mg every 48 h over 120 min.
Probability of target attainment: 39% (95% CI 37-42%, n=1500, seed 20260101).
```

The **information ratio** (1 − posterior variance / prior variance) is the number most
Bayesian tools hide: 0.52 means the two levels moved the clearance estimate about half
way from the population prior to the truth. It is not a confidence interval pretending
that two samples are a dataset. The re-dosing PTA then integrates the posterior
uncertainty as well as between-patient variability, and the attribution folds the
posterior in as its own exact factor (`tdm_posterior x0.77`) so the breakdown still
multiplies out to the clearance actually used.

## Install and run

```bash
python -m venv .venv && . .venv/bin/activate
pip install -e .            # no runtime dependencies
precision-dose library list
precision-dose validate     # published-regimen anchors, exits non-zero on drift
```

Commands:

| command | what it does |
| --- | --- |
| `precision-dose dose --patient p.json --drug vancomycin [--target auc24_mic] [--out r.json] [--trace]` | regimen, curve, PTA, safety margin, attribution, rationale |
| `precision-dose tdm --patient p.json --drug vancomycin --levels l.json` | MAP estimate from measured levels, then re-dose |
| `precision-dose library list\|show\|validate` | inspect and check the curated library |
| `precision-dose validate` | run the validation anchors |
| `precision-dose serve` | the same engine over a local HTTP API |

Exit codes: `0` ok, `2` invalid input, `3` refusal, `4` target unattainable (the best
candidate is still returned, never hidden).

`--trace` prints the calculation: every formula with the patient's numbers substituted.

### HTTP API

```bash
precision-dose serve --host 127.0.0.1 --port 8765
curl -s localhost:8765/v1/dose -H 'content-type: application/json' \
  -d '{"patient": {...}, "drug": "vancomycin", "options": {"pta_n": 5000}}'
```

`POST /v1/dose`, `POST /v1/tdm`, `GET /v1/library/{drug}`, `GET /healthz`. Stateless, no
storage, no auth. The JSON body is byte-identical to the CLI output for the same input —
there is one engine, and the HTTP layer is a transport.

## The model

Everything is mechanistic; the only fitted component is the population parameter set in
the library entry, which comes from a published population pharmacokinetic model.

```
CL = CL_ref · (renal/renal_ref)^θ · (WT/70)^0.75 · CP · PGx · DDI · fu
V  = V_ref  · (WT/70)^1.0   · CP
```

- **Structural models**: one- and two-compartment linear models plus phenytoin's
  Michaelis-Menten elimination. Linear systems are solved *exactly* with the matrix
  exponential, and a periodic regimen is evaluated at its fixed point with one linear
  solve: no step-size error, mass balance to machine precision, and no dependence on a
  sampling grid. The biexponential closed form is a unit test.
- **Renal function**: eGFR by CKD-EPI 2021 and Cockcroft-Gault CrCl with three weight
  metrics. Which one scales clearance is declared by the drug entry, because guidelines
  differ and the choice changes the answer.
- **Hepatic impairment**: Child-Pugh class scales clearance only where the entry supplies
  a cited per-drug multiplier. MELD is not used: it predicts mortality, not metabolic
  capacity. An incomplete Child-Pugh means no scaling plus a warning.
- **Pharmacogenomics (CPIC)**: diplotype → activity score → phenotype, from versioned CPIC
  tables. An allele without a consensus score is refused rather than defaulted to 1.0.
  Effects are applied to the fraction of clearance the enzyme accounts for, or — for
  phenytoin — to Vmax, which reproduces the CPIC 25%/50% dose guidance as a consequence of
  the mechanism rather than as a lookup.
- **Drug-drug interactions (DDI)**: competitive (Ki and unbound concentration),
  mechanism-based inactivation (kinact/kdeg) and induction (Emax) are supported and are
  *not* interchangeable. Effects flow through the fraction metabolised, so an inhibitor of
  a pathway the drug does not use changes nothing — and the output says so.
- **Binding**: free fraction from albumin where the entry declares it matters, with the
  free-drug-is-active distinction. Vancomycin's TDM target is a *total* concentration and
  the engine refuses to free-correct it.
- **Dialysis**: intradialytic clearance declared by the entry, with the between-session
  rebound emerging from the two-compartment model rather than being hard-coded. A renally
  cleared drug without an intradialytic clearance is refused for dialysis patients.

### PTA and regimen selection

Monte Carlo over the library's between-subject variability, fixed seed, **common random
numbers** across candidates (so the ranking is not Monte Carlo noise), Wilson 95%
intervals, and a documented inverse-normal transform (Wichura AS241) so results are
stable across Python builds. Candidates come from the entry's constraints — strengths,
permitted intervals, maximum single and daily doses — not from a free grid, and the ones
that were excluded are listed.

### Determinism

Same input, same bytes. Proven by tests across repeated runs, fresh processes, different
working directories, permuted input key order, and `PYTHONHASHSEED=0/1/random`. Output has
sorted keys, quantised floats, no timestamps, no hostnames, no locale-dependent
formatting.

## The drug library

Five entries, each chosen to exercise a different mechanism:

| drug | mechanism exercised | decision basis |
| --- | --- | --- |
| vancomycin | two-compartment, renal scaling, AUC24/MIC target | PK target |
| tacrolimus | hepatic CYP3A4/3A5, whole-blood trough, azole DDI | PK target |
| apixaban | label criteria, not clearance | label criteria |
| phenytoin | Michaelis-Menten saturation, CYP2C9, salt factor | PK target |
| lithium | purely renal, mmol units, 12-h monitoring level | PK target |

`apixaban` is included deliberately: its dose is set by three patient criteria and an
interaction class, not by a target exposure, so the tool reports criteria and refuses to
invent a concentration curve. Making that distinction visible is part of the point.

Each entry carries its citations, the target definition and its source, the safety
threshold, the monitoring strategy, the curator, the review date, and an honest
`known_gaps` list that is printed in every rationale. Entries whose numbers cannot be
sourced ship as `status = "unresolved"` and the engine refuses to dose them.

**Curation is enforced by code.** `precision-dose library validate` rejects an entry with
a parameter that names a citation that does not exist, a target without a source, a
non-positive variance, a fraction-metabolised sum above 1, or a leftover placeholder.
Add a drug by following `CURATION.md`.

## Validation

`precision-dose validate` reproduces published anchors and exits non-zero on drift:

```
ckd_epi_2021_reference_case        PASS       2/2    ok
cockcroft_gault_reference_case     PASS       3/3    ok
vancomycin_ashp_2020_auc_target    PASS       4/4    ok
apixaban_label_dose_reduction      PASS       5/5    ok
phenytoin_michaelis_menten         PASS       3/3    ok
lithium_12h_level                  PASS       3/3    ok
tacrolimus_cyp3a5_and_azole_ddi    PASS       4/4    ok
paediatric_refusal                 PASS       0/0    E_AGE_OUT_OF_SCOPE
8/8 cases passed
```

Plus 81 unit and property tests: analytic PK oracles, mass balance, monotonicity (AUC in
dose, CL in renal function, CL in inhibitor strength, V in weight), Bayesian parameter
recovery, determinism, no-network, and the edge cases (dialysis, severe hepatic
impairment, obesity, age 19 and 95, missing height, anuria).

## Limitations, stated plainly

- **Adults only.** Paediatric allometry and maturation are not implemented, so patients
  under 18 are refused rather than dosed with adult exponents.
- **One model per drug.** Population models disagree; each entry names the one it uses and
  lists the disagreement in `known_gaps`. Swapping the model is a curation decision, not a
  code change.
- **Some parameters are declared assumptions.** Tacrolimus IIV and absorption rate are
  assumptions, marked as such and bounded by tests; they widen the PTA but do not move the
  point estimate. The entry says so in every output.
- **No FHIR server client.** The input schema is deliberately flat so an HL7 FHIR
  `Patient`/`Observation`/`MedicationRequest` bundle maps onto it without a database; the
  adapter is the next milestone, not a shipped feature.
- **Not a medical device.** It computes and explains; a pharmacist decides. The output
  carries that statement, and the license says the same thing in fewer words.

## Repository layout

```
src/precision_dose/
  clin.py          patient validation, body size, eGFR/CrCl, Child-Pugh
  pgx.py           CPIC activity scores -> phenotype -> multiplier
  ddi.py           competitive / mechanism-based / induction algebra
  pk.py            exact linear compartment engine, Michaelis-Menten, dialysis timeline
  engine.py        individualisation, PTA, selection, output assembly
  bayes.py         MAP estimation, posterior covariance, information ratio
  attribution.py   exact log-space covariate decomposition
  report.py        deterministic rationale and calculation trace
  library/         curated TOML entries + interaction table (schema-validated)
  cli.py, serve.py, validation.py
validation/cases/  published anchors as JSON
tests/             unit and property tests, standard library only
demo/              the cases used above
```

## Keywords

model-informed precision dosing (MIPD) · Bayesian dose optimization · therapeutic drug
monitoring (TDM) · pharmacometrics · population pharmacokinetics · clinical decision
support (CDS) · HL7 FHIR · CYP450 · eGFR · Child-Pugh · pharmacogenomics (CPIC) ·
drug-drug interactions (DDI) · probability of target attainment (PTA) · vancomycin ·
tacrolimus · phenytoin · lithium · renal impairment · dose optimisation · clinical
pharmacology · precision medicine

## License

MIT. See `LICENSE` for the licence and `DISCLAIMER.md` for the clinical scope: this is a
computational aid for qualified clinicians, not a medical device, and it ships with no
warranty of fitness for clinical use.
