# precision-dose — design proposal (v1, for approval)

Status: PROPOSED. No code written. This document is the review surface.
Author of record for clinical content: the pharmacist/curator.
Engine author: implementation.

---

## 1. Purpose, scope, non-goals

Purpose: a deterministic, mechanistic, auditable model-informed precision dosing (MIPD)
engine for one patient, one drug, one answer, in seconds, with the reasoning and the
uncertainty attached to every number.

In scope (v1)
- Single-patient forward dosing with covariate-informed individualization.
- Regimen selection against a library-declared target exposure.
- Monte Carlo probability of target attainment (PTA) with a fixed seed.
- Safety margin against a library-declared toxicity threshold.
- Bayesian TDM: MAP estimation from measured levels, then re-dosing.
- Covariate attribution in pharmacist language.
- Curated, cited, versioned drug library; curation enforced by code.
- Reproducible validation against published regimens/guidelines.

Non-goals (v1)
- Not a platform: no batch mode, no job queue, no database, no user accounts.
- No network access at runtime, no LLM anywhere in the decision path.
- No pediatric/neonatal dosing (see 5.2 refusal policy) — maturation models are a v2 item.
- Not a prescribing system: it computes and explains; the pharmacist decides.
  Output carries an explicit CDS statement to that effect.

## 2. Decisions on the table (defaults I will implement unless you override)

D1. v1 drug library = all five: vancomycin, tacrolimus, apixaban, phenytoin, lithium.
    Rationale: each exercises a distinct mechanism — 2-cmt renal with AUC target;
    hepatic CYP3A4/CYP3A5 + P-gp with trough target; DOAC label-criteria decision
    (not CL-driven); nonlinear Michaelis-Menten with binding correction; purely renal
    narrow-therapeutic-index with dialysis. If any one of them cannot be sourced to a
    primary publication, its entry ships UNRESOLVED and the tool refuses to dose it
    (never a guessed number).

D2. Bayesian TDM is in v1, scoped: 1- and 2-compartment, library priors, deterministic
    optimizer with a convergence gate, posterior covariance, prior->posterior
    information ratio, and a hard refusal if the optimizer does not converge.

D3. The library owns target definitions. Each entry declares target type, range, unit,
    matrix (serum vs whole blood), the guideline source, and the valid `--target` ids.
    The CLI can only select a declared target; it cannot invent a range. Default target
    per drug is the entry's `default_target`.

D4. Library lives in the package as YAML (`src/precision_dose/library/drugs/*.yaml`),
    overridable with `--library DIR`. Every entry carries schema_version, citations,
    curator provenance, and a content hash (sha256 over canonicalized entry) printed
    in every output so a result can be traced to the exact data that produced it.

D5. Zero runtime dependencies (pure stdlib). Analytic/exact linear solutions, RK4 for
    the nonlinear model, in-house deterministic optimizer, stdlib `http.server` for
    `serve`. numpy/scipy appear only in a dev-only test extra, used to cross-check the
    hand-rolled numerics (never imported by the engine).

If you prefer different answers, say which and I will revise before implementing.

## 3. Repo layout

precision-dose/
  pyproject.toml            hatchling; requires-python >=3.11; no runtime deps
  LICENSE                   MIT
  README.md                 positioning, worked case (real output), keyword set
  DESIGN.md                 this document
  CURATION.md               how a pharmacist adds/updates a drug entry (checklist)
  CHANGELOG.md
  src/precision_dose/
    __init__.py             __version__
    __main__.py             python -m precision_dose
    cli.py                  argparse; subcommands: dose, tdm, library, validate, serve
    schemas.py              dataclasses + strict validators; error codes
    io.py                   JSON load/save; canonical output; deterministic formatting
    units.py                unit conversion + salt/base factors; ambiguity -> error
    demographics.py         CKD-EPI 2021, Cockcroft-Gault, IBW/ABW/BSA/BMI, Child-Pugh
    pgx.py                  diplotype -> activity score -> phenotype -> CL multiplier
    ddi.py                  competitive / noncompetitive / mechanism-based / induction
    binding.py              fu, albumin dependence, free-drug normalization
    pk/
      model.py              compartment spec, dosing history -> input rate profile
      analytic.py           exact 1-cmt and 2-cmt (matrix exponential) superposition
      nonlinear.py          Michaelis-Menten RK4 with step-halving error estimate
      dialysis.py           intradialytic clearance, rebound, post-HD dosing timing
      scaling.py            allometry, renal/hepatic CL partition
    dosing/
      targets.py            target metric definitions and attainment evaluation
      design.py             candidate enumeration, linearity reuse, deterministic pick
      pta.py                Monte Carlo, fixed seed, CRN, Wilson CI
      safety.py             Cmax/AUC vs threshold, margin ratio, P(exceed)
    tdm/
      bayes.py              MAP, OFV, covariance, information ratio, convergence gate
      records.py            dosing history + level records, timing validation
    attribution.py          log-space exact covariate decomposition
    rationale.py            deterministic template text (no LLM)
    library/
      schema.json           JSON Schema (draft 2020-12) for entries
      loader.py             load, validate, hash, reject unresolved placeholders
      drugs/vancomycin.yaml tacrolimus.yaml apixaban.yaml phenytoin.yaml lithium.yaml
  tests/                    stdlib unittest; optional pytest via dev extra
    fixtures/patients/*.json
    test_demographics.py  test_pk_analytic.py  test_nonlinear.py
    test_binding.py  test_ddi.py  test_pgx.py  test_dialysis.py
    test_dosing.py  test_pta.py  test_bayes.py  test_attribution.py
    test_determinism.py  test_mass_balance.py  test_edge_cases.py
    test_library_integrity.py  test_validation_cases.py  test_no_network.py
  validation/
    cases/*.yaml            published anchor cases (input + expected + citation + tol)
    run_validation.py       table of pass/fail + max deviation; non-zero exit on fail
    GOLDEN.md               what each anchor proves
  demo/
    case_vanc.json  case_tacro.json  demo.sh  expected_*.json
  .github/workflows/ci.yml  python 3.11/3.12/3.13; unittest; validation; determinism

## 4. Input schemas

### 4.1 patient.json (schema_version "1.0")

Units are fixed by field name (`_kg`, `_cm`, `_mg_dL`, `_g_dL`). Only SCr and weight
accept an alternative unit block, converted explicitly; anything ambiguous is an error.

    {
      "schema_version": "1.0",
      "patient_id": "DEMO-VANC-001",
      "age_years": 68,
      "sex": "male",                       // required: CKD-EPI, CG, IBW
      "weight_kg": 92.0,
      "height_cm": 172.0,
      "serum_creatinine_mg_dL": 2.1,
      "creatinine_assay": "enzymatic",     // enzymatic | jaffe  (affects reliability note)
      "creatinine_stable": true,           // false -> W_SCr_NOT_STEADY_STATE
      "albumin_g_dL": 2.8,
      "hepatic": {
        "child_pugh_components": {
          "bilirubin_mg_dL": 2.4, "inr": 1.6, "albumin_g_dL": 2.8,
          "ascites": "moderate",           // none | mild | moderate_severe
          "encephalopathy": "grade_1"      // none | grade_1_2 | grade_3_4
        }
      },
      "renal": {
        "dialysis": null,                  // or {"modality":"hd","sessions_per_week":3,
                                           //     "hours_per_session":4.0,
                                           //     "blood_flow_mL_min":300,
                                           //     "dialysate_flow_mL_min":500,
                                           //     "dialyzer_koA_mL_min":900}
        "urine_output_mL_day": null
      },
      "genotypes": { "CYP2C9": "*1/*3", "CYP3A5": "*3/*3" },
      "phenotypes": {},                    // explicit wins; conflict with genotype -> error
      "co_medications": [
        { "name": "fluconazole", "dose_mg": 400, "route": "po", "frequency": "qd",
          "start_day": -5, "steady_state": true,
          "effects": { "CYP3A4": {"mechanism":"competitive","ki_umol_L":7.0},
                       "CYP2C9": {"mechanism":"competitive","ki_umol_L":20.0} } }
      ],
      "indication": "mrsa_bacteremia",
      "mic_mg_L": 1.0,
      "hematocrit_percent": null,          // tacrolimus binding covariate
      "pregnancy": false,
      "smoking_status": "non_smoker"
    }

Derived values never overwrite input: output has `patient.derived` with each value's
`formula` id and `source`.

### 4.2 Dosing history + levels (for TDM), levels.json

    {
      "schema_version": "1.0",
      "doses": [
        {"dose_mg": 2500, "start_h": 0.0,  "infusion_min": 120, "route": "iv"},
        {"dose_mg": 1250, "start_h": 12.0, "infusion_min": 90,  "route": "iv"}
      ],
      "levels": [
        {"time_h": 14.0, "value": 18.4, "unit": "mg/L", "matrix": "serum",
         "assay": "immunoassay", "note": "peak 2h after end of infusion"}
      ]
    }

Rules: every level must be attributable to a time origin (t=0 = first dose start);
levels before the first dose or with unknown timing are rejected; a level drawn during
an infusion is flagged (W_LEVEL_DURING_INFUSION) and still usable if the model can
place it.

### 4.3 Drug entry (YAML) — the curation contract

    schema_version: "1.0"
    drug: vancomycin
    display_name: Vancomycin
    drug_class: glycopeptide
    entry_version: "2026.11.1"
    status: resolved            # resolved | unresolved (unresolved -> engine refuses)

    # ---- structural model
    pk:
      model: two_compartment
      concentration_matrix: serum
      parameterization: CL_V_Q
      reference_individual: {weight_kg: 70, renal_reference: 100}   # units for typical
      typical:
        CL_L_h:   {value: ___, source: [S1]}
        V1_L:     {value: ___, source: [S1]}
        V2_L:     {value: ___, source: [S1]}
        Q_L_h:    {value: ___, source: [S1]}
      between_subject_variability:            # CV% + correlation matrix
        CL_L_h: {cv_percent: ___, source: [S1]}
        V1_L:   {cv_percent: ___, source: [S1]}
        correlation: {CL_L_h_V1_L: ___}
      residual_error: {type: proportional, cv_percent: ___, source: [S1]}

    # ---- covariate effects (each one cited; absent = no effect claimed)
    covariates:
      renal:
        fraction_renal_clearance: {value: ___, source: [S1]}
        function_of: crcl_cg_abw            # which renal metric this model used
        reference_value: 100
        form: proportional                  # proportional | power(exponent)
      weight:
        metric_for_cl: total                # total | ibw | abw
        metric_for_v: total
        cl_exponent: ___
        v_exponent: ___
      hepatic:
        child_pugh_multipliers: {A: ___, B: ___, C: ___}   # per-drug, cited, or omitted
        source: [S2]
      albumin:
        applies: false                      # true only if entry declares it matters
      genotype:
        CYP2C9: {fraction_metabolized: ___, source: [S3]}

    # ---- binding
    protein_binding:
      fraction_bound: ___
      albumin_dependence: none              # none | power(exponent, alb_ref)
      free_drug_is_active: false
      tdm_measures: total
      note: "TDM targets are total concentrations; do not free-correct."

    # ---- formulation / administration
    products:
      - name: "vial 500 mg"   base_mg: 500  salt_factor: 1.0
    salt_factors: {default: 1.0}
    routes:
      iv: {infusion_min_allowed: [60, 90, 120, 180], max_rate_mg_min: ___}
    absorption: null

    # ---- targets and safety (the only place target numbers may exist)
    targets:
      - id: auc24_mic
        metric: auc24_over_mic
        unit: h
        range: [400, 600]
        matrix: serum
        default: true
        source: [S4]
        note: "AUC24 by trapezoid at steady state; MIC 1.0 mg/L unless measured."
      - id: trough
        metric: trough_total
        unit: mg_L
        range: [10, 20]
        source: [S4]
        note: "Retained for trough-reporting labs; not the primary target."
    safety:
      - id: auc24_toxicity
        metric: auc24
        unit: mg_h_L
        threshold: 600
        direction: above
        source: [S4]
    monitoring:
      tdm_strategy: two_level_auc
      sampling: "peak 1-2 h after end of infusion + trough before next dose"

    # ---- decision pathway
    decision_basis: pk_target               # pk_target | label_criteria
    label_rules: []                         # used by apixaban-style entries

    # ---- dosing feasibility
    constraints:
      max_single_dose_mg: ___
      min_interval_h: 6
      max_duration_days: null

    citations:
      S1: {ref: "___", doi: "___", note: "popPK model used for CL/V/IIV"}
      S2: {ref: "___", doi: "___"}
      S3: {ref: "CPIC guideline ___", doi: "___"}
      S4: {ref: "Rybak MJ et al. Am J Health-Syst Pharm 2020;77:835-864", doi: "10.1093/ajhp/zxaa036"}
    provenance:
      curator: "___"
      curated_on: "___"
      reviewed_against: ["S1 table 2", "S4 recommendation 1"]
      known_gaps: ["___"]

Every numeric leaf carries `source`. `library validate` fails on: missing source, an
unresolved `___`, a target without a citation, a fraction_metabolized sum > 1, a
Child-Pugh multiplier outside (0,1.5), a negative IIV, or a correlation matrix that is
not positive definite. This is how the code makes curation first-class.

## 5. Computation pipeline

One deterministic sequence; each stage records what it used.

  1 validate inputs            -> typed values or a machine-readable error
  2 derive body size           -> IBW, ABW, BSA, BMI, weight metrics
  3 derive renal function      -> CKD-EPI 2021 eGFR and CG CrCl (multiple weights)
  4 derive hepatic status      -> Child-Pugh class (+ MELD 3.0 display only)
  5 load + hash library entry  -> resolved parameters or refusal
  6 structural parameters      -> typical CL, V, Q at the patient's covariates
  7 covariate multipliers      -> renal, weight, hepatic, PGx, DDI, binding
  8 regimen candidates         -> enumerated from library constraints
  9 exposure per candidate     -> exact linear solution or RK4 nonlinear
 10 PTA per candidate          -> Monte Carlo, common random numbers
 11 select regimen            -> max PTA subject to safety + practicality
 12 uncertainty               -> prediction interval, PTA CI, parameter CVs
 13 attribution               -> exact log-space decomposition of CL and V
 14 rationale + trace         -> deterministic text

### 5.1 Demographics (formula ids are printed with each value)

F_EGFR_CKDEPI2021 (Inker LA et al., NEJM 2021;385:1737-49), mL/min/1.73 m2:
    k = 0.7 (F) | 0.9 (M);  a = -0.241 (F) | -0.302 (M)
    eGFR = 142 * min(SCr/k,1)^a * max(SCr/k,1)^-1.200 * 0.9938^age * (1.012 if F)

F_CRCL_CG (Cockcroft DW, Gault MH. Nephron 1976;16:31-41):
    CrCl = (140 - age) * weight / (72 * SCr) * (0.85 if F)
    computed three times: weight = TBW, IBW, ABW. The entry declares which one its
    renal covariate was estimated against (`function_of: crcl_cg_abw`); the engine
    uses that one for CL and reports the others.

F_IBW_DEVINE: M 50 + 2.3*(inches - 60); F 45.5 + 2.3*(inches - 60)
F_ABW: IBW + 0.4*(TBW - IBW)      F_BSA_MOSTELLER: sqrt(H_cm * W_kg / 3600)
F_BMI: W / (H_m)^2

F_CHILD_PUGH: points 1-3 per component (bili, INR, albumin, ascites, enceph) ->
A 5-6, B 7-9, C 10-15.
F_MELD3: display only, coefficients transcribed from Kim WR et al. Gastroenterology
2021 and unit-tested against that paper's worked examples; if the coefficients cannot
be verified from source, MELD is omitted rather than approximated. MELD is never
used to scale clearance — it predicts mortality, not metabolic capacity.

Reference-weight choice is a clinical decision the engine surfaces, not hides:
obesity (BMI >= 30) always prints TBW, IBW and ABW renal functions and the % spread.

### 5.2 Structural PK

Internal units: mg, L, h, mg/L.

Linear exact solution (no numerical integration, no step-size error, mass conserving
by construction):
  1-cmt:  dA1/dt = R(t) - (CL/V1) A1
          A1(t+h) = A1 e^-k h + (R/k)(1 - e^-k h),  k = CL/V1
  2-cmt:  M = [[-(k10+k12), k21], [k12, -k21]], k10 = CL/V1, k12 = Q/V1, k21 = Q/V2
          x(t+h) = e^(Mh) x(t) + M^-1 (e^(Mh) - I) u,  u = [R, 0]
          e^(Mh) for 2x2 by Cayley-Hamilton: e^(Mh) = a0 I + a1 M, with
          a0 = (a e^(l1 h) - b e^(l2 h))/(a-b), a1 = (e^(l1 h)-e^(l2 h))/(a-b),
          l1,l2 = eigenvalues of M (= -alpha, -beta); repeated-eigenvalue branch
          handled explicitly. Any piecewise-constant input (multi-dose, interrupted
          infusion, dialysis gap) is exact and independent of the time grid.

Oral/extravascular: absorption compartment with first-order ka and optional lag;
F multiplies the dose. ER formulations require a cited input-rate model in the entry;
absent that, the ER route is refused rather than approximated.

Nonlinear (phenytoin): dA/dt = R(t) - (Vmax C)/(Km + C), C = A/Vd, RK4 with h = 1/60 h,
step-halving comparison reported as an integration error estimate; if the estimate
exceeds 1e-4 relative, h is halved until it does not (bounded at h = 1/600 h, then a
refusal rather than an unconverged number).

Dialysis: two mechanisms available per entry.
  (a) empirical: CL_HD cited in the entry, applied only during sessions, plus an
      interdialytic rebound interval;
  (b) mechanistic: CL_HD = Q_b (1 - Hct) (1 - exp(-KoA/Q_b_eff)) scaled by fu, with
      the dialysate-flow limitation min(Qb_eff, Qd) — used only where the entry
      supplies KoA and declares the model. Post-HD dosing timing is explicit in output
      ("dose 2500 mg after session 1"), intradialytic dosing is not silently assumed.

### 5.3 Covariate scaling (multiplicative on the log scale, which is what makes
attribution exact)

    CL = CL_ref
       * (RF / RF_ref)^theta_renal                 [RF = entry's declared renal metric]
       * (WT / 70)^0.75                            [allometric; per-drug metric/exponent]
       * CP_mult(Child-Pugh)                       [per-drug, cited, or 1.0 if absent]
       * PGx_mult                                  [1 + sum_g fm_g (AS_g/AS_ref - 1)]
       * DDI_mult                                  [product over interacting pathways]
       * fu_mult                                   [only if entry says binding matters]

    V  = V_ref * (WT/70)^1.0 * CP_mult_V          [V only where data support it]

Renal floor: CL_R is floored at the entry's declared floor (usually 0 for the renal
fraction, never negative); dialysis adds CL_HD on top, it does not multiply eGFR.

PGx: diplotype -> activity score from a versioned CPIC allele table per gene (allele
absent from the table = error, never a default of 1). CYP2C9 example: *1 = 1, *2 = 0.5,
*3 = 0, *5 = 0, *6 = 0, *8 = 0, *11 = 0, *13 = 0; AS 0 = PM, 0.5-1 = IM, 1.5-2 = NM,
>2 = RM/UM (gene-specific thresholds live in the table, cited to CPIC). Effect is
applied to the fraction of clearance that enzyme accounts for (fm), never to total CL
wholesale. For phenytoin/CYP2C9 the CPIC dose guidance (IM -25%, PM -50%) is an
independent validation anchor, not a substitute for the mechanistic path.

DDI: [I]u = unbound systemic concentration at steady state, from the co-medication's
declared dose and a cited volume/clearance, or from a declared Cmax,fu. Then

    competitive:      CL_int,app = CL_int / (1 + [I]u/Ki)
    noncompetitive:   same algebraic form, but flagged as concentration-independent
    mechanism-based: activity = kdeg / (kdeg + kinact [I]u/(KI + [I]u)), requires
                      kinact, KI, kdeg from the entry; time-to-onset is reported
    induction:        CL_int,app = CL_int * (1 + Emax [I]/(EC50+[I])), onset t1/2 from
                      the entry; v1 applies the steady-state factor and prints
                      "steady state assumed; full induction in ~X days"

Applied through fm: dCL = CL * fm_pathway * (multiplier - 1). A DDI on a pathway the
drug does not use changes nothing, and the attribution says so explicitly.

Binding: fu = fu_ref * (alb_ref/alb)^n where the entry declares n; CL scales with fu
only for restrictively cleared drugs (entry declares `restrictively_cleared`). For
high-extraction drugs total CL is unchanged while measured total concentration shifts,
which is a target-interpretation issue, not a CL issue. Phenytoin additionally gets the
Winter-Tozer normalization (C_norm = C/(0.2 Alb + 0.1); renal failure 0.1 Alb + 0.1)
as an explicitly labeled display transform, never as a silent correction.

### 5.4 Regimen candidates

Enumerated deterministically from library constraints, not from a free grid:
- loading dose from the entry's loading rule (e.g. 25-30 mg/kg of a declared weight
  metric) cross available strengths, capped by max_single_dose,
- maintenance dose levels across available strengths within an entry-declared
  mg/kg band, intervals from the entry's allowed list,
- route/formulation combinations the entry permits.

Candidates that are not administrable (dose not divisible by product strength within
tolerance, infusion rate above max) are dropped, and the drop list is reported in
`--trace` so a pharmacist sees what was excluded and why.

### 5.5 PTA and regimen selection

Per candidate: simulate N patients (default 5000).

    CL_i = CL_pop * exp(eta_CL_i),  V_i = V_pop * exp(eta_V_i)
    [eta_CL, eta_V] ~ N(0, Omega), Omega from IIV% + correlation, Cholesky once.

Sampling uses a documented, version-stable path: MT19937 (stdlib `random.Random`)
-> uniform -> inverse normal CDF by Wichura AS241 implemented in-repo (~1e-15
accuracy). This avoids `random.gauss`'s implementation dependence and makes output
identical across Python versions, not just within one.

Steady state is computed by iterating the exact map to convergence (relative change
in A1 and A2 < 1e-9, max 60 doses) rather than by a closed-form shortcut; the number
of doses to steady state and the residual are reported.

Common random numbers: the same eta draws are reused across all candidates and the
same seed is used for every run, so candidate rankings are not Monte Carlo noise.
PTA = fraction with the target metric inside its range at steady state, with a Wilson
95% interval. PST (probability of successful therapy) = P(attained AND no safety
breach) is also reported, because those are different clinical questions.

Selection: maximize PTA subject to (a) P(exceed safety threshold) <= 0.05,
(b) dose not exceeding constraints, (c) interval >= min_interval. Tie-break, in order:
higher PST, lower P(exceed), lower total daily dose, longer interval, fewer
administrations/day, then lexicographic on (dose, interval). The full candidate table
with PTA/PST/P(exceed) is emitted, so the choice is auditable rather than asserted.

If no candidate reaches PTA >= 0.80, the tool still returns the best candidate with
`status: target_unattainable` and the reason (e.g. MIC too high for max dose), exit 4.
It never silently widens the target range.

### 5.6 Bayesian TDM (MAP)

Objective (log-normal residual, library priors):

    OFV(eta) = sum_j (ln Cobs_j - ln Cpred_j(eta))^2 / sigma_j^2 + eta^T Omega^-1 eta
    sigma_j^2 = (cv_prop^2) + (sd_add/C_j)^2        [entry declares the error model]

Minimized by Levenberg-damped Gauss-Newton with a finite-difference Jacobian
(central differences, step 1e-5 in log-parameter space). Deterministic start:
eta = 0 (no random restarts). Convergence gate: gradient infinity-norm < 1e-6 within
60 iterations, else `status: not_converged` and no regimen is recommended.

Reported: EBE (eta), posterior CL/V with CV% from (J^T W J + Omega^-1)^-1, OFV,
weighted residuals per level, prior->posterior information ratio per parameter
(1 - posterior_var/prior_var; < 0.2 -> W_WEAK_INFORMATION, "this answer is mostly the
prior"), and the predicted-vs-observed check per level.

Re-dosing after TDM uses the posterior distribution, not the point estimate: the PTA
for the new regimen integrates over posterior parameter uncertainty as well as IIV,
and the report states which of the two dominates. This is the difference between a
number and an honest number.

### 5.7 Attribution

Because every covariate effect is multiplicative on CL and on V, log CL is additive:

    ln CL = ln CL_ref + sum_k contribution_k

Toggle each effect from off to on one at a time, holding the others fixed, and report
the exact percentage contributions plus a residual check that they sum to the total
log-change (asserted < 1e-12 relative). This is not a heuristic ranking; it is the
model's own algebra. Interaction terms are reported separately where they exist
(e.g. the renal x weight interaction through ABW), never folded silently into a main
effect.

Example output line:
  "Maintenance CL 4.9 -> 2.6 L/h (-47%): eGFR 28 vs 100 (-34%), fluconazole CYP3A4
   competitive inhibition at [I]u 6.1 uM vs Ki 7.0 uM (-9%), Child-Pugh B (-4%)."

### 5.8 Rationale and trace

`rationale` is a deterministic template over the structured result (no LLM): what was
asked, which renal metric was used and why, what drove CL, what the regimen is, what
the exposure and PTA are, what the margin is, what is uncertain, what was missing.
`--trace` additionally prints every formula with its substituted numbers
(auditable "show your work" page). Same input -> same bytes.

## 6. Output contract (regimen.json)

    {
      "schema_version": "1.0",
      "status": "ok",                     // ok | target_unattainable | insufficient_data
      "engine": {"name":"precision-dose","version":"0.1.0","seed":20260101,
                 "rng":"mt19937+AS241","model_hash":"sha256:...","library_hash":"sha256:..."},
      "question": {"drug":"vancomycin","target_id":"auc24_mic","indication":"...",
                   "mic_mg_L":1.0},
      "patient": {"input": {...}, "derived": {"egfr_ckd_epi_2021": {"value":...,"formula":"F_EGFR_CKDEPI2021"}, ...}},
      "library_entry": {"drug":"vancomycin","entry_version":"...","content_hash":"sha256:...",
                        "citations":[...], "known_gaps":[...]},
      "parameters": {"CL_L_h":..., "V1_L":..., "V2_L":..., "Q_L_h":...,
                     "individualized":false, "posterior":null},
      "attribution": [{"covariate":"renal","delta_cl_pct":-34.1,"basis":"eGFR 28 vs 100, f_renal 0.85"}],
      "regimen": {"loading":{"dose_mg":2500,"infusion_min":120,"at_h":0},
                  "maintenance":{"dose_mg":1250,"interval_h":12,"infusion_min":90},
                  "decision_basis":"pk_target","candidates_considered":N},
      "exposure": {"auc24_mg_h_L":..., "cmax_mg_L":..., "trough_mg_L":...,
                   "t_half_h":..., "accumulation_ratio":..., "tss_doses":...},
      "pta": {"n":5000,"seed":20260101,"probability":0.87,"ci95":[0.861,0.879],
              "criterion":"auc24/mic 400-600","pst":0.84},
      "safety": {"metric":"auc24","threshold":600,"predicted":...,"margin_ratio":...,
                 "margin_sd":..., "p_exceed":0.02},
      "uncertainty": {"auc24_pi90":[...], "cl_cv_percent":..., "pta_ci95":[...]},
      "warnings": [{"code":"W_CG_EGFR_DISCORDANCE","severity":"info","message":"..."}],
      "rationale": "...",
      "curve": {"time_h":[...], "conc_mg_L":[...], "matrix":"serum"},
      "disclaimer": "Computational decision support. Requires pharmacist review..."
    }

Stable formatting: keys sorted, floats serialized by Python's shortest-roundtrip repr,
values rounded to 6 significant digits before serialization, no timestamps unless
`--timestamp` is passed, no absolute paths, no environment-dependent values.

## 7. CLI and HTTP contract

    precision-dose dose    --patient P.json --drug vancomycin
                           [--target auc24_mic] [--renal crcl_cg_abw|egfr_ckd_epi_2021]
                           [--weight-metric total|ibw|abw] [--levels levels.json]
                           [--seed 20260101] [--pta-n 5000] [--curve-points 200]
                           [--out regimen.json] [--format json|text] [--trace]
    precision-dose tdm     --patient P.json --drug vancomycin --levels levels.json
                           [--redose] [--out tdm.json]
    precision-dose library list|show <drug>|validate [--library DIR]
    precision-dose validate [--cases validation/cases] [--json]
    precision-dose serve   [--host 127.0.0.1] [--port 8765]

Exit codes: 0 ok; 2 invalid input (schema/units/unknown drug); 3 insufficient data or
refusal (unresolved entry, non-convergence, pediatric); 4 target unattainable (best
candidate returned, not hidden).

HTTP (stdlib ThreadingHTTPServer, localhost by default, stateless, no auth, no storage):
  POST /v1/dose    body = {patient, drug, target?, options?} -> regimen.json shape
  POST /v1/tdm     body = {patient, drug, dosing_history, levels} -> MAP result
  GET  /v1/library            -> entry index
  GET  /v1/library/{drug}     -> entry + hash
  GET  /healthz               -> {"status":"ok","version":"..."}
Error bodies: {"error":{"code":"E_MISSING_FIELD","field":"weight_kg","message":"..."}}
422 on validation failure, 404 unknown drug, 400 malformed JSON. Same engine, same
determinism: the HTTP layer is a transport, nothing more.

## 8. Determinism and auditability

- No wall-clock, no PID, no hostname, no locale-dependent formatting in output.
- One seed drives everything; CRN across candidates; AS241 inverse normal.
- Sorted keys, quantized floats, stable list orders (curated, not hash-order).
- Proof, not promise: `test_determinism.py` hashes output across repeated runs, across
  PYTHONHASHSEED=0/1/2, across input key order permutations, and across process
  restarts; all four hashes must be equal.
- `test_no_network.py` patches socket.socket to raise and runs a full dose computation;
  any network attempt fails the suite.
- Every output carries library_hash + model_hash so a regimen can be reproduced from
  the exact entry version that produced it.

## 9. Refusal and warning policy

Never impute silently. Three tiers:
- REFUSE (exit 3): missing required field; unresolved library entry; allele not in the
  CPIC table; optimizer not converged; age < 18; ER route without a cited input model;
  dialysis declared but modality parameters absent when the entry needs them.
- DEGRADE with warning: height missing (IBW/ABW paths unavailable, eGFR still valid);
  Child-Pugh incomplete (hepatic multiplier not applied, stated); creatinine not at
  steady state; albumin missing (binding correction not applied, stated).
- INFORM: CG vs CKD-EPI discordance > 30%; obesity weight ambiguity; MIC assumed 1.0;
  induction not at steady state; PTA CI width > 0.1.

Warning codes are stable strings so downstream CDS can render or suppress them.

## 10. Validation program (validation/cases/*.yaml, run by `precision-dose validate`)

Each case: id, citation, inputs, expected quantity, tolerance, and a one-line note on
what it proves. Passing is a gate in CI, not a report.

Anchors planned:
V1 ASHP/IDSA/PIDS 2020 vancomycin (Rybak 2020, doi:10.1093/ajhp/zxaa036):
   AUC24/MIC 400-600 for MRSA; loading 25-30 mg/kg actual body weight; the guideline's
   example regimens must land inside 400-600 with our model and default target.
V2 Vancomycin framework-vs-trough comparison: same regimen judged by AUC24 vs trough
   10-20 must reproduce the guideline's stated discordance, demonstrating why AUC is
   the primary target (the tool reports both).
V3 CKD-EPI 2021 worked examples (NEJM 2021) to 3 decimal places.
V4 Cockcroft-Gault against published worked examples, three weight metrics.
V5 Apixaban label criteria: >=2 of (age >=80, weight <=60 kg, SCr >=1.5) -> 2.5 mg BID;
   strong dual CYP3A4+P-gp inhibitor -> 2.5 mg BID or avoid; Child-Pugh C -> not
   recommended. Decision-pathway test: dose comes from label_rules, not from CL.
V6 Phenytoin Michaelis-Menten: published Vmax/Km reproduce a documented
   disproportionate-increase case; Winter-Tozer normalization against a published
   example.
V7 Lithium: published clearance-vs-CrCl relation and a published maintenance regimen
   reproduce 0.6-0.8 mmol/L at steady state.
V8 Tacrolimus: trough target by phase; CYP3A5*3/*3 vs *1/*1 dose requirement difference
   in the direction and magnitude of the cited model; azole DDI magnitude.
V9 KDIGO renal-dosing consistency: for each drug, the label/guideline renal cutoffs are
   encoded as tests, not prose.

Where a cited source does not supply a number the entry needs, the entry ships
`status: unresolved` and the drug is refused with a message naming the missing
parameter. No anchor is invented to make a test green; a missing anchor is reported as
missing in GOLDEN.md.

## 11. Test plan (beyond the anchors)

- Unit: demographics, unit conversion, salt factors, PGx tables, DDI algebra, binding.
- Analytic correctness: 1-cmt/2-cmt against closed-form superposition; 2-cmt infusion
  against the classic biexponential infusion equation; matrix exponential vs
  eigen-decomposition for random parameter sets.
- Mass balance (property-based, deterministic random parameter draws): amount in +
  amount eliminated = dose, to 1e-12 relative; concentration never negative.
- Monotonicity: AUC non-decreasing in dose; CL_R non-decreasing in renal function;
  CL non-increasing as Ki decreases; V non-decreasing in weight. Violations are bugs.
- Nonlinear: RK4 step-halving convergence; phenytoin at Km (half-saturation) matches
  the analytic MM steady state Vmax*C/(Km+C) = R.
- Bayes: simulate with known CL/V + noise, refit, recover CL within 15% with two
  levels; single-level case must report weak information rather than false precision;
  prior-only (no levels) must return the prior exactly.
- Edge cases: HD patient; Child-Pugh C with no data (refusal); BMI 45 (TBW/IBW/ABW
  spread printed); age 19 and 95; anuria; albumin 1.5; height missing; SCr 0.4
  (low muscle mass); MIC 2.0 with max dose (target_unattainable path).
- Real profiles: the five demo patients in tests/fixtures are complete, clinically
  coherent cases, each asserted against hand-computed exposure.
- Library integrity: schema, citations, hashes, unresolved placeholders, positive-
  definite correlation matrices, fm sums.
- Performance: `dose` with n=5000 under 3 s; `validate` under 10 s.

## 12. Acceptance criteria (v1 is done when all hold, verified by execution)

A1 `precision-dose dose` produces regimen + curve + PTA + safety margin + rationale for
   all five drugs (or a cited refusal for any unresolved entry).
A2 Byte-identical output across repeated runs, PYTHONHASHSEED values, input key order,
   and process restarts (four hashes equal).
A3 `precision-dose validate` passes every anchor V1-V9 that has a verified source, and
   exits non-zero if any verified anchor deviates beyond tolerance.
A4 Mass balance and monotonicity properties hold for every drug entry.
A5 Bayesian TDM recovers known parameters within stated tolerance and refuses to report
   an unconverged fit; weak-information cases are labeled.
A6 All five edge-case classes (dialysis, severe hepatic, obesity, extreme age, missing
   data) behave per the refusal/warning policy, with tests.
A7 README worked case is generated by the tool and pasted verbatim, and a golden test
   asserts the README block matches current output (docs cannot drift from behavior).
A8 `serve` answers /v1/dose with output identical to the CLI for the same input.
A9 Zero runtime dependencies: clean venv, `pip install .`, run all commands, no
   third-party import in the engine (asserted by an import test).
A10 No network: suite passes with socket disabled.
A11 MIT license, CHANGELOG, CURATION.md checklist, and library entries whose every
   number carries a source.

## 13. Curation workflow (CURATION.md)

For each drug: pick one primary popPK publication (not a review, not a dosing app);
transcribe parameters and IIV with table/section references; declare which renal
metric the model used; declare fm per pathway; declare targets only from the governing
guideline, with the exact recommendation number; declare safety thresholds with source;
add at least one validation anchor; run `precision-dose library validate`; record
known gaps honestly. Entry status stays `unresolved` until every `___` is filled and
every anchor passes. The curator's name and review date are part of the entry, and the
content hash makes any later edit visible in every output that used the old data.

Candidate sources to transcribe (you validate; I will fetch and quote the primary text
during implementation, and leave the entry unresolved if the source does not support a
number):
   vancomycin  - Rybak 2020 guideline (targets/loading) + one published popPK model
   tacrolimus  - a popPK model including CYP3A5 genotype and hematocrit
                 (e.g. Passey et al., Br J Clin Pharmacol 2011) + phase-based troughs
   apixaban    - US prescribing information criteria + one popPK model for exposure
                 reporting only (decision_basis: label_criteria)
   phenytoin   - classic Michaelis-Menten parameters + Winter-Tozer correction
   lithium     - published CL-vs-CrCl relation + maintenance target range

## 14. Risks and mitigations

R1 Hand-rolled numerics could be subtly wrong.
   Mitigation: exact linear solutions where possible, scipy cross-checks in a dev-only
   test extra, analytic identities (mass balance, MM half-saturation, biexponential
   infusion equation) as independent oracles.
R2 PopPK models disagree; choosing one is a clinical judgment.
   Mitigation: one model per entry, named and cited, disagreement surfaced in
   known_gaps; the schema supports multiple candidate models later without changing
   the engine.
R3 Monte Carlo PTA can look precise while being wrong.
   Mitigation: Wilson CI, MC standard error, CRN across candidates, and a reported
   convergence check on steady-state iteration.
R4 Child-Pugh scaling is coarse and often unsourced.
   Mitigation: per-drug multipliers only, absence is a refusal to apply, never a
   generic 0.75/0.5 table.
R5 Unit and matrix errors (whole blood vs serum; salt vs base).
   Mitigation: matrix is a required field on every target and on the output curve;
   salt_factor required on every product; tests assert both.
R6 Scope creep into a platform.
   Mitigation: no batch, no storage, no auth, one stateless request per process.

## 15. Implementation phases (after approval)

P1 skeleton: pyproject, MIT, schemas, units, io, demographics + tests (CKD-EPI, CG,
   IBW/ABW, Child-Pugh) — the parts that are pure arithmetic and independently
   checkable.
P2 PK core: 1-cmt/2-cmt exact solutions, nonlinear RK4, dialysis, scaling; mass-balance
   and monotonicity property tests.
P3 dosing: targets, candidates, PTA with CRN, selection, safety; attribution.
P4 library: schema, loader, validator, hashing, the five entries with citations and
   anchors; validation harness.
P5 TDM: MAP, covariance, information ratio, re-dosing; recovery tests.
P6 interface: CLI, text report, trace, serve, README (generated worked case), demo.
P7 hardening: determinism across environments, no-network test, performance, CI,
   docs, CHANGELOG.

Estimated v1 shape: ~4.5-6k lines including tests and library YAML; no file over ~400
lines except the library schema.

## 16. Open questions for the curator (non-blocking, defaults stated)

Q1 Vancomycin renal metric: the 2020 guideline's initial dosing uses CrCl (CG) while
   AUC targeting uses Bayesian estimation. Default: `function_of: crcl_cg_abw` for the
   CL covariate, with CKD-EPI reported alongside. Confirm.
Q2 Vancomycin target default: AUC24/MIC 400-600 with MIC 1.0 assumed, trough retained
   as a selectable secondary target. Confirm.
Q3 Tacrolimus formulation: IR only in v1 (ER needs a cited input-rate model). Confirm.
Q4 Lithium: accept the published CL = k * CrCl relation as the mechanistic core, with
   k from the primary source, plus Vd from the same source. Confirm.
Q5 Whether apixaban should be included at all in v1, given its dose is label-driven
   rather than PK-target-driven; my view is yes, because it demonstrates the
   `decision_basis: label_criteria` pathway and the honesty about which mechanism drove
   the answer. Confirm or drop to v1.1.

---

## As built (v0.1.0)

Delivered as designed, with the deviations below. Everything else — the pipeline stages,
the refusal policy, the PTA method, the MAP estimator, the attribution algebra, the CLI
and HTTP contracts — is implemented as specified.

**Deviations**

1. **Library format is TOML, not YAML.** Parsed with the standard library's `tomllib`.
   YAML would have required either a runtime dependency (breaking the zero-dependency
   promise) or a hand-rolled parser, and a hand-rolled parser is exactly the kind of
   unaudited component this tool exists to avoid. The schema and the validation rules are
   unchanged.
2. **MELD is not computed at all.** The design allowed a display-only MELD 3.0 with
   coefficients transcribed from the primary paper. Those coefficients could not be
   verified from a primary source in this build, so the score is omitted rather than
   approximated. Child-Pugh is the only hepatic covariate, and it never applies without a
   per-drug multiplier.
3. **Dialysis is an entry-declared intradialytic clearance**, not a KoA-based mechanistic
   model. Vancomycin carries 6.0 L/h (100 mL/min) from the recovery studies; the
   post-session rebound is predicted by the two-compartment model rather than assumed. A
   renally cleared drug with no declared dialysis clearance (lithium) is refused for a
   dialysis patient.
4. **`target_unattainable` is reserved for a genuinely unreachable target.** A fixed dose
   that attains the window in only 43% of patients is reported as `ok` with a W_LOW_PTA
   warning, because that is the quantitative case for TDM, not a failure of the target.
   Only a maximum-dose median below the target's lower bound produces
   `target_unattainable` (exit 4).
5. **Apixaban produces no PK output.** With no sourced population model for the exposure
   calculation, the entry declares `structure = "none"` and the engine returns criteria,
   the fired rule and the interaction class evaluation, with `exposure: null`. Inventing a
   curve would have been the easier and worse choice.
6. **Tacrolimus V/F, IIV, ka and lithium's absorption parameters are labelled
   assumptions** (`citations.S4`/`S3` = "modelling assumption") rather than presented as
   estimates. They are bounded by tests that show the trough is insensitive to ka, and
   they appear in `known_gaps` in every rationale.

**Verification performed**

- 81 unit/property tests pass (analytic PK oracles, mass balance, monotonicity, Bayesian
  recovery, determinism across four conditions, no-network, edge cases, library
  integrity, documentation drift).
- 8/8 validation anchors pass; `precision-dose validate` exits non-zero on drift.
- Wheel installs into a clean venv with no third-party packages and runs.
- HTTP output is byte-identical to the CLI for the same input.
- The README worked case is re-run by a test; its headline numbers cannot drift.
