# Changelog

## 0.1.1 — 2026-11-13

Library revision: the four numbers flagged at 0.1.0 as thin are replaced with
published values, and two engine capabilities were added to carry them faithfully.

### Library
- **Tacrolimus** rebuilt on a single model (Andrews 2019, Br J Clin Pharmacol
  85:601-615): two-compartment, CL/F 23.0 L/h with 38.6% IIV, V1/F 692 L, V2/F 5340 L,
  Q/F 11.6 L/h, ka 3.58/h, all from the same cohort. The previous entry mixed a CL/F
  equation from one model with a V/F from another, and declared IIV as an assumption;
  both are gone. CYP3A5 (1.63x for any *1 carrier) and CYP3A4*22 (0.80x) are applied as
  diplotype multipliers. The model's own V2 is large and poorly identified, so the
  engine now reports the accumulation-based effective half-life over the dosing
  interval rather than a terminal phase the data cannot support.
- **Vancomycin** now carries the non-renal clearance intercept S1 omitted:
  +3.66 mL/min (Matzke 1984). It is applied as an additive component and reported as
  one, because folding an intercept into a multiplier would invalidate the log-space
  attribution.
- **Lithium** clearance variability is 24% from Jermain 1991 rather than an invented
  22%; between-subject variability in volume is no longer invented at all (it is set to
  zero and the resulting over-confidence is stated).
- **Fluconazole → CYP3A4** uses the measured in vivo clearance ratio (0.38, Ikeda 2006,
  tacrolimus CL 38.8 → 14.7 L/h) instead of an in-vitro Ki extrapolation. The Ki route
  is retained in the note as the cross-check: it predicts 0.18, which is stronger than
  is observed, because systemic unbound concentration overstates the hepatic
  concentration available to the enzyme.

### Engine
- Continuous covariates from published equations: body surface area, albumin,
  haematocrit and serum creatinine, each as a ratio to the model's own reference, each
  refusing rather than defaulting when the value is absent.
- Additive non-renal clearance support, with the attribution reporting the scaled and
  additive components separately and the log residual computed against the scaled part.
- Entries can declare `requires` fields; the engine refuses early, naming the field.
- Effective half-life for multi-compartment models (accumulation-based).
- `empirical_fold` DDIs apply to total clearance directly instead of through fraction
  metabolised, where they would double count.
- Tie-break between regimens with identical PTA and daily dose now prefers the interval
  closest to the guideline's usual one, declared per entry as `preferred_interval_h`.

## 0.1.0 — 2026-11-13

First working version.

### Engine
- One- and two-compartment linear PK solved exactly (matrix exponential; periodic
  regimens at the cycle fixed point), Michaelis-Menten elimination for phenytoin
  (RK4 with the steady state evaluated algebraically).
- Individualisation from CKD-EPI 2021 eGFR and Cockcroft-Gault CrCl (three weight
  metrics), Child-Pugh, allometric weight scaling, CPIC genotype multipliers,
  competitive/mechanism-based/induction DDIs applied through fraction metabolised,
  albumin binding where the entry declares it matters.
- Monte Carlo PTA with common random numbers, Wilson intervals, AS241 inverse-normal
  sampling so results do not depend on the Python build.
- Bayesian MAP estimation with a convergence gate, posterior covariance, and a
  prior-to-posterior information ratio; re-dosing integrates posterior uncertainty.
- Exact log-space covariate attribution, asserted to close to zero residual.
- Dialysis: intradialytic clearance from the entry, post-session rebound predicted by
  the two-compartment model; refusal when a renally cleared drug has no dialysis data.
- Refusal policy: no silent imputation. Unresolved entries, unknown alleles, missing
  height/albumin/Child-Pugh and paediatric ages all have defined, coded behaviour.

### Library
- Five curated entries: vancomycin, tacrolimus, apixaban, phenytoin, lithium.
- `decision_basis` distinguishes PK-target drugs from label-criteria drugs (apixaban).
- Perpetrator interaction table: fluconazole, ketoconazole, rifampicin, omeprazole.
- Every entry content-hashed; the hash is printed in every result.

### Interface
- CLI: `dose`, `tdm`, `library`, `validate`, `serve`; `--trace` prints the derivation.
- Local HTTP API with the same engine and byte-identical output.
- Validation harness with eight anchors, non-zero exit on drift.

### Notes
- The library format is TOML parsed with the standard library (`tomllib`), not YAML as
  originally proposed: YAML would have required a third-party parser or a hand-rolled
  one, and a hand-rolled parser is an unaudited component in a tool whose value is
  auditability.
- Zero runtime dependencies.
