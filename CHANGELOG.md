# Changelog

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
