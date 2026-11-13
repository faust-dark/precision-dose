# Curation: how a drug enters the library

The library is the part of this project a clinical pharmacist owns. The code exists to
make that ownership mechanical, so that an error in a number cannot hide behind a
functioning pipeline.

## The contract

A drug entry is one TOML file in `src/precision_dose/library/drugs/`. It is valid only
if:

- every parameter used in a calculation names a citation key that exists in
  `[citations.*]`;
- every target and safety threshold names its source;
- `provenance.curator` and `provenance.curated_on` are filled in;
- no value is still the string `PLACEHOLDER`;
- variances are non-negative, fractions are in [0, 1], and every `ki_umol_L` is positive;
- `status = "resolved"` — otherwise the engine refuses to dose the drug and says why.

`precision-dose library validate` enforces all of it, and the test suite runs the same
check, so a bad entry cannot reach a patient through a green pipeline.

## Step by step

1. **Choose one primary population PK publication.** One model, not a review, not a
   dosing app, not a guideline table. Record the model's reference individual (weight,
   renal function, age) and its covariate structure. If the paper reports CL and V but not
   IIV, the IIV is an assumption and must be labelled as one.
2. **Transcribe the fixed effects** into `[model]`: `cl_ref`, `v1_ref`, and for a
   two-compartment model `v2_ref` and `q_ref`. These are the values for the *reference
   individual*, not for an average patient in your ward.
3. **Transcribe the covariate structure** into `[model.scaling]`. Declare which renal
   metric the model was estimated against (`crcl_cg_total`, `crcl_cg_ibw`,
   `crcl_cg_abw`, `egfr_ckd_epi_2021`). Declare exponents; do not "improve" them.
4. **Declare the fractions metabolised** in `[model.ddi]`. Without them the DDI engine
   cannot know which interactions matter, and an inhibitor of an unused pathway will be
   reported as inapplicable rather than silently applied.
5. **Targets and safety** come from the governing guideline, with the recommendation
   number in the note. Never invent a range; if the guideline gives two targets (AUC and
   trough), encode both and mark one default.
6. **Constraints**: strengths, permitted intervals, maximum single and daily doses,
   loading rule if the guideline has one. These bound the search; without them the engine
   will happily propose a 12 g dose.
7. **Validation anchors**: add at least one case to `validation/cases/` that reproduces a
   published regimen or a published parameter value, with a tolerance and a note saying
   what it proves.
8. **known_gaps**: write the limitations a colleague would want to know — which parameter
   is borrowed from another cohort, what is assumed, what the model does not represent.
   This text is printed in every rationale, so it is read by the person at the bedside,
   not just by a reviewer.

## What "unresolved" means

If any number cannot be sourced, leave `status = "unresolved"`. The engine refuses the
drug with `E_LIBRARY_UNRESOLVED` and lists the gaps. A refusal is a defensible clinical
answer; a plausible number that nobody can trace is not.

## Editing an existing entry

Entries are content-hashed (sha256 over the canonical form) and the hash is printed in
every output. Changing a number changes the hash, which means every regimen produced
afterwards is distinguishable from those produced before. Bump `entry_version` and record
what changed and why in `CHANGELOG.md`.

## Adding an interaction

Perpetrators live in `library/interactions.toml`: mechanism, Ki (or kinact/kdeg, or Emax),
the unbound systemic concentration model, and a confidence label. A perpetrator without a
source is rejected by the same validator.

## Review checklist before merging a library change

- [ ] One named primary source per parameter, with a table or equation pointer.
- [ ] Reference individual stated.
- [ ] Renal metric declared and consistent with the source's own covariate.
- [ ] Targets traceable to a guideline recommendation number.
- [ ] At least one validation anchor added or updated.
- [ ] `known_gaps` updated honestly.
- [ ] `precision-dose validate` passes.
- [ ] `python -m unittest discover -s tests -t .` passes.
