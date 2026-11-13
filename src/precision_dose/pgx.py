"""CYP genotype -> activity score -> phenotype -> clearance multiplier.

Consistent with CPIC activity-score tables. An allele that is not in the table
is an error, never a default of 1.0: guessing an activity score is how a
pharmacogenomic tool silently becomes wrong.
"""

from .errors import E_BAD_VALUE, EngineError

# Activity scores by allele, per gene. Sources are the CPIC activity score
# tables (CYP2C9: CPIC 2020 update; CYP2C19: CPIC 2022; CYP3A5: CPIC 2015;
# CYP2D6: CPIC 2021). Curation note: only alleles with a well-established
# consensus score are listed; anything else refuses rather than defaults.
ALLELE_SCORES = {
    "CYP2C9": {
        "*1": 1.0, "*2": 0.5, "*3": 0.0, "*5": 0.0, "*6": 0.0,
        "*8": 0.0, "*11": 0.0, "*13": 0.0,
    },
    "CYP2C19": {
        "*1": 1.0, "*2": 0.0, "*3": 0.0, "*4": 0.0, "*5": 0.0, "*6": 0.0,
        "*7": 0.0, "*8": 0.0, "*9": 0.5, "*10": 0.5, "*17": 1.5,
    },
    "CYP2D6": {
        "*1": 1.0, "*2": 1.0, "*4": 0.0, "*5": 0.0, "*9": 0.5, "*10": 0.25,
        "*17": 0.5, "*41": 0.5,
    },
    "CYP3A5": {"*1": 1.0, "*3": 0.0, "*6": 0.0, "*7": 0.0},
}

# Phenotype thresholds on the activity score, per gene, as (min_as, phenotype)
# ordered descending. CPIC calls these normal/intermediate/poor/rapid/ultrarapid
# metaboliser.
PHENOTYPE_THRESHOLDS = {
    "CYP2C9": [(2.0, "normal_metabolizer"), (1.5, "normal_metabolizer"),
               (0.5, "intermediate_metabolizer"), (0.0, "poor_metabolizer")],
    "CYP2C19": [(3.0, "ultrarapid_metabolizer"), (2.5, "rapid_metabolizer"),
                (1.5, "normal_metabolizer"), (1.0, "normal_metabolizer"),
                (0.5, "intermediate_metabolizer"), (0.0, "poor_metabolizer")],
    "CYP2D6": [(2.5, "ultrarapid_metabolizer"), (1.5, "normal_metabolizer"),
               (1.0, "normal_metabolizer"), (0.5, "intermediate_metabolizer"),
               (0.0, "poor_metabolizer")],
    "CYP3A5": [(2.0, "normal_metabolizer"), (1.0, "intermediate_metabolizer"),
               (0.0, "poor_metabolizer")],
}

# Reference activity for "normal" used to scale intrinsic clearance.
REFERENCE_SCORE = {"CYP2C9": 2.0, "CYP2C19": 2.0, "CYP2D6": 2.0, "CYP3A5": 2.0}


def diplotype_score(gene, diplotype):
    table = ALLELE_SCORES.get(gene)
    if table is None:
        raise EngineError(
            E_BAD_VALUE,
            f"no activity-score table for {gene}; supported: {sorted(ALLELE_SCORES)}",
            "genotypes",
        )
    parts = [p.strip() for p in diplotype.split("/")]
    if len(parts) != 2:
        raise EngineError(
            E_BAD_VALUE, f"{gene} diplotype {diplotype!r} must have two alleles", "genotypes"
        )
    total = 0.0
    for allele in parts:
        if allele not in table:
            raise EngineError(
                E_BAD_VALUE,
                f"allele {allele} for {gene} has no CPIC activity score in this "
                "library; refusing to guess (add it to pgx.ALLELE_SCORES with a "
                "citation first)",
                "genotypes",
            )
        total += table[allele]
    return total


def phenotype_for(gene, score):
    for threshold, name in PHENOTYPE_THRESHOLDS.get(gene, []):
        if score >= threshold:
            return name
    return "indeterminate"


def resolve(patient, gene):
    """Return {gene, source, activity_score, phenotype, multiplier} or None.

    Explicit phenotype wins over genotype; if both are present and disagree the
    call is refused (a silent override is a clinical hazard).
    """
    genotypes = patient.get("genotypes") or {}
    phenotypes = patient.get("phenotypes") or {}
    explicit = phenotypes.get(gene)
    if gene in genotypes:
        score = diplotype_score(gene, genotypes[gene])
        derived_ph = phenotype_for(gene, score)
        if explicit is not None and explicit != derived_ph:
            raise EngineError(
                E_BAD_VALUE,
                f"{gene}: supplied phenotype {explicit!r} contradicts genotype "
                f"{genotypes[gene]!r} (activity score {score} -> {derived_ph})",
                "phenotypes",
            )
        ref = REFERENCE_SCORE[gene]
        return {
            "gene": gene,
            "diplotype": genotypes[gene],
            "activity_score": score,
            "phenotype": derived_ph,
            "relative_activity": score / ref,
            "source": f"CPIC activity score table ({gene})",
            "basis": "genotype",
        }
    if explicit is not None:
        # Phenotype without genotype: map phenotype to a canonical activity
        # fraction. Documented, not inferred from an assumed diplotype.
        mapping = {
            "poor_metabolizer": 0.0,
            "intermediate_metabolizer": 0.5,
            "normal_metabolizer": 1.0,
            "rapid_metabolizer": 1.25,
            "ultrarapid_metabolizer": 1.5,
        }
        if explicit not in mapping:
            raise EngineError(
                E_BAD_VALUE, f"unknown {gene} phenotype {explicit!r}", "phenotypes"
            )
        return {
            "gene": gene,
            "diplotype": None,
            "activity_score": None,
            "phenotype": explicit,
            "relative_activity": mapping[explicit],
            "source": "phenotype supplied by caller",
            "basis": "phenotype",
        }
    return None
