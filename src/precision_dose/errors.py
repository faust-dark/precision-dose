"""Machine-readable error and warning codes.

Every refusal carries a stable code so downstream CDS can branch on it without
parsing prose. Codes are part of the public contract: do not rename them.
"""


class EngineError(Exception):
    """A refusal. Carries a stable code, the offending field, and a message."""

    def __init__(self, code, message, field=None, details=None):
        super().__init__(message)
        self.code = code
        self.message = message
        self.field = field
        self.details = details or {}

    def as_dict(self):
        out = {"code": self.code, "message": self.message}
        if self.field is not None:
            out["field"] = self.field
        if self.details:
            out["details"] = self.details
        return out


# --- input / schema ---------------------------------------------------------
E_SCHEMA = "E_SCHEMA"
E_MISSING_FIELD = "E_MISSING_FIELD"
E_BAD_VALUE = "E_BAD_VALUE"
E_BAD_UNIT = "E_BAD_UNIT"
E_UNKNOWN_DRUG = "E_UNKNOWN_DRUG"
E_UNKNOWN_TARGET = "E_UNKNOWN_TARGET"
E_AGE_OUT_OF_SCOPE = "E_AGE_OUT_OF_SCOPE"

# --- library ----------------------------------------------------------------
E_LIBRARY_INVALID = "E_LIBRARY_INVALID"
E_LIBRARY_UNRESOLVED = "E_LIBRARY_UNRESOLVED"
E_LIBRARY_NO_MODEL = "E_LIBRARY_NO_MODEL"

# --- computation ------------------------------------------------------------
E_NOT_IDENTIFIABLE = "E_NOT_IDENTIFIABLE"
E_NOT_CONVERGED = "E_NOT_CONVERGED"
E_TARGET_UNATTAINABLE = "E_TARGET_UNATTAINABLE"
E_UNSUPPORTED = "E_UNSUPPORTED"

# --- warnings (informational, never blocking) -------------------------------
W_SCr_NOT_STEADY = "W_SCr_NOT_STEADY"
W_CG_EGFR_DISCORDANCE = "W_CG_EGFR_DISCORDANCE"
W_OBESITY_WEIGHT_AMBIGUITY = "W_OBESITY_WEIGHT_AMBIGUITY"
W_MODEL_NO_WEIGHT_V = "W_MODEL_NO_WEIGHT_V"
W_LOW_ALBUMIN_BINDING = "W_LOW_ALBUMIN_BINDING"
W_HEPATIC_DATA_MISSING = "W_HEPATIC_DATA_MISSING"
W_GENOTYPE_INCOMPLETE = "W_GENOTYPE_INCOMPLETE"
W_DDI_PRESENT = "W_DDI_PRESENT"
W_DIALYSIS_ACTIVE = "W_DIALYSIS_ACTIVE"
W_LEVEL_NEAR_DOSE = "W_LEVEL_NEAR_DOSE"
W_WEAK_INFORMATION = "W_WEAK_INFORMATION"
W_TARGET_UNATTAINABLE = "W_TARGET_UNATTAINABLE"
W_PTA_WIDE_CI = "W_PTA_WIDE_CI"
W_LABEL_CRITERIA = "W_LABEL_CRITERIA"
W_INDUCTION_NOT_STEADY = "W_INDUCTION_NOT_STEADY"
W_MIC_ASSUMED = "W_MIC_ASSUMED"
W_SALT_FACTOR = "W_SALT_FACTOR"
W_CURVE_INTERPOLATED = "W_CURVE_INTERPOLATED"
W_NONLINEAR_NO_STEADY_STATE = "W_NONLINEAR_NO_STEADY_STATE"

SEVERITY = {
    W_SCr_NOT_STEADY: "info",
    W_CG_EGFR_DISCORDANCE: "info",
    W_OBESITY_WEIGHT_AMBIGUITY: "info",
    W_MODEL_NO_WEIGHT_V: "caution",
    W_LOW_ALBUMIN_BINDING: "info",
    W_HEPATIC_DATA_MISSING: "caution",
    W_GENOTYPE_INCOMPLETE: "info",
    W_DDI_PRESENT: "info",
    W_DIALYSIS_ACTIVE: "info",
    W_LEVEL_NEAR_DOSE: "info",
    W_WEAK_INFORMATION: "caution",
    W_TARGET_UNATTAINABLE: "caution",
    W_PTA_WIDE_CI: "caution",
    W_LABEL_CRITERIA: "info",
    W_INDUCTION_NOT_STEADY: "caution",
    W_MIC_ASSUMED: "info",
    W_SALT_FACTOR: "info",
    W_CURVE_INTERPOLATED: "info",
    W_NONLINEAR_NO_STEADY_STATE: "caution",
}
