from typing import List, Dict, Any
import numpy as np
import pandas as pd


def rank_candidates(
    df: pd.DataFrame,
    affinity_col: str = "binding_affinity_kcal_mol",
    lipinski_col: str = "lipinski_violations",
    similarity_col: str = "tanimoto_morgan_ecfp4",
    le_col: str = "ligand_efficiency",
) -> pd.DataFrame:
    """
    Ranks candidates using Multi-Criteria Decision Analysis (MCDA).
    Balances binding affinity (40%), structural similarity to reference (30%),
    ligand efficiency (15%), and drug-likeness/Lipinski adherence (15%).
    """
    ranked = df.copy()

    # Normalize binding affinity (lower/more negative kcal/mol is better)
    affinities = ranked[affinity_col].fillna(0.0).values
    min_aff = np.min(affinities) if len(affinities) > 0 else -10.0
    max_aff = np.max(affinities) if len(affinities) > 0 else 0.0
    aff_range = max_aff - min_aff if (max_aff - min_aff) != 0 else 1.0
    norm_affinity = (max_aff - affinities) / aff_range  # 1.0 = best affinity, 0.0 = worst

    # Similarity score (0.0 to 1.0)
    norm_sim = ranked[similarity_col].fillna(0.0).values

    # Ligand efficiency score
    le_vals = ranked[le_col].fillna(0.0).values
    max_le = np.max(le_vals) if np.max(le_vals) > 0 else 1.0
    norm_le = le_vals / max_le

    # Lipinski penalty (0 violations = 1.0, 1 violation = 0.75, 2 violations = 0.5, >=3 = 0.0)
    lip_violations = ranked[lipinski_col].fillna(4).values
    norm_lipinski = np.maximum(0.0, 1.0 - (lip_violations * 0.25))

    # Composite Score calculation (0 - 100)
    composite_score = (
        (0.40 * norm_affinity)
        + (0.30 * norm_sim)
        + (0.15 * norm_le)
        + (0.15 * norm_lipinski)
    ) * 100.0

    ranked["composite_score"] = np.round(composite_score, 2)
    ranked["priority_rank"] = ranked["composite_score"].rank(ascending=False, method="dense").astype(int)

    def assign_category(score: float, violations: int) -> str:
        if score >= 75.0 and violations <= 1:
            return "Tier 1: Priority Lead"
        elif score >= 50.0 and violations <= 1:
            return "Tier 2: Validated Analog"
        elif score >= 35.0:
            return "Tier 3: Moderate Potential"
        else:
            return "Tier 4: Deprioritized / Inactive"

    ranked["lead_classification"] = [
        assign_category(s, v)
        for s, v in zip(ranked["composite_score"], ranked[lipinski_col])
    ]

    return ranked.sort_values(by="priority_rank", ascending=True).reset_index(drop=True)
