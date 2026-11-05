import pytest
import pandas as pd
from src.ranking import rank_candidates


def test_ranking_logic():
    data = {
        "compound_id": ["C1", "C2", "C3"],
        "name": ["PotentLead", "WeakAnalog", "NonDrugControl"],
        "binding_affinity_kcal_mol": [-9.5, -6.0, -2.0],
        "tanimoto_morgan_ecfp4": [0.85, 0.40, 0.05],
        "ligand_efficiency": [0.38, 0.25, 0.12],
        "lipinski_violations": [0, 1, 3],
    }
    df = pd.DataFrame(data)
    ranked = rank_candidates(df)

    assert ranked.iloc[0]["name"] == "PotentLead"
    assert ranked.iloc[0]["priority_rank"] == 1
    assert "Tier 1" in ranked.iloc[0]["lead_classification"]
    assert ranked.iloc[-1]["name"] == "NonDrugControl"
