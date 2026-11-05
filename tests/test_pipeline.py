import pytest
from pathlib import Path
import pandas as pd
from src.pipeline import run_pipeline


def test_full_pipeline_run(tmp_path):
    project_root = Path(__file__).resolve().parent.parent
    data_path = project_root / "data" / "egfr_compounds.csv"
    pdb_path = project_root / "data" / "1m17.pdb"

    results = run_pipeline(
        input_csv=data_path,
        receptor_pdb=pdb_path,
        output_dir=tmp_path,
    )

    assert isinstance(results, pd.DataFrame)
    assert len(results) > 0
    assert "priority_rank" in results.columns
    assert "binding_affinity_kcal_mol" in results.columns
    assert (tmp_path / "screening_results.csv").exists()
    assert (tmp_path / "figures" / "lipinski_chemical_space.png").exists()
    assert (tmp_path / "figures" / "affinity_ranking.png").exists()
