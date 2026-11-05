import sys
from pathlib import Path

# Add project root to sys.path
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from src.pipeline import run_pipeline

if __name__ == "__main__":
    data_path = project_root / "data" / "egfr_compounds.csv"
    pdb_path = project_root / "data" / "1m17.pdb"
    output_dir = project_root / "results"

    print("Running Computational Drug Discovery Screening Demo...")
    results = run_pipeline(
        input_csv=data_path,
        receptor_pdb=pdb_path,
        output_dir=output_dir,
    )
    print("\nScreening Complete! Results preview:")
    cols = ["priority_rank", "name", "molecular_weight", "logp", "binding_affinity_kcal_mol", "composite_score", "lead_classification"]
    print(results[cols].head(8).to_string(index=False))
