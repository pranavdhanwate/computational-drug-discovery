import argparse
import sys
from pathlib import Path
from typing import Optional
import pandas as pd
from rdkit import Chem

from .data_loader import MoleculeDataLoader
from .descriptors import add_descriptors_to_dataframe
from .fingerprints import screen_similarity_against_reference
from .molecular_prep import prepare_ligand_3d, clean_receptor_pdb
from .docking import DockingEngine
from .ranking import rank_candidates
from .visualize import generate_structure_grid, plot_lipinski_chemical_space, plot_affinity_ranking


def run_pipeline(
    input_csv: str | Path,
    receptor_pdb: Optional[str | Path] = None,
    reference_smiles: str = "COCCOc1cc2c(Nc3cccc(C#C)c3)ncnc2cc1OCCOC",  # Erlotinib
    output_dir: str | Path = "results",
) -> pd.DataFrame:
    """Executes the full computational drug discovery virtual screening workflow."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = out_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    mols_3d_dir = out_dir / "3d_structures"
    mols_3d_dir.mkdir(parents=True, exist_ok=True)

    # 1. Dataset Loading and SMILES Standardization
    loader = MoleculeDataLoader()
    df = loader.load_csv(input_csv)
    valid_df = df[df["is_valid"]].copy()

    # 2. Molecular Descriptors & Lipinski Rules
    df_desc = add_descriptors_to_dataframe(valid_df, smiles_col="canonical_smiles")

    # 3. Structural Similarity against Reference
    df_sim = screen_similarity_against_reference(
        df_desc,
        reference_smiles=reference_smiles,
        smiles_col="canonical_smiles",
    )

    # 4. Receptor Preparation
    if receptor_pdb and Path(receptor_pdb).exists():
        cleaned_pdb_path = out_dir / "prepared_receptor.pdb"
        clean_receptor_pdb(receptor_pdb, cleaned_pdb_path)

    # 5. 3D Conformation Generation and Docking Scoring
    docking_engine = DockingEngine()
    affinities = []
    kis = []
    ligand_efficiencies = []

    for _, row in df_sim.iterrows():
        smi = row["canonical_smiles"]
        cid = row.get("compound_id", "MOL")
        mol_3d, _ = prepare_ligand_3d(smi, cid, output_dir=mols_3d_dir)

        dock_res = docking_engine.evaluate_affinity_empirical(mol_3d, target_name="EGFR")
        affinities.append(dock_res["binding_affinity_kcal_mol"])
        kis.append(dock_res["estimated_ki_um"])
        ligand_efficiencies.append(dock_res["ligand_efficiency"])

    df_sim["binding_affinity_kcal_mol"] = affinities
    df_sim["estimated_ki_um"] = kis
    df_sim["ligand_efficiency"] = ligand_efficiencies

    # 6. Multi-Criteria Lead Prioritization
    final_ranked = rank_candidates(df_sim)

    # 7. Scientific Outputs & Visualizations
    csv_out = out_dir / "screening_results.csv"
    final_ranked.to_csv(csv_out, index=False)

    generate_structure_grid(final_ranked, figures_dir / "top_leads_grid.png")
    plot_lipinski_chemical_space(final_ranked, figures_dir / "lipinski_chemical_space.png")
    plot_affinity_ranking(final_ranked, figures_dir / "affinity_ranking.png")

    return final_ranked


def main():
    parser = argparse.ArgumentParser(description="Computational Drug Discovery Pipeline")
    parser.add_argument("--input", "-i", type=str, default="data/egfr_compounds.csv", help="Input compounds CSV")
    parser.add_argument("--receptor", "-r", type=str, default="data/1m17.pdb", help="Receptor PDB file")
    parser.add_argument("--output", "-o", type=str, default="results", help="Output directory")
    args = parser.parse_args()

    results = run_pipeline(args.input, args.receptor, output_dir=args.output)
    print(f"Pipeline executed successfully. Top {min(5, len(results))} prioritized candidates:")
    print(results[["priority_rank", "name", "binding_affinity_kcal_mol", "composite_score", "lead_classification"]].head(5).to_string(index=False))


if __name__ == "__main__":
    main()
