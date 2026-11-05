from pathlib import Path
from typing import List, Optional
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw


def generate_structure_grid(
    df: pd.DataFrame,
    output_path: str | Path,
    max_mols: int = 12,
    smiles_col: str = "canonical_smiles",
    name_col: str = "name",
) -> Path:
    """Generates a 2D chemical structure grid PNG for top ranked molecules."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    mols = []
    legends = []
    subset = df.dropna(subset=[smiles_col]).head(max_mols)

    for _, row in subset.iterrows():
        smi = str(row[smiles_col])
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mols.append(mol)
            name = str(row.get(name_col, "Mol"))
            rank = row.get("priority_rank", "")
            legends.append(f"#{rank} {name}" if rank != "" else name)

    if mols:
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=4,
            subImgSize=(300, 250),
            legends=legends,
            useSVG=False,
        )
        img.save(str(output_path))

    return output_path


def plot_lipinski_chemical_space(df: pd.DataFrame, output_path: str | Path) -> Path:
    """Plots Molecular Weight vs LogP with Lipinski Rule of 5 drug-likeness boundary."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)

    # Lipinski boundary box: MW <= 500, LogP <= 5.0
    ax.axhspan(-2.0, 5.0, xmin=0.0, xmax=500.0 / 700.0, color="#E8F5E9", alpha=0.6, label="Lipinski Compliant Space")
    ax.axvline(500, color="#C62828", linestyle="--", alpha=0.7, label="MW Limit (500 Da)")
    ax.axhline(5.0, color="#E65100", linestyle="--", alpha=0.7, label="LogP Limit (5.0)")

    # Scatter plot
    scatter = ax.scatter(
        df["molecular_weight"],
        df["logp"],
        c=df["composite_score"],
        cmap="viridis",
        s=80,
        edgecolors="black",
        linewidth=0.8,
        alpha=0.9,
    )
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Composite Candidate Score", fontsize=10)

    # Annotate top 5 candidates
    for _, row in df.head(5).iterrows():
        ax.annotate(
            row["name"],
            (row["molecular_weight"], row["logp"]),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=8,
            fontweight="bold",
            color="#0D47A1",
        )

    ax.set_xlabel("Molecular Weight (g/mol)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Calculated LogP (Crippen)", fontsize=11, fontweight="bold")
    ax.set_title("Chemical Space & Drug-Likeness Evaluation (Lipinski Filter)", fontsize=12, fontweight="bold")
    ax.set_xlim(100, 650)
    ax.set_ylim(-1.5, 6.5)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="upper left", fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    return output_path


def plot_affinity_ranking(df: pd.DataFrame, output_path: str | Path) -> Path:
    """Plots binding affinity ranking bar chart across all evaluated candidates."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sorted_df = df.sort_values(by="binding_affinity_kcal_mol", ascending=True)

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    colors = ["#1B5E20" if cat == "Tier 1: Priority Lead" else
              "#2E7D32" if cat == "Tier 2: Validated Analog" else
              "#F57C00" if cat == "Tier 3: Moderate Potential" else
              "#C62828" for cat in sorted_df["lead_classification"]]

    bars = ax.barh(sorted_df["name"], sorted_df["binding_affinity_kcal_mol"], color=colors, edgecolor="black", alpha=0.85)

    ax.axvline(-7.0, color="gray", linestyle="--", linewidth=1.0, label="High Affinity Threshold (-7.0 kcal/mol)")
    ax.set_xlabel("Estimated Binding Free Energy ΔG (kcal/mol)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Candidate Molecule", fontsize=11, fontweight="bold")
    ax.set_title("Predicted Binding Affinity Ranking at EGFR Kinase Domain", fontsize=12, fontweight="bold")
    ax.grid(axis="x", linestyle=":", alpha=0.6)
    ax.legend(loc="lower right", fontsize=9)

    for bar in bars:
        val = bar.get_width()
        ax.text(val - 0.25, bar.get_y() + bar.get_height()/2, f"{val:.1f}", va="center", ha="right", color="white", fontsize=8, fontweight="bold")

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    return output_path
