from typing import List, Dict, Any, Optional
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys


def get_morgan_fingerprint(mol: Chem.Mol, radius: int = 2, n_bits: int = 2048):
    """Generates Morgan / Circular fingerprint (ECFP4 equivalent for radius=2)."""
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)


def get_maccs_keys(mol: Chem.Mol):
    """Generates 166-bit MACCS structural key fingerprint."""
    if mol is None:
        return None
    return MACCSkeys.GenMACCSKeys(mol)


def calculate_tanimoto_similarity(fp1, fp2) -> float:
    """Calculates Tanimoto coefficient between two RDKit bit vectors."""
    if fp1 is None or fp2 is None:
        return 0.0
    return float(DataStructs.TanimotoSimilarity(fp1, fp2))


def screen_similarity_against_reference(
    df: pd.DataFrame,
    reference_smiles: str,
    smiles_col: str = "canonical_smiles",
    radius: int = 2,
    n_bits: int = 2048,
) -> pd.DataFrame:
    """Computes Morgan and MACCS Tanimoto similarity of every molecule against a reference drug."""
    ref_mol = Chem.MolFromSmiles(reference_smiles)
    if ref_mol is None:
        raise ValueError(f"Invalid reference SMILES: {reference_smiles}")

    ref_morgan = get_morgan_fingerprint(ref_mol, radius=radius, n_bits=n_bits)
    ref_maccs = get_maccs_keys(ref_mol)

    morgan_sims = []
    maccs_sims = []

    for _, row in df.iterrows():
        smi = row.get(smiles_col)
        if pd.isna(smi) or not smi:
            morgan_sims.append(0.0)
            maccs_sims.append(0.0)
            continue

        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            morgan_sims.append(0.0)
            maccs_sims.append(0.0)
            continue

        mol_morgan = get_morgan_fingerprint(mol, radius=radius, n_bits=n_bits)
        mol_maccs = get_maccs_keys(mol)

        morgan_sims.append(round(calculate_tanimoto_similarity(ref_morgan, mol_morgan), 4))
        maccs_sims.append(round(calculate_tanimoto_similarity(ref_maccs, mol_maccs), 4))

    result_df = df.copy()
    result_df["tanimoto_morgan_ecfp4"] = morgan_sims
    result_df["tanimoto_maccs"] = maccs_sims
    return result_df
