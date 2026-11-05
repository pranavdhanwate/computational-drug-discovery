from typing import Dict, Any, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
import pandas as pd


def compute_molecular_descriptors(mol: Chem.Mol) -> Dict[str, Any]:
    """Calculates key physicochemical descriptors for a given RDKit molecule."""
    if mol is None:
        return {}

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rot_bonds = Lipinski.NumRotatableBonds(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    mr = Crippen.MolMR(mol)

    # Lipinski Rule of Five criteria
    lipinski_violations = 0
    if mw > 500.0:
        lipinski_violations += 1
    if logp > 5.0:
        lipinski_violations += 1
    if hbd > 5:
        lipinski_violations += 1
    if hba > 10:
        lipinski_violations += 1

    lipinski_pass = lipinski_violations <= 1

    # Veber bioavailability criteria
    veber_pass = (rot_bonds <= 10) and (tpsa <= 140.0)

    return {
        "molecular_weight": round(float(mw), 2),
        "logp": round(float(logp), 2),
        "hbd": int(hbd),
        "hba": int(hba),
        "tpsa": round(float(tpsa), 2),
        "rotatable_bonds": int(rot_bonds),
        "heavy_atoms": int(heavy_atoms),
        "fraction_csp3": round(float(fsp3), 3),
        "molar_refractivity": round(float(mr), 2),
        "lipinski_violations": int(lipinski_violations),
        "lipinski_pass": bool(lipinski_pass),
        "veber_pass": bool(veber_pass),
    }


def add_descriptors_to_dataframe(df: pd.DataFrame, smiles_col: str = "canonical_smiles") -> pd.DataFrame:
    """Computes descriptors for all valid SMILES in a dataframe and appends columns."""
    records = []
    for _, row in df.iterrows():
        smiles = row.get(smiles_col)
        if pd.isna(smiles) or not smiles:
            records.append({})
            continue

        mol = Chem.MolFromSmiles(str(smiles))
        desc = compute_molecular_descriptors(mol)
        records.append(desc)

    desc_df = pd.DataFrame(records)
    return pd.concat([df.reset_index(drop=True), desc_df.reset_index(drop=True)], axis=1)
