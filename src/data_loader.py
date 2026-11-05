from pathlib import Path
from typing import List, Dict, Any, Optional
import pandas as pd
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover


class MoleculeDataLoader:
    """Loads, validates, and standardizes molecular datasets."""

    def __init__(self):
        self.salt_remover = SaltRemover()

    def load_csv(self, filepath: str | Path, smiles_col: str = "smiles", id_col: str = "compound_id") -> pd.DataFrame:
        path = Path(filepath)
        if not path.exists():
            raise FileNotFoundError(f"Dataset not found at: {path}")

        df = pd.read_csv(path)
        if smiles_col not in df.columns:
            raise ValueError(f"Column '{smiles_col}' not found in CSV.")

        standardized_records = []
        for idx, row in df.iterrows():
            raw_smiles = str(row[smiles_col]).strip()
            mol_id = str(row[id_col]) if id_col in df.columns else f"MOL_{idx+1}"

            clean_smiles, mol, is_valid, err_msg = self.validate_and_clean_smiles(raw_smiles)
            record = row.to_dict()
            record["canonical_smiles"] = clean_smiles
            record["is_valid"] = is_valid
            record["validation_error"] = err_msg
            standardized_records.append(record)

        return pd.DataFrame(standardized_records)

    def validate_and_clean_smiles(self, smiles: str) -> tuple[Optional[str], Optional[Chem.Mol], bool, Optional[str]]:
        if not smiles or smiles.lower() == "nan":
            return None, None, False, "Empty or NaN SMILES"

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None, None, False, "Invalid chemical structure syntax"

            # Remove counterions / salts and keep the parent molecule
            stripped_mol = self.salt_remover.StripMol(mol)
            if stripped_mol.GetNumAtoms() == 0:
                stripped_mol = mol

            canonical_smiles = Chem.MolToSmiles(stripped_mol, canonical=True, isomericSmiles=True)
            clean_mol = Chem.MolFromSmiles(canonical_smiles)

            return canonical_smiles, clean_mol, True, None
        except Exception as e:
            return None, None, False, str(e)
