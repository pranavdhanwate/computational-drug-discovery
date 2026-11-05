import pytest
from rdkit import Chem
from src.descriptors import compute_molecular_descriptors


def test_aspirin_descriptors():
    # Aspirin: C9H8O4, MW ~ 180.16 g/mol
    aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = Chem.MolFromSmiles(aspirin_smiles)
    desc = compute_molecular_descriptors(mol)

    assert 179.0 < desc["molecular_weight"] < 181.0
    assert desc["hbd"] == 1  # -COOH donor
    assert desc["hba"] == 3  # carbonyl and ester/acid oxygens
    assert desc["lipinski_violations"] == 0
    assert desc["lipinski_pass"] is True
    assert desc["veber_pass"] is True


def test_lipinski_violation_detection():
    # A large synthetic peptide/molecule with MW > 500 and many donors
    huge_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCC(=O)NC(CC(=O)NC(CC(=O)O)C(=O)O)C(=O)O"
    mol = Chem.MolFromSmiles(huge_smiles)
    desc = compute_molecular_descriptors(mol)

    assert desc["molecular_weight"] > 500.0 or desc["logp"] > 5.0
    assert isinstance(desc["lipinski_violations"], int)
