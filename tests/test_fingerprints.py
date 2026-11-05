import pytest
from rdkit import Chem
from src.fingerprints import get_morgan_fingerprint, calculate_tanimoto_similarity


def test_identical_molecule_tanimoto():
    # Identical molecules must yield Tanimoto coefficient exactly 1.0
    smi = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"  # Gefitinib
    mol1 = Chem.MolFromSmiles(smi)
    mol2 = Chem.MolFromSmiles(smi)

    fp1 = get_morgan_fingerprint(mol1)
    fp2 = get_morgan_fingerprint(mol2)

    similarity = calculate_tanimoto_similarity(fp1, fp2)
    assert pytest.approx(similarity, abs=1e-5) == 1.0


def test_dissimilar_molecule_tanimoto():
    # Gefitinib vs Metformin
    smi1 = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"
    smi2 = "CN(C)C(=N)N=C(N)N"

    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)

    fp1 = get_morgan_fingerprint(mol1)
    fp2 = get_morgan_fingerprint(mol2)

    similarity = calculate_tanimoto_similarity(fp1, fp2)
    assert 0.0 <= similarity < 0.25
