import pytest
from pathlib import Path
from src.molecular_prep import prepare_ligand_3d, clean_receptor_pdb


def test_prepare_ligand_3d(tmp_path):
    smi = "CC(=O)Nc1ccc(O)cc1"  # Paracetamol
    mol, sdf_path = prepare_ligand_3d(smi, "PARACETAMOL", output_dir=tmp_path)

    assert mol is not None
    assert mol.GetNumConformers() > 0
    assert sdf_path.exists()
    assert sdf_path.stat().st_size > 0


def test_clean_receptor_pdb(tmp_path):
    sample_pdb_text = (
        "ATOM      1  N   MET A   1      20.154  12.350  14.210  1.00 20.00           N\n"
        "ATOM      2  CA  MET A   1      21.230  13.200  15.100  1.00 20.00           C\n"
        "HETATM 1001  O   HOH A 201      25.100  10.200  18.400  1.00 30.00           O\n"
        "END\n"
    )
    raw_pdb = tmp_path / "raw.pdb"
    raw_pdb.write_text(sample_pdb_text)

    clean_pdb = tmp_path / "clean.pdb"
    summary = clean_receptor_pdb(raw_pdb, clean_pdb)

    assert summary["waters_removed"] == 1
    assert summary["protein_atom_lines"] == 3  # ATOM 1, ATOM 2, END
    assert clean_pdb.exists()
