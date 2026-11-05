from pathlib import Path
from typing import Optional, Tuple, Dict, Any
from rdkit import Chem
from rdkit.Chem import AllChem


def prepare_ligand_3d(
    smiles: str,
    compound_id: str,
    output_dir: Optional[Path] = None,
    random_seed: int = 42,
) -> tuple[Optional[Chem.Mol], Optional[Path]]:
    """Converts 2D SMILES into energy-minimized 3D conformation using ETKDGv3 and MMFF94."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None

    # Add explicit hydrogens required for 3D stereochemistry and force-field physics
    mol_h = Chem.AddHs(mol)

    # Generate 3D coordinates using ETKDG v3 algorithm
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    embed_status = AllChem.EmbedMolecule(mol_h, params)

    if embed_status != 0:
        # Fallback embedding
        AllChem.EmbedMolecule(mol_h, randomSeed=random_seed, useRandomCoords=True)

    # Energy minimize conformation using Merck Molecular Force Field (MMFF94)
    if AllChem.MMFFHasAllMoleculeParams(mol_h):
        AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)

    output_path = None
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{compound_id}_3d.sdf"
        writer = Chem.SDWriter(str(output_path))
        writer.write(mol_h)
        writer.close()

    return mol_h, output_path


def clean_receptor_pdb(input_pdb: str | Path, output_pdb: str | Path) -> Dict[str, Any]:
    """Cleans receptor PDB by stripping crystallographic water and non-protein heteroatoms."""
    input_path = Path(input_pdb)
    output_path = Path(output_pdb)

    if not input_path.exists():
        raise FileNotFoundError(f"Receptor PDB not found at: {input_path}")

    cleaned_lines = []
    heteroatoms_removed = 0
    waters_removed = 0
    ca_coordinates = []

    with open(input_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ATOM"):
                cleaned_lines.append(line)
                if line[12:16].strip() == "CA":
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        ca_coordinates.append((x, y, z))
                    except ValueError:
                        pass
            elif line.startswith("HETATM"):
                res_name = line[17:20].strip()
                if res_name in ["HOH", "WAT", "H2O", "DOD"]:
                    waters_removed += 1
                else:
                    heteroatoms_removed += 1
            elif line.startswith("TER") or line.startswith("END"):
                cleaned_lines.append(line)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as out:
        out.writelines(cleaned_lines)

    center = (0.0, 0.0, 0.0)
    if ca_coordinates:
        import numpy as np
        arr = np.array(ca_coordinates)
        center = tuple(np.round(np.mean(arr, axis=0), 2))

    return {
        "output_path": str(output_path),
        "protein_atom_lines": len(cleaned_lines),
        "waters_removed": waters_removed,
        "heteroatoms_removed": heteroatoms_removed,
        "geometric_center": center,
    }
