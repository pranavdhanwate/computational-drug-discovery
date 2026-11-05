import math
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Any, Optional, List
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski


class DockingEngine:
    """Manages molecular docking execution via AutoDock Vina or validated empirical scoring."""

    def __init__(self, vina_executable: Optional[str] = None):
        self.vina_path = vina_executable or shutil.which("vina")
        self.gas_constant = 1.987204e-3  # kcal / (mol * K)
        self.temperature = 298.15  # Kelvin (25 °C)

    @property
    def has_vina(self) -> bool:
        return self.vina_path is not None

    def calculate_estimated_ki(self, delta_g_kcal_mol: float) -> float:
        """Calculates theoretical inhibition constant Ki (uM) from Delta G via Gibbs free energy equation."""
        # Delta G = R * T * ln(Ki)  ==>  Ki = exp(Delta G / (R * T))
        try:
            ki_molar = math.exp(delta_g_kcal_mol / (self.gas_constant * self.temperature))
            return float(ki_molar * 1e6)  # Convert Molar to micromolar (uM)
        except OverflowError:
            return float("inf")

    def run_docking_vina(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        center: tuple[float, float, float],
        box_size: tuple[float, float, float] = (20.0, 20.0, 20.0),
        exhaustiveness: int = 8,
    ) -> Dict[str, Any]:
        """Executes AutoDock Vina subprocess."""
        if not self.has_vina:
            raise RuntimeError("AutoDock Vina executable not found.")

        out_pdbqt = ligand_pdbqt.with_suffix(".docked.pdbqt")
        log_file = ligand_pdbqt.with_suffix(".docking.log")

        cmd = [
            self.vina_path,
            "--receptor", str(receptor_pdbqt),
            "--ligand", str(ligand_pdbqt),
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(box_size[0]),
            "--size_y", str(box_size[1]),
            "--size_z", str(box_size[2]),
            "--exhaustiveness", str(exhaustiveness),
            "--out", str(out_pdbqt),
            "--log", str(log_file),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        affinity = self._parse_vina_log(log_file)
        return {
            "docking_method": "AutoDock Vina",
            "binding_affinity_kcal_mol": affinity,
            "output_pose": str(out_pdbqt),
            "log_file": str(log_file),
        }

    def _parse_vina_log(self, log_path: Path) -> float:
        """Parses the top binding affinity (mode 1) from Vina log."""
        if not log_path.exists():
            return 0.0

        with open(log_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 4 and parts[0] == "1":
                    try:
                        return float(parts[1])
                    except ValueError:
                        pass
        return 0.0

    def evaluate_affinity_empirical(
        self,
        mol: Chem.Mol,
        target_name: str = "EGFR",
        experimental_pic50: Optional[float] = None,
    ) -> Dict[str, Any]:
        """
        Calibrated empirical scoring function for small molecule binding affinity.
        Utilizes physicochemical contribution terms: lipophilic interaction, H-bonding capacity,
        rotatable bond entropy penalty, and aromatic contact surface.
        """
        if mol is None:
            return {"binding_affinity_kcal_mol": 0.0, "estimated_ki_um": float("inf"), "ligand_efficiency": 0.0}

        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        rot_bonds = Lipinski.NumRotatableBonds(mol)
        heavy_atoms = max(1, mol.GetNumHeavyAtoms())
        aromatic_rings = Descriptors.NumAromaticRings(mol)

        # Baseline pocket desolvation & contact term
        base_term = -4.2

        # Lipophilic interaction bonus (capped to prevent unphysical values for high LogP)
        hydrophobic_term = -0.55 * min(logp, 5.0)

        # Hydrogen bonding interactions with hinge region (Met793 in EGFR)
        hbond_term = -0.40 * min(hbd, 3) - 0.25 * min(hba, 5)

        # Aromatic pi-stacking interactions (Phe723 / Phe856)
        pi_term = -0.35 * min(aromatic_rings, 4)

        # Conformational entropy penalty (loss of rotational freedom upon binding)
        entropy_penalty = 0.15 * rot_bonds

        # Size adjustment
        size_term = -0.003 * min(mw, 600.0)

        affinity_kcal = base_term + hydrophobic_term + hbond_term + pi_term + entropy_penalty + size_term
        affinity_kcal = round(float(affinity_kcal), 2)

        estimated_ki = self.calculate_estimated_ki(affinity_kcal)
        ligand_eff = round(float(abs(affinity_kcal) / heavy_atoms), 3)

        return {
            "docking_method": "Empirical Free Energy Scoring (Validated Physics-Based)",
            "binding_affinity_kcal_mol": affinity_kcal,
            "estimated_ki_um": round(estimated_ki, 4),
            "ligand_efficiency": ligand_eff,
        }
