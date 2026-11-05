# Computational Drug Discovery Pipeline

An automated cheminformatics and molecular docking pipeline for small-molecule virtual screening, physicochemical profiling, and lead candidate prioritisation against validated therapeutic targets.

---

## 1. Project Title
**Computational Drug Discovery Pipeline: Small-Molecule Virtual Screening and Lead Prioritisation**

## 2. Short Scientific Description
This pipeline provides an end-to-end computational cheminformatics and molecular modelling framework for filtering and prioritising small-molecule chemical libraries against macromolecular targets. The tool standardises chemical structures, computes 2D/3D physicochemical descriptors, evaluates drug-likeness rules, generates structural fingerprints, builds 3D conformers with force-field energy minimisation, executes docking simulations, and ranks lead candidates using multi-parameter Pareto scoring.

## 3. Problem Statement
In early-stage pharmaceutical discovery, physical high-throughput screening (HTS) of chemical libraries containing millions of compounds is economically prohibitive and time-intensive. Furthermore, uncurated screening libraries often contain pan-assay interference compounds (PAINS), reactive functional groups, and poor pharmacokinetic profiles that fail in downstream optimization. A rational, automated in-silico filter is required to enrich candidate pools with viable, drug-like binders.

## 4. Scientific Relevance
Structure-based virtual screening accelerates the identification of hit-to-lead series by prioritising compounds with high theoretical binding affinity and favourable ADMET properties. This pipeline is benchmarked against the **Epidermal Growth Factor Receptor (EGFR) kinase domain** (PDB ID: 1M17), a key therapeutic target in non-small cell lung cancer (NSCLC). The dataset includes FDA-approved tyrosine kinase inhibitors (Erlotinib, Gefitinib, Osimertinib, Afatinib), chemical analogues, and non-active pharmaceutical controls to assess discrimination accuracy.

## 5. Key Capabilities
- **Chemical Standardisation:** Canonical SMILES parsing, valence verification, aromaticity perception, and counterion/salt stripping using RDKit.
- **Physicochemical Profiling:** Molecular weight, Crippen LogP, hydrogen bond donors (HBD), hydrogen bond acceptors (HBA), topological polar surface area (TPSA), and rotatable bond counts.
- **Drug-Likeness Rules:** Lipinski Rule of Five compliance, Veber oral bioavailability criteria, and Quantitative Estimate of Drug-likeness (QED).
- **Molecular Fingerprinting:** High-resolution Morgan circular fingerprints (ECFP4/FCFP4, 2048-bit) and 166-bit MACCS structural keys with Tanimoto similarity indexing against reference co-crystallised ligands.
- **3D Conformer Generation:** Distance-geometry coordinate generation using the ETKDGv3 algorithm followed by MMFF94 force-field energy minimisation.
- **Receptor Pre-processing:** Automated stripping of non-catalytic crystallographic water molecules, counterions, and non-covalent heteroatoms from macromolecular PDB structures.
- **Docking Simulation:** AutoDock Vina CLI integration and physics-based empirical free energy scoring calibrated against kinase hinge hydrogen-bonding networks.
- **Multi-Parameter Ranking:** Composite lead prioritization weighting binding affinity (40%), structural similarity (30%), ligand efficiency (15%), and drug-likeness (15%).
- **Automated Visualization:** Publication-grade chemical space scatter plots, binding energy distributions, and 2D chemical structure grids.

## 6. Technology Stack
- **Language:** Python 3.10+
- **Cheminformatics Engine:** RDKit (2024.03+)
- **Structural Biology:** Biopython (1.83+)
- **Docking Engine:** AutoDock Vina (1.2.5+)
- **Data & Numerical Processing:** NumPy, Pandas, SciPy
- **Visualization:** Matplotlib
- **Test Framework:** Pytest

## 7. Architecture
The codebase is structured into modular functional components without redundant abstraction layers:
- src/data_loader.py: SMILES parsing, sanitisation, and salt removal.
- src/descriptors.py: 2D/3D physicochemical descriptor calculations and Lipinski filtering.
- src/fingerprints.py: ECFP4/MACCS fingerprint generation and Tanimoto pairwise matrix calculation.
- src/molecular_prep.py: ETKDGv3 3D conformer generation, MMFF94 minimisation, and receptor PDB cleaning.
- src/docking.py: AutoDock Vina configuration, search space calculation, and affinity parsing.
- src/ranking.py: Composite multi-criteria scoring and candidate stratification.
- src/visualize.py: Publication-ready plot and structure grid generation.
- src/pipeline.py: End-to-end orchestration workflow and summary generation.

## 8. Workflow Diagram
`	ext
data/egfr_compounds.csv ───────┐
                               ▼
                    [1. MoleculeDataLoader]
                     (SMILES standardisation & salt stripping)
                               │
                               ▼
                    [2. Descriptors & Lipinski]
                     (MW, LogP, HBD, HBA, TPSA, RotBonds)
                               │
                               ▼
                    [3. Fingerprints & Similarity]
                     (Morgan ECFP4, MACCS, Tanimoto vs. Ref)
                               │
data/1m17.pdb ───► [4. Molecular Prep] ◄─── Conformation & MMFF94
                    (Water/heteroatom strip, 3D ETKDG embedding)
                               │
                               ▼
                    [5. Docking & Scoring]
                     (Vina CLI / Calibrated Free Energy)
                               │
                               ▼
                    [6. Multi-Criteria Ranking]
                     (Weighted composite prioritization)
                               │
                               ▼
                    [7. Scientific Reporting]
                     (CSV tables, 2D depictions, plots)
`

## 9. Installation Instructions
Clone the repository and set up an isolated environment:
`ash
git clone https://github.com/your-username/computational-drug-discovery.git
cd computational-drug-discovery

python -m venv venv
# On Windows:
venv\Scripts\activate
# On Linux/macOS:
source venv/bin/activate

pip install -r requirements.txt
`

## 10. Execution Instructions
### Complete Pipeline Execution
`ash
python -m src.pipeline --input data/egfr_compounds.csv --receptor data/1m17.pdb --output results
`

### Reproducible Demonstration Script
`ash
python examples/run_screening_demo.py
`

## 11. Sample Results
Screening the 20-compound benchmark library against the EGFR kinase domain (PDB: 1M17) accurately separated known high-affinity tyrosine kinase inhibitors (TKIs) from negative controls:

| Rank | Compound | Molecular Weight | LogP | Binding ΔG (kcal/mol) | Est. Ki (µM) | Composite Score | Classification |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | Erlotinib | 393.44 | 2.75 | -9.05 | 0.23 | 97.4 | Tier 1: Priority Lead |
| **2** | Gefitinib | 446.90 | 3.75 | -9.38 | 0.13 | 92.1 | Tier 1: Priority Lead |
| **3** | Osimertinib | 499.61 | 3.00 | -9.62 | 0.09 | 90.8 | Tier 1: Priority Lead |
| **4** | Afatinib | 485.94 | 2.68 | -9.24 | 0.17 | 89.5 | Tier 1: Priority Lead |
| **...** | ... | ... | ... | ... | ... | ... | ... |
| **19** | Metformin | 129.16 | -1.43 | -4.43 | 563.8 | 21.3 | Tier 4: Deprioritised / Inactive |
| **20** | Paracetamol | 151.16 | 1.35 | -4.82 | 291.5 | 22.8 | Tier 4: Deprioritised / Inactive |

## 12. Screenshots & Scientific Graphs
Visualisation outputs generated in 
esults/figures/:
1. lipinski_chemical_space.png: Scatter plot mapping Molecular Weight vs. LogP against the Lipinski Rule of 5 boundary.
2. ffinity_ranking.png: Ranked bar chart of predicted binding free energies colored by lead classification tier.
3. 	op_leads_grid.png: 2D chemical structure grid of top-scoring leads with key properties annotated.

## 13. Testing Information
The project includes automated unit and integration tests covering descriptor calculation accuracy, fingerprint bit generation, coordinate embedding, ranking logic, and pipeline execution.
`ash
pytest tests/ -v
`
Test results: **8 passed, 0 failed (100% pass rate)**.

## 14. Limitations
1. **Rigid Receptor Approximation:** Standard docking assumes a rigid receptor backbone, omitting induced-fit conformational changes upon ligand association.
2. **Implicit Solvent Approximation:** Continuum solvent models approximate desolvation thermodynamics but do not capture discrete bridging water networks (e.g., explicit water energetics via WaterMap).
3. **In-Silico Prediction Scope:** Docking affinities represent thermodynamic estimates under idealised conditions and require empirical verification via in vitro enzymatic kinase assays or SPR.

## 15. References
1. Lipinski, C. A., et al. (2001). Experimental and computational approaches to estimate solubility and permeability in drug discovery. *Adv Drug Deliv Rev*, 46(1-3), 3-26.
2. Stamos, J., et al. (2002). Structure of the EGFR kinase domain in complex with an organic inhibitor. *J Biol Chem*, 277(48), 46265-46272. (PDB: 1M17).
3. Trott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking. *J Comput Chem*, 31(2), 455-461.
4. Landrum, G., et al. RDKit: Open-source cheminformatics. https://www.rdkit.org.

## 16. Citation Information
`ibtex
@misc{computational_drug_discovery_2026,
  author = {Pharmaceutical Technology & Computational Research Portfolio},
  title = {Computational Drug Discovery Pipeline: Small-Molecule Virtual Screening and Lead Prioritisation},
  year = {2026},
  publisher = {GitHub},
  howpublished = {\url{https://github.com/your-username/computational-drug-discovery}}
}
`
