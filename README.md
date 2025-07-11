# All-atom molecular dynamics simulation and binding free energy calculation for liang-bindng RNA
This is a prototype of and end-to-end, reproducible pipeline for molecular dynamics simulations and binding free energy calculations.

This pipeline performs the following tasks:
- Builds an all-atom simulation system for a ligand binding DNA or RNA given a publicly available PDB entry
- Defines the relevant binding site region
- Estimates absolute binding free energy

## How to reporduce the results in a Linux system
First navigate to the root directory of this pipeline.
1. (Optional) install conda if it's not already installed. I prefer the minimal installer, [miniforge](https://github.com/conda-forge/miniforge) 
2. Create the conda environment: `mamba env create -f env.yml`
3. Activate the conda environment: `conda activate md_env`
4. System setup: `python build.py --pdb 1O15 --lig TEP`
5. Open the jupyter notebook, `mdpipeline.ipynb`, and run through all cells
6. (Optional) Run the last section of `mdpipeline.ipynb` for a hybrid MD/MMPBSA pipeline where the ligand's vdW and bonded terms are parameterized by ANI-2x.

## Runtime analysis
1. Minimization + two-step equilibration runs time: ~ 5 minutes 
2. Production run time: ~8.5 ns / day (or ~2.8 h / ns) on a 8-core CPU 
3. MMPBSA calculation time: ~ 15 min

## Software versions
After creating the conda environment, the exact versions of software can be exported by `mamba env export`

Versions of some key software:
1. ambertools=22.5
2. MDAnalysis=2.7.0
3. numpy=1.26.4
4. openmm=8.2.0
5. openmmforcefields=0.12.0
6. openmm-torch=1.5.1 
7. openmm-ml=1.2 
8. parmed=4.3.0 
9. pdbfixer=1.11 
10. python=3.9 
11. torchani=2.2.4

## Rationale behind the chosen protocol
The protocol is designed to ensure a physically realistic, reproducible workflow from structure preparation through free‐energy estimation.
PDBFixer rapidly fills in missing atoms and hydrogens while Amber (for protein/RNA) and GAFF2 (for ligands) provide well‐validated force fields.
TIP3P solvation with appropriate padding and ionic strength mimics the aqueous environment.
The first equilibration with heavy‐atom restraints maintains the binding pose.
The unrestrained equilibration allows gentle relaxation.
The production run under NPT yields an ensemble of conformations for the subseqquent calculations.
The MMPBSA calculation gives an efficient, end‐point estimate of ligand binding free energy.






