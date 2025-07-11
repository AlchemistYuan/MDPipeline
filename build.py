import sys
import os
import argparse
import urllib
from pdbfixer import PDBFixer
from openmm import unit
from openmm.app import Modeller, PDBFile, ForceField, PME
from openmm import XmlSerializer
from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit import Molecule, Topology
import MDAnalysis as mda
from rdkit import Chem
from rdkit.Chem import AllChem
import parmed as pmd
from parmed.tools.actions import changeRadii


"""
This script builds an OpenMM-ready system given a PDB ID.
It adds missing atoms and hydrogens, assigns AMBER ff14SB + OL15
parameters for macromolecules (e.g. RNA/DNA or protein) and GAFF2 for ligands.
It then solvates the system by adding water and neutralizes the system
by adding ions.
In the end, it saves the serialized system and coordinate file.

Usage:
    python build.py --pdb [PDB ID] --lig [LIG RESNAME]
"""


def read_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Build OpenMM system from a PDB ID')
    parser.add_argument('--pdb', help='PDB ID to fetch and build', required=True)
    parser.add_argument('--lig', help='The residue name of the ligand', default='LIG')
    parser.add_argument('--pH', type=float, default=7.0,
                        help='Protonation pH for adding hydrogens')
    parser.add_argument('--padding', type=float, default=1.0,
                        help='Water padding distance in nm')
    parser.add_argument('--ionic_strength', type=float, default=0.15,
                        help='Salt concentration in molar')
    parser.add_argument('--inputdir', type=str, default='input',
                        help='The directory for the input PDB file')
    args = parser.parse_args()
    return args


def fetch_pdb(pdbid: str, inputdir: str) -> str:
    """
    Download PDB file from RCSB and save to input_dir.

    Parameters
    ----------
    pdbid : str
        The PDB ID to be searched or downloaded
    inputidr : str
        The directory for the downloaded PDB file

    Returns
    -------
    pab_path : str
        The path to the pdb file
    """
    os.makedirs(inputdir, exist_ok=True)
    pdb_path = os.path.join(inputdir, f"{pdbid}.pdb")
    url = f"https://files.rcsb.org/download/{pdbid}.pdb"
    print(f"Fetching {pdbid} from {url}...")
    urllib.request.urlretrieve(url, pdb_path)
    print(f"Saved {pdbid} to {pdb_path}")
    return pdb_path


def fetch_ligand(ligname: str, inputdir: str) -> str:
    """
    Download SDF file from RCSB and save to inputdir.

    Parameters
    ----------
    ligname : str
        The ligand residue name to be searched or downloaded
    inputidr : str
        The directory for the downloaded PDB file

    Returns
    -------
    pab_path : str
        The path to the pdb file
    """
    os.makedirs(inputdir, exist_ok=True)
    sdf_path = os.path.join(inputdir, f"{ligname}_real.sdf")
    #url = f"https://files.rcsb.org/ligands/download/{ligname}_real.sdf"
    #print(f"Fetching {ligname} from {url}...")
    os.system(f'curl "https://files.rcsb.org/ligands/download/{ligname}_ideal.sdf" -o {sdf_path}')
    #urllib.request.urlretrieve(url, sdf_path)
    print(f"Saved {ligname} to {sdf_path}")
    return sdf_path


def write_lig_only_sdf_from_pdb(pdb_file: str, ligname: str, inputdir: str) -> Molecule:
    sdf_path = f'{inputdir}/{ligname}_real.sdf'
    if not os.path.isfile(sdf_path):
        sdf_path = fetch_ligand(ligname, inputdir)
    lig_template = Chem.MolFromMolFile(sdf_path, removeHs=False)
    u = mda.Universe(pdb_file)
    u_lig = u.select_atoms(f'resname {ligname}')
    u_lig.atoms.write(f'{pdb_file.rsplit(".", 1)[0]}_{ligname}.pdb')
    mol = Chem.MolFromPDBFile(f'{pdb_file.rsplit(".", 1)[0]}_{ligname}.pdb',
                              removeHs=False, sanitize=False)
    mol = Chem.AllChem.AssignBondOrdersFromTemplate(lig_template, mol)
    Chem.SanitizeMol(mol)  # standard valence/aromaticity checks
    Chem.rdmolops.AssignAtomChiralTagsFromStructure(
        mol)  # uses 3-D to mark R/S stereocentres :contentReference[oaicite:0]{index=0}
    Chem.rdmolops.DetectBondStereoChemistry(mol,
                                            mol.GetConformer())  # sets E/Z on double bonds :contentReference[oaicite:1]{index=1}
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    w = Chem.SDWriter(f'{pdb_file.rsplit(".", 1)[0]}_{ligname}.sdf')
    w.write(mol)
    w.close()

    # Load the definition of the small molecule in the system from an SDF file
    ligand = Molecule.from_file(f'{pdb_file.rsplit(".", 1)[0]}_{ligname}.sdf')
    print('Total charge of the ligand:', ligand.total_charge)
    return ligand


def main() -> int:
    args = read_arguments()

    # Determine PDB file path
    # check for an existing PDB in an input directory.
    # only fetch from RCSB if it’s missing
    # saving it locally before reading into OpenMM
    pdb_file = os.path.join(args.inputdir, f"{args.pdb}.pdb")
    if os.path.isfile(pdb_file):
        print(f"Using local PDB file: {pdb_file}")
        fixer = PDBFixer(filename=pdb_file)
    else:
        pdb_file = fetch_pdb(args.pdb, args.inputdir)
        fixer = PDBFixer(filename=pdb_file)
    ligand = write_lig_only_sdf_from_pdb(pdb_file, args.lig, args.inputdir)

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=args.pH)

    modeller = Modeller(fixer.topology, fixer.positions)

    # Define force fields: protein/DNA, ligand (GAFF2), and water
    forcefield = ForceField("amber/protein.ff14SB.xml",
        "amber/RNA.OL3.xml",
        "amber/tip3p_standard.xml",
        "amber/tip3p_HFE_multivalent.xml")
    gaff = GAFFTemplateGenerator(forcefield='gaff-2.11', cache=f'{args.lig}.xml',
                                 molecules=ligand)
    forcefield.registerTemplateGenerator(gaff.generator)
    # Solvate & neutralize
    modeller.addSolvent(
        forcefield,
        model='tip3p',
        padding=args.padding * unit.nanometer,
        ionicStrength=args.ionic_strength * unit.molar,
        neutralize=True,
    )

    # Create the OpenMM System
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        constraints=None,
        rigidWater=False,
        removeCMMotion=False
    )

    # Write Amber prmtop and inpcrd for MMPBSA calculations
    structure = pmd.openmm.load_topology(modeller.topology, system=system, xyz=modeller.positions)
    action = changeRadii(structure, 'mbondi3')
    action.execute()
    structure.save(f'output/{args.pdb}_solvated.prmtop', overwrite=True)
    structure.save(f'output/{args.pdb}_solvated.inpcrd', format='rst7', overwrite=True)
    print(f"Wrote Amber system to: {args.pdb}_solvated.prmtop and {args.pdb}_solvated.inpcrd")
    complex_strip = structure['!:HOH,NA,CL']
    complex_strip.save(f'output/complex_strip.prmtop', overwrite=True)
    complex_strip.save(f'output/complex_strip.inpcrd', format='rst7', overwrite=True)
    print(f"Wrote stripped complex to: complex_strip.prmtop and complex_strip.inpcrd")
    receptor = complex_strip[f'!:{args.lig}']
    receptor.save(f'output/receptor.prmtop', overwrite=True)
    receptor.save(f'output/receptor.inpcrd', format='rst7', overwrite=True)
    print(f"Wrote receptor-only system to: receptor.prmtop and receptor.inpcrd")
    lig_amber = structure[f':{args.lig}']
    lig_amber.save(f'output/ligand.prmtop', overwrite=True)
    lig_amber.save(f'output/ligand.inpcrd', format='rst7', overwrite=True)
    print(f"Wrote ligand-only system to: ligand.prmtop and ligand.inpcrd")
    
    # Serialize and save the system to XML
    xml = XmlSerializer.serialize(system)
    xml_file = f"{args.pdb}_system.xml"
    with open(f'output/{xml_file}', 'w') as f:
        f.write(xml)
    print(f"Wrote system XML to: {xml_file}")

    # Write out the solvated structure as PDB
    pdb_file = f"{args.pdb}_solvated.pdb"
    with open(f'output/{pdb_file}', 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f, keepIds=True)
    print(f"Wrote solvated structure to: {pdb_file}")
    return 0


if __name__ == '__main__':
    sys.exit(main())

