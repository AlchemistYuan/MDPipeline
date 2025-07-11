import argparse
import sys
import json
from typing import List
import numpy as np
import yaml
import MDAnalysis as mda


def read_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Indentify the ligand binding site from a PDB file. '
                    'The binding site is found by selecting residues within certain distacne around the ligand'
                    ' or by user input')
    parser.add_argument('--pdb', help='PDB file to be read', default=None)
    parser.add_argument('--lig', help='The ligand resname', default=None)
    parser.add_argument('--cut', help='The cutoff distance around a ligand (Å)',
                        type=float, default=6.0)
    parser.add_argument('--site', help='The manual input of the binding site. e.g. ["123", "45"]',
                        type=List[str], default=[])
    parser.add_argument('--outname', help='The name of the output file (json or yaml)',
                        type=str, default='binding_site.json')
    parser.add_argument('--outputdir', type=str, default='output',
                        help='The directory for the output json/yaml file')
    args = parser.parse_args()
    return args


def save_binding_site(residues_list: list, outname: str) -> None:
    """
    Write a binding‐site definition from a list of residue ID strings.

    Parameters
    ----------
    residues_list : list 
        Each string must identify a residue, e.g. "123", "45", etc.
    outname : str
        Name of the output file. Extension must be .json, .yaml, or .yml.
    format : str
        The format of the output file. It must be .json, .yaml, or .yml.

    Raises
    ------
    ValueError
        If the format isn’t one of .json/.yaml/.yml.
    """
    outformat = outname.rsplit('.', 1)[1]
    data = {'binding_site': np.array(residues_list, dtype=np.int64).tolist()}
    print(data)
    with open(outname, 'w') as f:
        if outformat == 'json':
            json.dump(data, f, indent=4)
        else:
            yaml.dump(data, f)


def main() -> int:
    args = read_arguments()
    outformat = args.outname.rsplit('.', 1)[1]
    if outformat not in ['json', 'yaml', 'yml']:
        raise ValueError("Output format must be json, yaml, or yml")
    # If a list of residues are provided, just convert them to the yaml file.
    # PDB file will be skipped in this case
    if args.site:
        print(f'A set of binding site residues are provided. Will convert the list to a {outformat} file.')
        save_binding_site(args.site, f'{args.outputdir}/{args.outname}')
        return 0
    # If neither the binding site residues nor the PDB file is provided, raise an error.
    elif args.site == [] and (args.pdb is None):
        raise ValueError('Please provide a PDB file to identify the binding site or manually provide the binding site.')
    # If only the PDB file is provided but the ligand name is not, raise an error.
    elif args.site == [] and (args.pdb is not None) and (args.lig is None):
        raise ValueError('A PDB file is provided but te residue name of the ligand is not specified.')
    # If both the PDB file and the ligand name are provided, find the binding site from the structure.
    elif args.site == [] and (args.pdb is not None) and (args.lig is not None):
        print('A PDB file and the resname of the ligand are provided. Will find the binding site from the structure')
        u = mda.Universe(args.pdb)
        binding_site_group = u.select_atoms(f'(nucleic or protein) and around {args.cut} (resname {args.lig})', periodic=False)
        binding_site_residues = binding_site_group.residues
        site = [res.resid for res in binding_site_residues]
        save_binding_site(site, f'{args.outputdir}/{args.outname}')
        return 0


if __name__ == '__main__':
    sys.exit(main())
