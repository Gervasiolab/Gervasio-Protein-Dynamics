#!/bin/python

'''
(c) 2019 Antonija Kuzmanic

'''

import argparse as ap
import textwrap
import numpy as np
import pandas as pd

parser = ap.ArgumentParser(
    add_help=False,
    formatter_class=ap.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
        The script creates a PLUMED contact map using representative atoms
        from Timescapes (see 10.1021/ct900229u or Timescapes manual)
        located in alpha-helices and beta-sheets. The map is meant to be
        used as a restraint during SWISH simulations.
        The map is composed of short-range distances that are away from the
        site of interest.

        The script requires the following inputs:
        - pdb or gro file with the reference structure
        - STRIDE secondary structure assignment
        - Cartesian coordinates of the pocket centre

        The output name is optional (default - cmap.dat).
        The minimal distance from the pocket centre and the maximum
        distance between atoms can also be defined.
        Residues can also be specifically included or excluded in the
        contact map through PyMOL-like syntax (e.g. 10-17+24+12).

        Copyright (C) 2022  Francesco Luigi Gervasio

        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with this program.  If not, see <https://www.gnu.org/licenses/.
        '''),
    epilog=textwrap.dedent('''\
        WARNING:

        The script will not work for very large systems where
        atom number counter resets to 1 after reaching the space
        limit of the gro file.
        '''))

# Define required arguments.
required_args = parser.add_argument_group('required arguments')
required_args.add_argument('-f', help=('reference structure file - ' +
                                       '.gro or .pdb format'),
                           dest='fstruct', required=True, type=open)
required_args.add_argument('-s', help=('secondary structure file - ' +
                                       'STRIDE format'),
                           dest='fss', required=True, type=open)
required_args.add_argument('-com', help=('x, y, z coordinates of the pocket' +
                                         ' centre (nm) (0 0 0)'),
                           dest='com', nargs=3, required=True, type=float)

# Define optional arguments.
optional_args = parser.add_argument_group('optional arguments')
optional_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
optional_args.add_argument('-o', help='output PLUMED cmap file', dest='fout',
                           default='cmap.dat',
                           type=ap.FileType('w'))
optional_args.add_argument('-vis', dest='fpml', default='show_dists.pml',
                           help='output PyMOL script to visualise distances',
                           type=ap.FileType('w'))
optional_args.add_argument('-cutoff', help='contact distance cutoff (nm)',
                           dest='cutoff', default=0.6, type=float)
optional_args.add_argument('-mindist', help='minimal distance from the ' +
                                            'pocket centre (nm)',
                           dest='mindist', default=0.6, type=float)
optional_args.add_argument('-exclude', help='residues to exclude from the' +
                                            'contact map',
                           dest='exclude', type=str)
optional_args.add_argument('-include', help='residues to include in the' +
                                            'contact map',
                           dest='include', type=str)
optional_args.add_argument('-bckbn', help='include backbone N and O atoms in' +
                                          'the contact map',
                           dest='bckbn', default=False, type=bool)

args = parser.parse_args()  # parse arguments

# Define Timescapes representative atoms.
ts_atms = {
    "ALA": "CA", "ARG": "CZ", "ARGN": "CZ", "ASP": "CG", "ASPH": "CG",
    "ASPP": "CG", "ASH": "CG", "ASN": "CG", "ASN1": "CG", "CYS": "CB",
    "CYM": "CB", "CYSH": "CB", "CYN": "CB", "CYX": "CB", "CYS1": "CB",
    "CYS2": "CB", "GLN": "CD", "GLU": "CD", "GLH": "CD", "GLUH": "CD",
    "GLUP": "CD", "GLY": "CA", "HIS": "CG", "HIP": "CG", "HIE": "CG",
    "HID": "CG", "HISP": "CG", "HISE": "CG", "HISD": "CG", "HIS1": "CG",
    "HISH": "CG", "HISB": "CG", "HISA": "CG", "HSC": "CG", "HSP": "CG",
    "HSE": "CG", "HS2": "CG", "HSD": "CG", "ILE": "CG1", "LEU": "CG",
    "LYS": "CE", "LYN": "CE", "LSN": "CE", "LYP": "CE", "LYSH": "CE",
    "MET": "SD", "PHE": "CG", "PHEU": "CG", "PRO": "CG", "SER": "CB",
    "THR": "CB", "TRP": "CE2", "TRPU": "CE2", "TYR": "CG", "TYRU": "CG",
    "VAL": "CB"
}
# Put representative atoms in a data frame format.
atm_lib = pd.DataFrame(ts_atms.items(), columns=['res_name', 'atm_name'])
hb_bckbn = ['N', 'O']


def process_reslist(res_str):
    """
    Transform PyMOL-like syntax for residue selection into a list of residue
    numbers.
    """
    residues = []
    try:
        if '+' in res_str:
            # Separate individual values.
            reslist = res_str.split('+')
        else:
            # Set a single value as a list. Otherwise, the string is split into
            # characters.
            reslist = [res_str]
        for res in reslist:
            if '-' in res:
                # Split the residues on a hyphen and create a range of values.
                temp = res.split('-')
                residues.extend(np.arange(int(temp[0]), int(temp[1])+1))
            else:
                residues.append(int(res))
        return list(set(residues))
    except ValueError:
        print('Wrong residue selection syntax. See --help.')


def process_files(struct_file, ss_file, com, mindist):
    """
    Read in the structure and secondary structure files, and return
    representative atoms of residues/atoms defined in atm_lib which have an
    alpha-helix or beta-sheet secondary structure assignment, and which are
    at least at a selected value from the given centre of the pocket.
    """
    ext = struct_file.name.split('.')[-1].lower()  # get file extension
    if ext in set(['pdb', 'gro']):
        print('Processing ' + ext + ' file.')
        atoms = struct_file.readlines()
        if ext == 'pdb':
            # make sure you only take ATOM entries
            atoms = [x[:-1] for x in atoms if x[0:4] == 'ATOM']
            n1, n2 = 6, 11  # atom number
            a1, a2 = 12, 16  # atom name
            i1, i2 = 21, 26  # residue number
            r1, r2 = 17, 21  # residue name
            x1, x2 = 30, 38  # x coordinate
            factor = 0.1  # conversion factor for A to nm
        if ext == 'gro':
            atoms = atoms[2:-1]  # skip non-atom lines
            n1, n2 = 15, 20  # atom number
            a1, a2 = 11, 15  # atom name
            r1, r2 = 5, 9   # residue name
            i1, i2 = 0, 5  # residue number
            x1, x2 = 20, 28  # x coordinate
            factor = 1  # gro file is in nm
        # Read the file into a data frame.
        df = pd.DataFrame(zip(map(lambda x: int(x[n1:n2].strip()), atoms),
                          map(lambda x: x[a1:a2].strip(), atoms),
                          map(lambda x: int(x[i1:i2].strip()), atoms),
                          map(lambda x: x[r1:r2].strip(), atoms),
                          map(lambda x: float(x[x1:x2].strip())*factor, atoms),
                          map(lambda x: float(x[x2:x2+8].strip())*factor,
                              atoms),
                          map(lambda x: float(x[x2+8:x2+16].strip())*factor,
                              atoms)),
                          columns=['atm_num', 'atm_name', 'res_num',
                                   'res_name', 'x', 'y', 'z'])
        # Use an inner join so that among all atoms, only the ones in the
        # database are selected.
        df_join = (pd.merge(df, atm_lib, on=['res_name', 'atm_name'],
                            how='inner')
                   .sort_values('atm_num')
                   .reset_index(drop=True))
        if args.bckbn:
            df_bckbn = df[(df['atm_name'].isin(hb_bckbn)) &
                          (df['res_num'].isin(df_join['res_num'].tolist()))]
            df_join = (pd.concat([df_bckbn, df_join])
                       .sort_values('atm_num')
                       .reset_index(drop=True))
        # Read in the secondary structure assignment.
        ss = ss_file.readlines()
        ss = [x for x in ss if x[0:3] == 'ASG']
        # Check if all ok with ssfile.
        message = 'No secondary structure assignment in the provided file.'
        assert len(ss) > 0, message
        ss_df = pd.DataFrame(zip(map(lambda x: int(x[11:15].strip()), ss),
                             map(lambda x: x[24:25].strip(), ss)),
                             columns=['res_num', 'ss'])
        # Join the ss assignment to the atom details.
        df_join = pd.merge(df_join, ss_df, on=['res_num'], how='left')
        # Exclude residues from the list if provided.
        if args.exclude is not None:
            excluded = process_reslist(args.exclude)
            print('Excluding residues {}.'.format(args.exclude))
            df_join = (df_join[~df_join['res_num'].isin(excluded)]
                       .reset_index(drop=True))
        # Take only residues forming alpha-helix or beta-sheet or are on the
        # include list.
        if args.include is not None:
            included = process_reslist(args.include)
            print('Including residues {}.'.format(args.include))
        else:
            included = []
        df_join = df_join[(df_join['ss'].isin(['H', 'E'])) |
                          (df_join['res_num'].isin(included))]
        # Calculate the distance from the pocket centre and filter out atoms
        # that are too close.
        df_join = (df_join
                   .assign(com_dist=lambda x: np.sqrt((x.x-com[0])**2 +
                           (x.y-com[1])**2 + (x.z-com[2])**2)))
        df_join = df_join[df_join['com_dist'] > mindist].reset_index(drop=True)
        print('Found {} atoms that are >{} nm away from the pocket centre.'
              .format(len(df_join), mindist))
        return df_join
    else:
        raise ValueError('Wrong input file format. See --help.')


def get_pairwise_distances(xyz1, xyz2):
    """
    Get fast pairwise Euclidean distances between coordinates.
    xyz1, xyz2 : numpy arrays of shapes (N, 3) and (M, 3)
    The result is a numpy matrix array (N, M).
    """
    print('Calculating pairwise distances.')
    A2 = np.broadcast_to([np.dot(r, r.T) for r in xyz1],
                         (len(xyz2), len(xyz1))).T
    AB = np.dot(xyz1, xyz2.T).astype('float')
    B2 = np.broadcast_to([np.dot(r, r.T) for r in xyz2],
                         (len(xyz1), len(xyz2)))
    # necessary to add the absolute part due to small negative values
    return np.sqrt(np.absolute(A2 - 2*AB + B2))


def process_distances(df, cutoff):
    """
    Calculate pairwise distances for representative atoms away from the pocket.
    Return only those that are within a given cutoff.
    """
    pairwise_dist = get_pairwise_distances(df.loc[:, "x":"z"].values,
                                           df.loc[:, "x":"z"].values)
    # Create the data frame with atomic pairs and their calculated distances.
    d = {'atm1': np.tile(df['atm_num'].values, len(df['atm_num'].values)),
         'atm2': np.repeat(df['atm_num'].values, len(df['atm_num'].values)),
         'dist': pairwise_dist.reshape(len(pairwise_dist)**2,)}
    df_dist = pd.DataFrame(d)
    # Remove the duplicates and zeros and take only those within a cutoff.
    df_dist = (df_dist[(df_dist['atm1'] < df_dist['atm2']) &
               (df_dist['dist'] < cutoff)]
               .sort_values(['atm1', 'atm2'])
               .reset_index(drop=True))
    print('Keeping distances shorter than {} nm.'.format(cutoff))
    # Get the additional atom info.
    df_dist = pd.merge(df_dist, df.loc[:, 'atm_num':'res_name'],
                       left_on=['atm1'], right_on=['atm_num'], how='left')
    df_dist = pd.merge(df_dist, df.loc[:, 'atm_num':'res_name'],
                       left_on=['atm2'], right_on=['atm_num'], how='left',
                       suffixes=('1', '2'))
    # Exclude pairs that are separated by less than 5 residues and get rid of
    # the original atm_num columns
    df_dist = (df_dist[abs(df_dist['res_num1'] - df_dist['res_num2']) >= 5]
               .drop(['atm_num1', 'atm_num2'], axis=1)
               .reset_index(drop=True))
    # Separate the distance data frame into two data frames based on the atom
    # type (backbone or other).
    df_bckbn = df_dist[(df_dist['atm_name1'].isin(hb_bckbn)) |
                       (df_dist['atm_name2'].isin(hb_bckbn))]
    df_nobckbn = df_dist[(~df_dist['atm_name1'].isin(hb_bckbn)) &
                         (~df_dist['atm_name2'].isin(hb_bckbn))]
    # Allow the backbone atoms to form only one contact.
    temp = df_bckbn.groupby(['res_num1'])['dist'].min().reset_index()
    df_bckbn = pd.merge(df_bckbn, temp, on=['res_num1', 'dist'], how='right')
    # Merge the data frames back into one and return it.
    df_dist = (pd.concat([df_bckbn, df_nobckbn])
               .sort_values(['atm1', 'atm2'])
               .reset_index(drop=True))
    return df_dist


def switch(r, r0, d0=0, n=6, m=12):
    """
    Calculate the rational switching function where d0 is the point around
    which the function is 1, while r0 represents an inflexion point of the
    function and for values r>r0>d0 and r<r0<d0, the function quickly decays to
    zero.
    https://www.plumed.org/doc-v2.5/user-doc/html/switchingfunction.html
    """
    if r0 == (r-d0):
        # to avoid divisions with 0
        r = r + 0.000000001
    frac = (r-d0)/r0
    s = (1-frac**n)/(1-frac**m)
    return s


def create_cmap(df_dist):
    """
    Create a contact map following Federico's approach in his JCTC paper
    (DOI: 10.1021/acs.jctc.8b00263). He used d0=0, n=6, m=12, and r0=0.35 for
    backbone contacts and r0=0.55 for the rest. Also output a PyMOL script for
    the visualisation of distances.
    """
    sticks_res = []
    args.fout.write('CONTACTMAP ...\n')
    args.fpml.write('set cartoon_side_chain_helper, off\n')
    args.fpml.write('load {}\n'.format(args.fstruct.name))
    # Loop through the pairwise distances data frame.
    for i, row in enumerate(df_dist.iterrows()):
        d = row[1]  # get the data frame row
        j = i+1
        # Set the r0 according to the atom type.
        if (d['atm_name1'] in hb_bckbn or d['atm_name2'] in hb_bckbn):
            r0 = 0.35
        else:
            r0 = 0.55
        # Calculate and write out the contact map.
        args.fout.write(
            '{:18}   {:>18}   {:>9}{}  # {:>7}-{:3} {:>7}-{:3}\n'
            .format('ATOMS{}={},{}'.format(j, d['atm1'], d['atm2']),
                    'REFERENCE{}={:.3f}'.format(j, switch(d['dist'], r0)),
                    'SWITCH{}'.format(j),
                    '={{RATIONAL R_0={} D_0=0.00 NN=6 MM=12}}'.format(r0),
                    '{}{}'.format(d['res_name1'], d['res_num1']),
                    d['atm_name1'],
                    '{}{}'.format(d['res_name2'], d['res_num2']),
                    d['atm_name2']))
        # Write out the distance part of the PyMOL visualisation script.
        args.fpml.write('distance resi {} and name {}, resi {} and '
                        'name {}\n'.format(d['res_num1'], d['atm_name1'],
                                           d['res_num2'], d['atm_name2']))
        # Write out the stick part of the PyMOL visualisation script.
        for inf in range(1, 3):
            res = d['res_num{}'.format(inf)]
            if res not in sticks_res:
                sticks_res.append(res)
    print('Creating the contact map with {} distances.'.format(j))
    if j > 200:
        print('WARNING! The number of distances is larger than 200. ' +
              'Consider reducing their number to speed up the simulations.')
    # Finishing bits of both outfiles.
    args.fout.write('LABEL=cmap\nCMDIST\n... CONTACTMAP')
    sticks_res.sort()
    args.fpml.write('show sticks, resi {} and not name h*\n'
                    .format('+'.join([str(res) for res in sticks_res])))
    print('Done!')


if __name__ == '__main__':
    # Read and process files.
    df_atoms = process_files(args.fstruct, args.fss, args.com, args.mindist)
    message = 'No atoms selected for the map. Check your cutoff values.'
    assert len(df_atoms) > 0, message
    # Get the distances.
    df_dist = process_distances(df_atoms, args.cutoff)
    message = 'No distances selected for the map. Check your cutoff values.'
    assert len(df_dist) > 0, message
    # Create the map.
    create_cmap(df_dist)
