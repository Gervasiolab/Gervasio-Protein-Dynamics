#!/bin/python

import numpy as np
from time import gmtime, strftime, localtime
import argparse as ap
import textwrap
from collections import defaultdict

parser = ap.ArgumentParser(
    add_help=False,
    formatter_class=ap.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
        The script creates GROMACS topology files for running simulations using
        the SWISH method (DOI:10.1021/jacs.6b05425).
        It modifies Lennard-Jones pair potentials between water oxygen and
        apolar (C and S) atoms by linear (default) or custom (using -scales
        argument) scaling in replicas.
        
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

        The script is based on Vladas Oleinikovas' script for preparing GROMACS
        topology files for SWISH.

        Havva Yalinca modified it (21.02.2019) to:
          (i) adapt it for topology files which already contain
              nonbonded_param line(s),
         (ii) add argument parsing (based on Vladas' SWISH script),
        (iii) make it compatible with python 3 as well as python 2.7,
         (iv) check that the input was preprocessed by grompp.

        Antonija Kuzmanic modified the script further (2.10.2019) to:
          (i) bring it up to PEP8 standards,
         (ii) allow the users to choose whether they want to scale the ligand
              atoms and how,
        (iii) make it force-field agnostic (as much as possible),
         (iv) allow the users to scale only preferred residues.

        Ladislav Hovan made some tweaks (21.01.2020) to:
          (i) add the option of automatically including nonpolar hydrogens
              attached to aliphatic atoms (AMBER only),
         (ii) request confirmation after atom types are identified unless
              explicitly turned off,
        (iii) convert of boolean arguments to flags (store_true),
         (iv) fix minor bugs (now runs without -preferred).
        
        Alberto Borsatto (05.06.2021):
          (i) fix minor bugs
        
        Francesco Luigi Gervasio:
          (i) initial idea design and general supervision

        Running the script on the command line:
        python SWISH_GMX.py -f <input>.top -scale 1.0 1.1 1.2 1.3 1.4 1.5
        python SWISH_GMX.py -f <input>.top -smin 1.0 -smax 1.5 -nreps 6
        '''),
    epilog=textwrap.dedent('''\
        WARNING:

        Make sure the atom types match your input topology and the desired
        outcome. Ligand keywords -all and -heavy rely on the fact that typical
        GAFF names are lowercase. If this is not the case, please provide a
        full atom type list.
        '''))

# Define required arguments.
required_args = parser.add_argument_group('required arguments')
required_args.add_argument('-f', type=str, required=True,
                           help=('input Gromacs topology file'))

# Define optional arguments.
optional_args = parser.add_argument_group('optional arguments')
optional_args.add_argument('-h', '--help', action='help',
                           help='show this help message and exit')
optional_args.add_argument('-o', type=str, default=None,
                           help=('output Gromacs topology file name ' +
                                 '(default: <input><X>.top)'))
optional_args.add_argument('-y', action='store_true',
                           help=('skip manual confirmation after atom ' +
                                 'type detection'))

# Define protein arguments.
protein_args = parser.add_argument_group(
    'protein arguments',
    textwrap.dedent('''\
    arguments to control the scaling of protein atoms
    If the atom types are not defined, they will be determined from
    the input topology.
    '''))
protein_args.add_argument('-nreps', type=int, default=6,
                          help='number of replicas (default: %(default)s)')
protein_args.add_argument('-smin', type=float, default=1.0,
                          help='min scaling factor (default: %(default)s)')
protein_args.add_argument('-smax', type=float, default=1.5,
                          help='max scaling factor (default: %(default)s)')
protein_args.add_argument('-scale', type=float, nargs='+', default=[],
                          help=('list of scaling factors; overrides -nreps,' +
                                ' -smax, -smin'))
# Define atom types
protein_args.add_argument('-prot_atms', type=str, nargs='+', default=[],
                          help='list of protein atoms for scaling')
protein_args.add_argument('-carbonyl', type=str, default=None,
                          help='carbonyl atom type')
protein_args.add_argument('-water', type=str, default=None,
                          help='water oxygen atom type')
protein_args.add_argument('-prot_h_amber', action='store_true',
                          help='include nonpolar protein hydrogens from AMBER')
# Preferred residues
protein_args.add_argument('-preferred', action='store_true',
                          help=('scale only the atoms in preferred residues ' +
                                '(default: %(default)s)'))
protein_args.add_argument(
    '-pref_res', type=str, nargs='+', default=['TYR', 'PHE', 'TRP', 'MET',
                                               'HIS', 'ILE', 'LEU', 'GLY',
                                               'CYS', 'VAL'],
    help=('list of preferred residues for scaling (default: %(default)s)'))

# Ligands
ligand_args = parser.add_argument_group(
    'ligand arguments',
    'arguments to control the scaling of ligand atoms')
ligand_args.add_argument('-ligands', action='store_true',
                         help='scale ligand atoms (default: %(default)s)')
ligand_args.add_argument('-smin_lig', type=float, default=1.0,
                         help=('min scaling factor for the ligand atoms ' +
                               '(default: %(default)s)'))
ligand_args.add_argument('-smax_lig', type=float, default=1.5,
                         help=('max scaling factor for the ligand atoms ' +
                               '(default: %(default)s)'))
ligand_args.add_argument('-scale_lig', type=float, nargs='+', default=[],
                         help=('list of scaling factors for the ligand ' +
                               'atoms; overrides -nreps, -smax, -smin'))
ligand_args.add_argument('-lig_atms', type=str, nargs='+', default=['all'],
                         help=('list of ligand atoms for scaling; accepts ' +
                               'keywords all and heavy, or an explicit ' +
                               'list of atoms (default: %(default)s)'))

args = parser.parse_args()  # parse arguments

# give them manageable names
inname = args.f
outname = args.o
# protein args
nreps = args.nreps
smax = args.smax
smin = args.smin
scale = [sf for sf in args.scale]
# atom types
apolar = [atom for atom in args.prot_atms]
carbonyl_atom = args.carbonyl
waterO = args.water
use_amber_h = args.prot_h_amber
# preferred residues
preferred = args.preferred
preferred_res = [res for res in args.pref_res]
preferred_res.sort()
# ligand args
ligands = args.ligands
smax_lig = args.smax_lig
smin_lig = args.smin_lig
scale_lig = [sf for sf in args.scale_lig]
lig_atms = args.lig_atms

# Print the swish sign (mostly for show).
swish_sign = 'S W I S H'
row_len = 80
swish_pos = int(row_len/2) - len(swish_sign)
print('{}\n{}{}\n{}\n'.format('#'*row_len, ' '*swish_pos, swish_sign,
                              '#'*row_len))

amber_np_h = ['H1', 'H4', 'H5', 'HA', 'HC', 'HP', 'HS']
# assign the scaling factors
if len(scale) == 0:
    sfactors = list(np.around(np.linspace(smin, smax, nreps), decimals=3))
else:
    sfactors = scale
    nreps = len(scale)
print('Using the following scaling factors for the protein:\n{}\n'
      .format(', '.join([str(sf) for sf in sfactors])))
# check if the protein has a replica w/o scaling
message = 'No replica with the scaling factor 1.0 for protein atoms!'
assert 1 in sfactors, message

# check if the number of replicas is even
assert nreps % 2 == 0, 'Number of replicas should be even!'

if ligands:
    if len(scale_lig) == 0:
        sfactors_lig = list(np.around(np.linspace(smin_lig, smax_lig, nreps),
                                      decimals=3))
    else:
        sfactors_lig = scale_lig
        nreps = len(scale_lig)
        # make sure the number of replicas is the same
        message = ('Number of scaling factors for the protein and the ligand' +
                   ' atoms has to be the same!')
        assert len(sfactors) == len(scale_lig), message
    print('Using the following scaling factors for the ligands:\n{}\n'
          .format(', '.join([str(sf) for sf in sfactors_lig])))
    # check if the protein has a replica w/o scaling
    message = 'No replica with scaling factor 1.0 for the ligand atoms!'
    assert 1 in sfactors_lig, message
    # check if the replica w/o scaling matches for the protein and ligands
    message = ('There is no replica where the scaling factor is 1.0 for both' +
               ' protein and ligand atoms!')
    assert sfactors.index(1) == sfactors_lig.index(1), message


def read_topology(inname):
    """
    Read the Gromacs topology and checks if it has any include statements.
    """
    with open(inname, 'r') as topol:
        top = [line for line in topol.readlines()]
    processed = [line for line in top if line.startswith('#include ')]
    # Make sure there are no include statements in the topology file.
    message = ('The topology file contains include statements.\n' +
               'It needs to be processed using gmx grompp first.')
    assert len(processed) == 0, message
    return top


def get_atomtypes(top, carbonyl_atom):
    """
    Gets all the apolar C/S atom types that will be scaled.
    """
    # First read all the atom types for each residue.
    read_res = False
    res_atomtypes = {}
    for line in top:
        if line.startswith('\n') or line.startswith('['):
            read_res = False
        if line.startswith('; residue'):
            if line.split()[3] not in res_atomtypes:
                res_atomtypes[line.split()[3]] = []
            read_res = True
            continue
        if read_res:
            tmp = line.split()
            res = tmp[3]  # get residue name
            atype = tmp[1]  # get atom type
            mass = int(tmp[7].split('.')[0])  # get mass
            if atype.islower():
                read_res = False
                del res_atomtypes[res]
            else:
                if (atype not in res_atomtypes[res] and
                        (mass == 12 or mass == 32 or (mass == 1 and
                         use_amber_h and atype in amber_np_h))):
                    res_atomtypes[res].append(atype)

    # Remove entries for ASN, ASP, GLN, GLU (and their protonated parts)
    # to exclude their polar carbonyl carbons.
    excl_keys = [key for key in res_atomtypes.keys()
                 if key[:2] in ['AS', 'GL'] and key != 'GLY']
    for key in excl_keys:
        del res_atomtypes[key]

    if preferred:
        # Remove all but the preferred residues
        excl_keys = list(set(res_atomtypes.keys())
                         .difference(set(preferred_res)))
        for key in excl_keys:
            del res_atomtypes[key]

    # Figure out which atom is the backbone carbonyl carbon
    if carbonyl_atom is None:
        # Find the common atom in all residues.
        common = set()
        for k, v in res_atomtypes.items():
            if len(common) == 0:
                common = set(v)
            else:
                common = common.intersection(set(v))
        if len(common) == 1:
            carbonyl_atom = list(common)[0]
        elif 'C' in common:
            carbonyl_atom = 'C'
        else:
            message = ('Backbone carbonyl atom cannot be determined.\n' +
                       'Please restart the script and define it using the ' +
                       '-carbonyl flag.')
            raise ValueError(message)
        print('{} determined as the backbone carbonyl atom type.'
              .format(carbonyl_atom))
    # Get unique atom types
    atom_types = []
    for v in res_atomtypes.values():
        atom_types.extend(v)
    atom_types = list(set(atom_types))
    atom_types.remove(carbonyl_atom)
    atom_types.sort()
    return atom_types


def get_waterO(top):
    """
    Gets the water oxygen atom type.
    """
    water = False
    read_res = False
    atype = None
    for line in top:
        if line.startswith('; Include water'):
            water = True
        if water and line.startswith('[ atoms ]'):
            read_res = True
        if (line.startswith('\n') or line.startswith('[') or
                line.startswith(';')):
            continue
        if water and read_res:
            tmp = line.split()
            mass = int(tmp[7].split('.')[0])
            if mass == 16:
                atype = tmp[1]
                print('{} determined as the water oxygen atom type.'
                      .format(atype))
                return atype
    message = ('Water oxygen atom cannot be determined.\n' +
               'Please restart the script and define it using the ' +
               '-water flag.')
    assert atype is not None, message


def comb_nb2(s_i, e_i, s_j, e_j):
    """
    Calculates the LJ interaction using the combination rule 2 for sigma and
    epsilon.
    """
    s_ij = (s_i+s_j)/2.0
    e_ij = (e_i*e_j)**0.5
    return s_ij, e_ij


def format_line(natms, atoms, line):
    """
    Formats a line from one of the types sections (bonds, angles, dihedrals)
    if it appears in the protein atom selection for scaling.
    """
    changed = False
    outline = None
    tmp = line.split()
    line_atoms = []
    for atom in tmp[:natms]:
        if atom in atoms:
            atom = f'X{atom}'
            changed = True
        line_atoms.append('{:4s}'.format(atom))
    if changed:
        outline = (' {}{} ; Added for SWISH.\n'
                   .format(''.join(line_atoms), line[natms*4+1:-1]))
    return outline


def write_scaled_aOW_parm_file(inname, outname, topology, atoms, scale, rep,
                               atoms_lig=None, scale_lig=None):
    """
    Scales apolar atom interactions with water (only for preferred residues)
    and writes out Gromacs topologies.
    """
    d = defaultdict(str)  # dict for atom types parameters
    if not outname:
        outname = '{}{}.top'.format(
            '.'.join([el for el in inname.split('.') if el != 'top']), rep)
    else:
        outname = '{}{}.top'.format(
            '.'.join([el for el in outname.split('.') if el != 'top']), rep)
    outfile = open(outname, 'w')
    current_time = strftime('%Y-%m-%d %H:%M:%S', localtime())
    outfile.write(textwrap.dedent('''\
        ;
        ;\tFile \'{}\' was modified from
        ;\tdefault \'{}\'
        ;\tusing an adaptation (by HY & AK) of VO's Python script
        ;\tTime: {}
        ;\tChanges comprise added new atom types (present only in preferred
        ;\tresidues) and all associated parameters, as well as
        ;\t[ nonbond_parms ] for apolar +/- lig vs OW.
        ;\tScaling: {:.3f}
        '''.format(outname, inname, current_time, scale)))

    atomtypes = False  # flag for atomtypes section
    nonbondPresent = False  # flag for nonbonded param section
    scaled = False  # flag to indicate new nonbond param section was made
    halt_print = False  # flag to stop printing lines (for dihedrals)
    if preferred:
        # define how many atoms need to be checked in each types section
        flag_dict = {'bondtypes': 2, 'angletypes': 3, 'dihedraltypes': 4}
        flag = None  # flag for one of the types section above
        constraints = False  # flag for constrainttypes section
        change_res = False  # flag for protein residues
        lines = []  # for blocks of dihedrals with old atom types
        outlines = []  # for blocks of dihedrals with new atom types

    for ind, line in enumerate(topology):
        # calculates the scaling and prints out the nonbonded section
        # if atomtypes flag is True
        if atomtypes:
            if line.startswith('\n') or line.startswith('['):
                if len(d) > 0:
                    outfile.write(textwrap.dedent('''
                        [ nonbond_params ]
                        ; type_i type_j  f.type   sigma   epsilon
                        ; f.type=1 means LJ (not Buckingham)
                        ; sigma&eps since mixing-rule = 2
                        ; rescale apolar {} LJ interactions by {:.3f}
                        '''.format(waterO, d[list(d.keys())[0]][2])))
                    scaled = True
                    # always sort the dictionary for easier comparison
                    for atom in sorted(d.keys()):
                        parms = d[atom]
                        s_i, e_i, scale_f = parms[0], parms[1], parms[2]
                        s_iOW, e_iOW = comb_nb2(s_i, e_i, s_OW, e_OW)
                        # scale only the epsilon
                        s_iOW, e_iOW = s_iOW, e_iOW * scale_f
                        outfile.write('{:10s}   {:10s}  1   {:8.6e}   '
                                      '{:8.6e} ; Added for SWISH.\n'
                                      .format(atom, waterO, s_iOW, e_iOW))
                    d = defaultdict(str)  # reset the atom dictionary
                atomtypes = False  # reset the keyword
                if line.startswith('\n') and scaled:
                    # removes the annoying empty line before the already
                    # existing nonbonded parameters
                    continue
            # this is the line with LJ data
            if len(line) > 1 and not line.startswith(';'):
                tmp = line.split()
                if tmp[0] in atoms:  # take protein atom parameters
                    if preferred:
                        d[f'X{tmp[0]}'] = (float(tmp[5]), float(tmp[6]), scale)
                        outline = ('X{}{} ; Added for SWISH.\n'
                                   .format(tmp[0], line[len(tmp[0])+1:-1]))
                        line = line + outline
                    else:
                        d[tmp[0]] = float(tmp[5]), float(tmp[6]), scale
                if tmp[0] == waterO:  # take water parameters
                    s_OW, e_OW = float(tmp[5]), float(tmp[6])
                if ligands:
                    if 'all' in atoms_lig:
                        # if 'all', check if lowercase and take all
                        if tmp[0] == tmp[0].lower():
                            d[tmp[0]] = (float(tmp[5]), float(tmp[6]),
                                         scale_lig)
                    elif 'heavy' in atoms_lig:
                        # if 'heavy', check if lowercase and doesn't start
                        # with h (so no hydrogens)
                        if (tmp[0] == tmp[0].lower() and
                                not tmp[0].startswith('h')):
                            d[tmp[0]] = (float(tmp[5]), float(tmp[6]),
                                         scale_lig)
                    else:
                        # if given a list of atoms, take those in the list
                        if tmp[0] in atoms_lig:
                            d[tmp[0]] = (float(tmp[5]), float(tmp[6]),
                                         scale_lig)
        elif line.startswith('[ atomtypes ]'):
            atomtypes = True
            scaled = False  # restart the flag

        if preferred:
            # Change the constraints section.
            if constraints:
                # constrainttypes format a bit different from other types
                # sections (e.g. allows empty lines)
                if line.startswith('['):
                    constraints = False
                if len(line) > 1 and not line.startswith(';'):
                    changed = False
                    tmp = line.split()
                    line_atoms = []
                    for atom in tmp[:2]:
                        if atom in atoms:
                            atom = f'X{atom}'
                            changed = True
                        line_atoms.append('{:6s}'.format(atom))
                    if changed:
                        outline = (' {} {} ; Added for SWISH.\n'
                                   .format(' '.join(line_atoms), line[15:-1]))
                        line = line + outline
            elif line.startswith('[ constrainttypes ]'):
                constraints = True

            # Change one of the types sections (defined in flag_dict).
            if flag is not None:
                if line.startswith('\n') or line.startswith('['):
                    flag = None
                if len(line) > 1 and not line.startswith(';'):
                    outline = format_line(flag_dict[flag], atoms, line)
                    if outline is not None:
                        # Dihedrals with function 9 have to be written out in
                        # blocks (Gromacs complains a lot otherwise). Both 4
                        # and 9 will be written out in blocks.
                        if flag == 'dihedraltypes':
                            if line[:18] == topology[ind+1][:18]:
                                halt_print = True
                                lines.append(line)
                                outlines.append(outline)
                            else:
                                lines.append(line)
                                outlines.append(outline)
                                line = ''.join(lines) + ''.join(outlines)
                                lines = []
                                outlines = []
                                halt_print = False
                        else:
                            line = line + outline
            elif line.startswith('[') and 'types' in line:
                section = line[2:-3]
                if section in flag_dict:
                    flag = section

            # Change the atom types in the [ atoms ] entry.
            if change_res:
                if line.startswith('\n') or line.startswith('['):
                    change_res = False
                if len(line) > 1 and not line.startswith(';'):
                    tmp = line.split()
                    if tmp[3] in preferred_res and tmp[1] in atoms:
                        line = ('{} {:>3s} {} ; Changed for SWISH.\n'
                                .format(line[:13], f'X{tmp[1]}', line[18:-1]))
            elif line.startswith('; residue'):
                change_res = True

        # Takes care of the case when the topology already includes the
        # [ nonbond_params ] section.
        if nonbondPresent:
            if line.startswith(';'):
                continue
            # end of the original nonbond_params section
            if line.startswith('\n') or line.startswith('['):
                nonbondPresent = False
            else:
                ATOM1, ATOM2, _, SIGMA, EPSILON = line.split()[:5]
                outfile.write('{:10s}   {:10s}  1   {:8.6e}   {:8.6e}'
                              ' ; In the original topology \n'
                              .format(ATOM1, ATOM2, float(SIGMA),
                                      float(EPSILON)))
                # adding a comment to the end of the line to indicate it
                # was in the original topology
                continue
        elif line.startswith('[ nonbond_params ]') and scaled:
            nonbondPresent = True
            continue  # do not write the nonbond_params header again

        if not halt_print:
            outfile.write(line)
    outfile.close()
    print('Finished with {}.'.format(outname))
    return


if __name__ == '__main__':
    # Read the topology.
    top = read_topology(inname)
    # Get the atom types that appear in protein residues if not user-defined.
    # Assumes that atom type names are in upper case.
    if len(apolar) == 0:
        apolar = get_atomtypes(top, carbonyl_atom)
    if waterO is None:
        waterO = get_waterO(top)
    # Warn about atom types used in scaling.
    print(textwrap.dedent(
        '''
        WARNING:
        Interactions between {} and the following protein atoms will be scaled:
        {}

        The backbone carbonyl carbon (often also the atom type of CZ atom in
        TYR) and side-chain carbonyl carbons are excluded from the scaling
        (typically, CG in ASN/ASP and CD in GLN/GLU).
        '''
        .format(waterO, ', '.join(apolar))))
    if preferred:
        print(textwrap.dedent('''\
            The scaling will be done only for preferred protein residues:
            {}

            The following atom types will be created:
            {}
            '''.format(', '.join(preferred_res),
                       ', '.join([f'X{atom}' for atom in apolar]))))
    if ligands:
        print(textwrap.dedent('''
            Interactions between {} and these ligand atoms will be scaled:
            {}'''.format(waterO, ', '.join(lig_atms))))

    # Interactive part
    if not args.y:
        while (1):
            reply = input('Are you sure you want to continue? (y/n)')
            if reply.lower() in ['y', 'yes']:
                print('Writing the topology files...')
                break
            elif reply.lower() in ['n', 'no']:
                print('Quitting')
                quit(0)
            else:
                print('Input not recognised:', reply)

    # Write the topology files.
    for rep in np.arange(nreps):
        if ligands:
            write_scaled_aOW_parm_file(inname=inname, outname=outname,
                                       topology=top, atoms=apolar,
                                       scale=sfactors[rep], rep=rep,
                                       atoms_lig=lig_atms,
                                       scale_lig=sfactors_lig[rep])
        else:
            write_scaled_aOW_parm_file(inname=inname, outname=outname,
                                       topology=top, atoms=apolar,
                                       scale=sfactors[rep], rep=rep)
