#!/bin/env python                                                                                                                                                             
# Author: Brett Savoie (brettsavoie@gmail.com)

import sys, argparse, os, json, codecs
import numpy as np
from copy import deepcopy
from utilities.PeriodictTable import Tables
from utilities.Excpetions import TAFFIException
from utilities.transify import transify
from utilities.util import GetInchi, find_files
from utilities.writers import write_xyz, write_modelist
from utilities.parse import parse_FF_params
from model.molecule import canon_dihedral
from model.molecule import Molecule


class FragGenAllException(TAFFIException):
    pass

def main(argv):
    parser = argparse.ArgumentParser(
        description='Reads in .xyz file(s) and generates smallest atomtype-consistent fragments for parametrizing '
                    'each missing force-field mode. Model compounds are placed in individual folders named by '
                    'inchi-key, and the generated file, dependency.txt, contains the dependency graph.')

    # required (positional) arguments
    parser.add_argument('coord_files', dest='xyz',
                        help='Space delimited string holding the names of the .xyz input geometries of one or several '
                             'molecules in need of parametrization. Wildcards are supported. ')

    DEFAULT = ''
    parser.add_argument('-FF', dest='ff', default=DEFAULT,
                        help='This variable holds the filename(s) of the .db files holding existing force-field '
                             'parameters to be avoided in the fragment generation (default: {}})'.format(DEFAULT))

    DEFAULT = ''
    parser.add_argument('-xyz_avoid', dest='xyz_avoid', default=DEFAULT,
                        help='This variable holds the space-delimited filename(s) of xyz files. Any atomtypes for '
                             'which these xyz files are canonical parametrization fragments for are avoided when '
                             'generating fragments (default: {})'.format(DEFAULT))

    DEFAULT = 'frags'
    parser.add_argument('-o', dest='outputname', default=DEFAULT,
                        help='Controls the output folder name for the fragments. (default: {}})'.format(DEFAULT))

    DEFAULT = 2
    parser.add_argument('-gens', type=int, dest='gens', default=DEFAULT,
                        help='Controls the bond search depth for identifying unique atom types (default: {})'.format(DEFAULT))

    parser.add_argument('-R', dest='recursive_opt', default=False, action='store_const', const=True,
                        help='When this flag is present, filenames and wildcards will be searched for recursively.')

    DEFAULT = 0
    parser.add_argument('-q', type=int, dest='q_tot', default=0,
                        help='The total charge on the generated fragments. (default: {})'.format(DEFAULT))

    parser.add_argument('--keep_terminal', dest='keep_terminal_opt', default=False, action='store_const', const=True,
                        help='When this flag is present, terminal atoms at the edges of fragments are retained. Default is to perform a full hydrogenation consistent with the hybridization of the atomtype. (default: off)')

    parser.add_argument('--avoid_frags', dest='avoid_frags_opt', default=False, action='store_const', const=True,
                        help='When this flag is present, the xyzs are used directly for the parameterization (default: False)')

    # Make relevant inputs lowercase
    args = parser.parse_args(argv)

    sys.stdout = Logger()
    print("PROGRAM CALL: python {}\n".format(' '.join([i for i in sys.argv[1:]])))

    args.xyz = args.xyz.split()
    args.xyz_avoid = args.xyz_avoid.split()
    args_dict = vars(args)

    # Generate the model compounds
    return run(args_dict)


def run(config):
    # Find files
    config['xyz'] = find_files(config['xyz'], recursive=config['recursive_opt'])

    # Perform some consistency checks
    if False in [i.split('.')[-1] == 'xyz' for i in config['xyz']]:
        raise FragGenAllException("ERROR: Check to ensure that the input file(s) are in .xyz format.")
    if False in [os.path.isfile(i) for i in config['xyz']]:
        raise FragGenAllException("ERROR: Could not find file {}. Exiting...".format(
            next([i for i in config['xyz'] if not os.path.isfile(i)])))
    for i in config['ff']:
        if os.path.isfile(i) == False:
            raise FragGenAllException("ERROR: FF file {} couldn't not be found. Exiting...".format(i))
    for i in config['xyz_avoid']:
        if os.path.isfile(i) == False:
            raise FragGenAllException("ERROR: xyz_avoid file {} couldn't not be found. Exiting...".format(i))

    # Generate a dictionary of the various modes contained in the supplied database files
    FF_db, modes_from_FF = parse_FF_params(config['ff'])

    # Generate model comopounds.
    mc_dict, mode_instances, instance_types = model_compound_generator(config, FF_db)

    # Generate dependency lists and assign generations. These are assigned to geo_dict, so there is no return value. 
    assign_generations(mc_dict, mode_instances)

    # Map redundant model compounds, write mode files associated with each model compound,
    # and write the transified geometries to dedicated folders (based on inchi key).
    geo_dict_out = generate_final(config, mc_dict, mode_instances, instance_types)

    # Does anything actually use this? If not, then remove and simplify generate_final. 7/17/20
    return geo_dict_out


def model_compound_generator(config, FF_db):

    """
    recursively generates model compounds for all modes present in config['xyz'] but not in the supplied force-field files (FF_db). A new model compound is generated for each mode that is discovered. This potentially results in many equivalent model compounds being generated, however redundant models are later removed by the generate_final function.

    :param config: user option configuration
    :type config: dict
    :param FF_db: existing force field data bank file
    :type FF_db: str

    :returns:
        *geo_dict* -       dictionary containing the  model compounds and their properties
        *mode_instances* - a list of tuples holding the key in geo_dict of the correponding compound and the atom
                        indices in the model compound for the mode (e.g., (5,(1,0,2)) would be an angle
                        corresponding to atoms 1, 0, and 2 in model compound geo_dict[5]. Later, assign_generations
                        inserts the generation for parameterization into the second element.\
        *instance_types* - a list of mode types, indexed to mode_instances. Each mode is specified as a tuple of the
                        atomtypes in the mode. Dihedrals have an extra term corresponding to the potential type.
    """
    # Generate a dictionary of the various atomtypes whose canonical fragments match the supplied xyz_avoid files
    xyz_avoid_set = parse_avoided_types(config["xyz_avoid"], config["gens"])

    # Initialize dictionary and lists for holding modes and geometry information

    bonds = []
    bond_types = []
    bond_keys = []
    angles = []
    angle_types = []
    angle_keys = []
    dihedrals = []
    dihedral_types = []
    dihedral_keys = []

    # Loop over the input structures, generate the all trans-minimized conformer,  and parse the modes
    print("\n" + "*" * 120)
    print("* {:^116s} *".format("Parsing atomtypes and modes in need of parameterization for the following files"))
    print("*" * 120 + "\n")
    for i in config["xyz"]: print("  {}".format(i))

    geo_dict = {}
    # Initialize geo_dict with the supplied coord_file information
    for count_i, i in enumerate(config["xyz"]):
        geo_dict[count_i] = Molecule(config['gens'])
        geo_dict[count_i].get_atom_types_pipe(i, method='complicate')

    # Find atom types in need of parameterization
    # TODO Another headache of using depth 2 angles and dihedrals is that we have to generate model compounds for the atom types first
    #     then generate the model compounds for the rest of the modes. 
    geo_dict_0 = deepcopy(geo_dict)
    # the followings are being updated in update_mc_dict
    geo_keys = list(geo_dict.keys())
    mc_dict = {}
    atoms = []
    atom_types = []
    atom_keys = []
    for i in geo_keys:

        # Find first instance of each atom and mode in need of parameterization
        # Generate their model compounds and propagate the information to the atom/bond/angle/dihedral lists
        # mc_dict holds the model compound information and gets updated by the helper function. 
        # For atoms
        for count_j, j in enumerate(geo_dict[i]["atom_types"]):
            if j not in atom_types:
                # Generate model
                update_mc_dict([count_j], j, geo_dict, mc_dict, config, geo_keys, i, atoms, atom_types, atom_keys)

        del geo_dict[i]

    # Determine atom type dependencies between model compounds
    keys = mc_dict.keys()
    for i in keys:
        tmp = []
        for j in keys:
            if any([k in mc_dict[j].mc_prop["min_types"] for k in mc_dict[i].mc_prop["dep_types"]]):
                tmp += [j]
        mc_dict[i]["dep_keys"] = tmp

    # TODO a duplicate section in assign_generation, will need to refactor this
    # Assign the generation to the model compounds for atom types (at_mc)
    # Description: Algorithm starts all calculations out at the same generation and moves individual 
    # fragments to later generations based on their dependences. Break occurs once no generation update occur
    for i in keys: mc_dict[i].gens = 1
    flag = True
    count = 0
    while flag is True:
        flag = False
        count += 1

        # Loop over each fragment and if it is in the same generation as one of its dependencies, increment its generation upward
        for i in keys:
            tmp_max = [mc_dict[j].gens for j in mc_dict[i].mc_prop["dep_keys"]]
            if len(tmp_max) > 0:
                tmp_max = max(tmp_max)
                if mc_dict[i].gens <= tmp_max:
                    mc_dict[i].gens = tmp_max + 1
                    flag = True

    # 1) Restore geo_dict to its original input structures,
    # 2) add a "gen" element that comes after the atom type model compounds, 
    # 3) add atomtype model compounds (current mc_dict) to restored geo_dict. 
    max_gen = max([mc_dict[i].gens for i in mc_dict.keys()])
    geo_dict = geo_dict_0  # restore original dictionary of input structures.
    for i in geo_dict.keys(): geo_dict[i].gens = max_gen + 1  # add gen argument
    for i in mc_dict.keys():  # add atomtype model compounds to geo_dict
        geo_dict[len(geo_dict.keys())] = mc_dict[i]

    # Find atoms and modes in need of parameterization starting from lowest generation up
    # This ordering is required to avoid angles and dihedrals being parameterized to parent-dependent atoms at positions more than two bonds away.
    # Model compounds for parameterization are added to geo_dict at each iteration
    # The loop only terminates when no new model compounds are generated
    atoms = []
    atom_types = []
    atom_keys = []
    _, geo_keys = map(list, zip(*sorted([(geo_dict[i].gens, i) for i in geo_dict.keys()])))
    mc_dict = {}  # Holds the actual model compounds information
    for i in geo_keys:

        # Find modes (returns all modes in canonicalized format)
        (bonds_tmp, angles_tmp, dihedrals_tmp, one_fives_tmp, bond_types_tmp, angle_types_tmp, dihedral_types_tmp,
         one_five_types_tmp) = geo_dict[i].find_modes(return_all=1)

        # Find first instance of each atom and mode in need of parameterization
        # Generate their model compounds and propagate the information to the atom/bond/angle/dihedral lists
        # mc_dict holds the model compound information and gets updated by the helper function. 
        # For atoms
        for count_j, j in enumerate(geo_dict[i].atom_types):
            if j not in atom_types and j not in xyz_avoid_set and (
                    (j, j) not in FF_db["vdw"].keys() or j not in FF_db["charges"].keys()):
                # Generate model
                update_mc_dict([count_j], j, geo_dict, mc_dict, config, geo_keys, i, atoms, atom_types, atom_keys)

        # For bonds
        for count_j, j in enumerate(bond_types_tmp):
            if j not in bond_types and j not in FF_db["bonds"]:
                # Generate model
                update_mc_dict(bonds_tmp[count_j], j, geo_dict, mc_dict, config, geo_keys, i, bonds, bond_types,
                               bond_keys)

        # For angles
        for count_j, j in enumerate(angle_types_tmp):
            if j not in angle_types and j not in FF_db["angles"]:
                # Generate model
                update_mc_dict(angles_tmp[count_j], j, geo_dict, mc_dict, config, geo_keys, i, angles, angle_types,
                               angle_keys)

        # For dihedrals
        for count_j, j in enumerate(dihedral_types_tmp):
            if j not in dihedral_types and j not in FF_db["dihedrals"]:
                # Generate model
                update_mc_dict(dihedrals_tmp[count_j], j, geo_dict, mc_dict, config, geo_keys, i, dihedrals,
                               dihedral_types, dihedral_keys)

        # Remove the current compound after all of its model compounds have been generated and added to mc_dict
        del geo_dict[i]

    # Define mass dictionary for determining which dihedrals are explicitly scanned
    mass_dict = Tables.MASSDICT

    # Collect mode_instances and instance_types
    mode_instances = []  # list of ( geo_key, mode ) tuples, where mode is a tuple of the atom indices of the mode
    instance_types = []  # list of  mode_types indexed to mode_instances
    print("\n" + "*" * 120)
    print("* {:^116s} *".format("Force-field terms in need of parameterization"))
    print("*" * 120 + "\n")
    if len(atom_types) > 0:
        print("    Atomtypes ({}):\n\n\t".format(len(set(atom_types))) + "\n\t".join(atom_types) + "\n")
        mode_instances += [(atom_keys[count_i], i) for count_i, i in enumerate(atoms)]
        instance_types += [(i,) for i in atom_types]
    if len(bond_types) > 0:
        print("    Bonds ({}):\n\n\t".format(len(set(bond_types))) + "\n\t".join(
            ["{:<40s} {:<40s}".format(*_) for _ in bond_types]) + "\n")
        mode_instances += [(bond_keys[count_i], i) for count_i, i in enumerate(bonds)]
        instance_types += bond_types
    if len(angle_types) > 0:
        print("    Angles ({}):\n\n\t".format(len(set(angle_types))) + "\n\t".join(
            ["{:<40s} {:<40s} {:<40s}".format(*_) for _ in angle_types]) + "\n")
        mode_instances += [(angle_keys[count_i], i) for count_i, i in enumerate(angles)]
        instance_types += angle_types
    if len(dihedral_types) > 0:
        print("    Dihedrals ({}) (*=implicit):\n".format(len(set(dihedral_types))))

        implicit_dihedrals = []
        assigned_dihedrals = []
        for count_i, i in enumerate(dihedral_types):

            # Skip the dihedral if it has already been assigned
            if i in assigned_dihedrals: continue

            # Harmonic dihedrals (TAFFI treats 2-3 atoms connected by a double bond explicitly as harmonic functions. 
            if i[4] == "harmonic":
                assigned_dihedrals += [i]
                mode_instances += [(dihedral_keys[count_i], dihedrals[count_i])]
                instance_types += [i]
                print("\t{:<40s} {:<40s} {:<40s} {:<40s} {:<10s}".format(*i))
                continue

            # Gather coincident dihedrals (dihedrals that share the 2-3 bond with this dihedral)
            left_atoms = [count_k for count_k, k in
                          enumerate(mc_dict[dihedral_keys[count_i]].adj_mat[dihedrals[count_i][1]]) if
                          k == 1 and count_k != dihedrals[count_i][2]]
            right_atoms = [count_k for count_k, k in
                           enumerate(mc_dict[dihedral_keys[count_i]].adj_mat[dihedrals[count_i][2]]) if
                           k == 1 and count_k != dihedrals[count_i][1]]
            coincident_dihedrals = [tuple([k]) + dihedrals[count_i][1:3] + tuple([m]) for k in left_atoms for m in
                                    right_atoms]
            coincident_types = [tuple([mc_dict[dihedral_keys[count_i]].atom_types[m] for m in k] + ["opls"]) for k in
                                coincident_dihedrals]

            # Canonicalize the dihedrals 
            for count_k, k in enumerate(coincident_dihedrals):
                coincident_types[count_k], coincident_dihedrals[count_k] = canon_dihedral(coincident_types[count_k],
                                                                                          coincident_dihedrals[count_k])

            # Find the dihedral that will be explicitly scanned (the dihedral with the largest combined mass on the 1 and 4 atoms gets explicitly scanned)
            max_mass = 0.0
            for count_k, k in enumerate(coincident_dihedrals):
                current_mass = mass_dict[mc_dict[dihedral_keys[count_i]]["elements"][k[0]]] + mass_dict[
                    mc_dict[dihedral_keys[count_i]]["elements"][k[3]]]
                if current_mass > max_mass:
                    max_mass = current_mass
                    keep = count_k

            # Add the implicitly treated dihedrals to the implicit_dihedrals list and update the mode_instances and instance_types lists with the explicit dihedral
            implicit_dihedrals += [k for k in coincident_types if
                                   (k != coincident_types[keep] and k != coincident_types[keep][::-1])]
            mode_instances += [(dihedral_keys[count_i], coincident_dihedrals[keep])]
            instance_types += [coincident_types[keep]]
            assigned_dihedrals += coincident_types
            print("\t{:<40s} {:<40s} {:<40s} {:<40s} {:<10s}".format(*coincident_types[keep]))

        # Print the dihedrals that are being implicitly scanned
        if implicit_dihedrals:
            for i in set(implicit_dihedrals):
                print("\t{:<40s} {:<40s} {:<40s} {:<40s} {:<10s} *".format(*i))
    print("\n\tTotal modes being explicitly parametrized: {}".format(len(mode_instances)))

    return mc_dict, mode_instances, instance_types


# TODO If this breaks, we may need to Revise in future to sort dependence based on atoms, then bonds, then angles, then dihedrals.
def assign_generations(mc_dict, mode_instances):
    """
    # Description: determines the dependencies among the model compounds based on the atom types.
    #              "dep_keys" holds a list of model compound indices that each model compound depends on
    #              "gen" holds the generation that the model compound belongs to. The generation info
    #              also added to each mode_instance as the second element of the tuple.
    # Returns:     None, data is added to the supplied mc_dict dictionary.
    """

    # Find the model compounds corresponding to atom types
    at_mc = [i[0] for i in mode_instances if len(i[1]) == 1]
    others = [_ for _ in mc_dict.keys() if _ not in at_mc]

    # Determine atom type dependencies between model compounds
    keys = mc_dict.keys()
    for i in keys:
        tmp = []
        for j in keys:
            if any([k in mc_dict[j].mc_prop["min_types"] for k in mc_dict[i].mc_prop["dep_types"]]):
                tmp += [j]
        mc_dict[i].mc_prop["dep_keys"] = tmp

    # Assign the generation to the model compounds for atom types (at_mc)
    # Description: Algorithm starts all calculations out at the same generation and moves individual
    # fragments to later generations based on their dependences. Break occurs once no generation update occur
    for i in at_mc: mc_dict[i].gens = 1
    flag = True
    count = 0
    while flag is True:
        flag = False
        count += 1

        # Loop over each fragment and if it is in the same generation as one of its dependencies, increment its generation upward
        for i in at_mc:
            tmp_max = [mc_dict[j].gens for j in mc_dict[i].mc_prop["dep_keys"]]
            if len(tmp_max) > 0:
                tmp_max = max(tmp_max)
                if mc_dict[i].gens <= tmp_max:
                    mc_dict[i].gens = tmp_max + 1
                    flag = True

        # Check for an unreasonable number of generations
        if count >= 1000:
            raise FragGenAllException(
                "There appears to be a mutual dependence among model compounds. assign_generations failed after 1000 attempts.")

    # Assign the generation to the model compounds for modes (others)
    for i in others:
        mc_dict[i].gens = max([mc_dict[j].gens for j in mc_dict[i].mc_prop["dep_keys"]])

    # Add generation information to each mode_instance as the second element
    for count_i, i in enumerate(mode_instances):
        mode_instances[count_i] = (i[0], mc_dict[i[0]].gens, i[1])

    return


def generate_final(config, mc_dict, mode_instances, instance_types):
    """
    # Description: this function handles mapping redundant model compounds, writing the mode files, and outputting
    #              the json dictionary. Mapping is necessary because a model compound is generated for every mode
    #              by the model_compoound_generator function. Many of these model compounds are equivalent, and the
    #              map_frag function here takes care of determining equivalency and making sure that the mode indices
    #              are aligned.
    #
    # Returns:     mc_dict_out: a dictionary with the final model compound geometry atomtypes etc.
    """
    #### Instead of output a lot of folders and xyz file, now we simply save all the ouput into a dict: mc_dict_out
    #### similar to one used in frag_gen, but now mc_dict_out[file], the filename is now inchikey
    #### which contains: elements,geometry,atom_types,gen ; so that we can directly input it into frag_gen.py
    #### Main difference: has generation
    _, keys = zip(*sorted([(mc_dict[i].gens, i) for i in mc_dict.keys()]))

    mc_dict_out = {}

    # Save one instance of each unique fragment and a frag-*.modes file to keep track of what modes to parametrize for each fragment
    print("\n" + "*" * 120)
    print("* {:^116s} *".format("Beginning canonical fragment generation"))
    print("*" + "-" * 118 + "*")
    print("* {:<116s} *".format(
        "NOTE: Each atomtype and bond are parametrized using the smallest molecular fragment that is consistent with its"))
    print("* {:<116s} *".format(
        "      topological definition. For angles and dihedrals, the smallest molecular fragment is used that is consistent"))
    print("* {:<116s} *".format("      with their central atom, and bond, respectively."))
    print("*" * 120 + "\n")
    run_list = []
    N_frag = 0
    for i in keys:

        # Skip already assigned fragments
        if i in run_list: continue

        # Find equivalent fragments based on a comparison of the hash lists
        match_idx = []
        for j in keys:
            if mc_dict[i].hash_list == mc_dict[j].hash_list:

                # Add the current match, and remap the atomtypes and mode indices
                # Note: map_frag does the heavy lifting of mapping the atom indices in j to i (graph isomorphism problem!) the idea is that multiple modes use the same model compound
                # so only one is subjected to explicit calculations, but the distinct atomtypes for fitting each mode are retained in teh *.modes file.
                match_idx += [j]
                mc_dict[j]["atom_types"], mapping, status = map_frag(mc_dict[i], mc_dict[j])
                mode_instances = [(m[0], m[1], tuple([mapping[k] for k in m[2]])) if m[0] == j else m for m in
                                  mode_instances]

                # Print problem fragment(s) in the case of failure status being returned by map_frag
                if status == 1:
                    mc_dict[i].write_xyz('incomp_fixed')
                    mc_dict[j].write_xyz('incomp_mapped')
                    quit()

                # add this model compound key to run list because it is being assigned here
                run_list += [j]

        # Generate inchi key for this model compound
        opt_geo = transify(mc_dict[i].geometry, mc_dict[i].adj_mat, elements=mc_dict[i].elements,
                           opt_terminals=True, opt_final=True)
        inchikey = GetInchi(mc_dict[i].elements, opt_geo)

        # mode_inds
        mode_inds = [count_j for count_j, j in enumerate(mode_instances) if j[0] in match_idx]

        # Assemble the dictionary of modes being scanned and print diagnostics for the user
        # NOTE: dihedrals are treated differently from bonds and angles because in some cases two dihedrals have identical fragments and simply
        #       differ in the atomtypes/VDW/charges used when fitting the potential.
        print("\n\tModes associated with model compound ({:<s}) (inchi: {}) (gen in leading parentheses):".format(
            str(N_frag), inchikey))
        scanned_modes = {}
        placed = []
        for j in [k for k in mode_inds if len(mode_instances[k][2]) == 1]:
            scanned_modes[instance_types[j]] = {"gens": [mode_instances[j][1]], "modes": [mode_instances[j][2]],
                                                "atomtypes": [mc_dict[mode_instances[j][0]].atom_types]}
            print("\t\t({})    {}".format(mode_instances[j][1],
                                          " ".join(["{:40s}".format(m) for m in instance_types[j]])))
        for j in [k for k in mode_inds if len(mode_instances[k][2]) == 2]:
            scanned_modes[instance_types[j]] = {"gens": [mode_instances[j][1]], "modes": [mode_instances[j][2]],
                                                "atomtypes": [mc_dict[mode_instances[j][0]].atom_types]}
            print("\t\t({})    {}".format(mode_instances[j][1],
                                          " ".join(["{:40s}".format(m) for m in instance_types[j]])))
        for j in [k for k in mode_inds if len(mode_instances[k][2]) == 3]:
            scanned_modes[instance_types[j]] = {"gens": [mode_instances[j][1]], "modes": [mode_instances[j][2]],
                                                "atomtypes": [mc_dict[mode_instances[j][0]].atom_types]}
            print("\t\t({})    {}".format(mode_instances[j][1],
                                          " ".join(["{:40s}".format(m) for m in instance_types[j]])))
        for j in [k for k in mode_inds if len(mode_instances[k][2]) == 4]:
            print("\t\t({})    {}".format(mode_instances[j][1],
                                          " ".join(["{:40s}".format(m) for m in instance_types[j]])))
            if j in placed: continue
            scanned_modes[instance_types[j]] = {"gens": [], "modes": [], "atomtypes": []}

            # If j is a harmonic dihedral then redundancy is assessed on the basis of all four mode
            # indices matching (this is necessary beacuse of the way that harmonic dihedrals are scanned)
            if 2 in [n[mode_instances[j][2][1], mode_instances[j][2][2]] for n in
                     mc_dict[mode_instances[j][0]].bond_mat]:

                # Iterate over all dihedrals that have identical atomic indices in their definition
                for m in [n for n in mode_inds if len(mode_instances[n][2]) == 4]:

                    # Save dihedrals that have indentical atomic indices
                    if mode_instances[j][2] == mode_instances[m][2] or mode_instances[j][2] == mode_instances[m][2][
                                                                                               ::-1]:
                        placed += [m]
                        scanned_modes[instance_types[j]]["gens"] += [mode_instances[m][1]]
                        scanned_modes[instance_types[j]]["modes"] += [mode_instances[m][2]]
                        scanned_modes[instance_types[j]]["atomtypes"] += [mc_dict[mode_instances[m][0]].atom_types]

            # If j is a fourier dihedral then redundancy is assessed on the basic of the 2-3 atoms in the dihedral
            else:

                # Iterate over all dihedrals that have identical fragments
                for m in [n for n in mode_inds if len(mode_instances[n][2]) == 4]:

                    # Save dihedrals that have identical 2-3 atoms in their definition
                    if mode_instances[j][2][1:3] == mode_instances[m][2][1:3] or mode_instances[j][2][1:3] == \
                            mode_instances[m][2][1:3][::-1]:
                        placed += [m]
                        scanned_modes[instance_types[j]]["gens"] += [mode_instances[m][1]]
                        scanned_modes[instance_types[j]]["modes"] += [mode_instances[m][2]]
                        scanned_modes[instance_types[j]]["atomtypes"] += [mc_dict[mode_instances[m][0]].atom_types]

        # Find real atomtypes for this model compound (needed for subsequent parsers)
        tmp_mol = deepcopy(mc_dict[i])
        tmp_mol.update_atomtypes_extra()
        model_atomtypes = tmp_mol.atom_types
        del tmp_mol

        # If mc_dict_out is not needed for the return, then these objects should just be passed to the json dictionary and not go through this intermediate. 7/17/20
        mc_dict_out[inchikey] = {}
        mc_dict_out[inchikey]["gen"] = sorted(set([k for j in scanned_modes.keys() for k in scanned_modes[j]["gens"]]))
        mc_dict_out[inchikey]["atom_types"] = model_atomtypes
        mc_dict_out[inchikey]["q"] = mc_dict[i]["q_tot"]

        # Make inchi_folder
        if os.path.isdir("{}".format(inchikey)) is False:
            os.makedirs("{}".format(inchikey))

        # Write the frag-*.modes file associated with this fragment
        write_modelist("{}/{}.modes".format(inchikey, inchikey), scanned_modes, mc_dict[i]["bond_mat"])
        write_xyz("{}/{}".format(inchikey, inchikey), mc_dict[i]["elements"], opt_geo, additional_col=model_atomtypes)

        # Update the list of modes that have been associated with a fragment
        run_list += match_idx
        N_frag += 1

    print("\n" + "*" * 120)
    print("* {:^116s} *".format("Fragment generation complete!"))
    print("*" * 120 + "\n")

    # dump all.json for the taffi driver to utilize for running the parameterization jobs
    flag_keywords = ["geoopt", "modescan"]
    all = {_: False for _ in flag_keywords}
    all["gens"] = {i: {_: False for _ in ["initial_fit", "charges", "after_charge", "vdw", "final_fit"]} for i in
                   set(k for j in mc_dict_out.keys() for k in mc_dict_out[j]["gen"])}
    all["mc"] = {i: {"gen": mc_dict_out[i]["gen"], "atom_types": mc_dict_out[i]["atom_types"], "q": mc_dict_out[i]["q"]}
                 for i in mc_dict_out.keys()}
    json.dump(all, codecs.open('all.json', 'w', encoding='utf-8'), indent=10)

    return mc_dict_out


def parse_avoided_types(xyz_avoid, gens):
    """
    Description:   Returns the list of atomtypes whose canonical fragments are included in
                   the xyz_avoid geometries.
    Inputs:        xyz_avoid:     A list of xyz filenames
                   gens:          The bond depth for determining atom uniqueness
    Outputs:       xyz_avoid_set: A set of atomtypes whose canonical fragments are
                                   contained in xyz_avoid.
    """
    xyz_avoid_set = []

    for file in xyz_avoid:

        mol = Molecule(gens)
        mol.get_atom_types_pipe(file, method='simple')

        for types in set(mol.atom_types):

            if minimal_structure(types, mol) is True:
                xyz_avoid_set += [types]

    return set(xyz_avoid_set)


def update_mc_dict(M, M_type, geo_dict, mc_dict, config, geo_keys, current, M_list, M_types, M_keys):
    """
    Description: helper function for model_compound_generator that updates the model compound dictionary
              based on the supplied mode and the parent compound information in geo_dict. The logic is
              that the model compound is generated for each mode and added to mc_dict. This ends up
              producing a lot of redundant model compounds, but these are later removed. The algorithm is
              simpler this way rather than trying to generate the minimal set directly.
    """
    # geo_dict[current] is the parent structure

    # Generate model
    geo_id = max(geo_dict.keys()) + 1
    geo_keys += [geo_id]
    N_M, N_molecule = mode_frag(M, geo_dict[current], keep_terminal=config["keep_terminal_opt"], keep_types=True,
                                force_linear=True)
    geo_dict[geo_id] = N_molecule
    # TODO Needs to be updated to propagate the charge of ionic model compounds. maybe this needs to take place within mode_frag 7/10/20
    geo_dict[geo_id].parse_data(q_tot=geo_dict[current].q_tot)

    # Determine bonding matrix for the model compound
    geo_dict[geo_id].update_bondmat_fc()

    # Save model compound and original atom types to mc_dict
    mc_id = len(mc_dict.keys())
    mc_dict[mc_id] = deepcopy(geo_dict[geo_id])
    M_list += [N_M]
    M_types += [M_type]
    M_keys += [mc_id]

    # Determine true atomtypes for the model compound (note, we are using keep_types = True for mode_frag). geo_dict needs these for
    # correct mode assignment in subsequent iterations. And atomtype model compounds also need these in mc_dict.
    geo_dict[geo_id].update_atomtypes_extra()

    # Because angles and dihedrals use d=2 model compounds, calculated with respect to the central atoms, they are parameterized with the parent compound atom types
    # however, bond and atom type model compounds need to use the actual atomtypes corresponding to the model compound. XXX We will revisit this in the future by potentially
    # using depth=1 atom types for angle and dihedral fitting.
    if len(M) in [1, 2]:
        mc_dict[mc_id].parse_data(atom_types=deepcopy(geo_dict[geo_id].atom_types))

    # Find the atom types for which this is a model compound (empty if this is a model for an intramolecular mode) 
    if len(M) == 1:
        mc_dict[mc_id].mc_prop['min_types'] = set([j for j in mc_dict[mc_id].atom_types if
                                                   minimal_structure(j, mc_dict[mc_id]) is True])
    else:
        mc_dict[mc_id].mc_prop['min_types'] = set([])

    # Find the atom types that this mode depends on
    if len(M) == 1:
        mc_dict[mc_id].mc_prop["dep_types"] = set([j.replace("link-", "") for j in mc_dict[mc_id].atom_types if
                                                   minimal_structure(j, mc_dict[mc_id]) is False])
    else:
        mc_dict[mc_id].mc_prop["dep_types"] = set([j.replace("link-", "") for j in mc_dict[mc_id]["atom_types"]])
    return


def minimal_structure(atomtype, molecule):
    """
    Description:   Checks if the supplied geometry corresponds to the minimal structure of the molecule

    # Outputs:       Boolean:       True if geo is the minimal structure for the atomtype, False if not.
    """

    # Check if this is a ring type, if not and if there are rings
    # in this geometry then it is not a minimal structure.
    if "R" not in atomtype:
        if True in ["R" in i for i in molecule.atom_types]:
            return False

    # Check minimal conditions
    count = 0
    for count_i, i in enumerate(molecule.atom_types):

        # If the current atomtype matches the atomtype being searched for then proceed with minimal geo check
        if i == atomtype:
            count += 1

            # Initialize lists for holding indices in the structure within "gens" bonds of the seed atom (count_i)
            keep_list = [count_i]
            new_list = [count_i]

            # Carry out a "gens" bond deep search
            for j in range(molecule.gens):

                # Id atoms in the next generation
                tmp_new_list = []
                for k in new_list:
                    tmp_new_list += [count_m for count_m, m in enumerate(molecule.adj_mat[k]) if
                                     m == 1 and count_m not in keep_list]

                # Update lists
                tmp_new_list = list(set(tmp_new_list))
                if len(tmp_new_list) > 0:
                    keep_list += tmp_new_list
                new_list = tmp_new_list

            # Check for the minimal condition
            keep_list = set(keep_list)
            if False in [molecule.elements[j] == "H" for j in range(len(molecule.elements)) if j not in keep_list]:
                minimal_flag = False
            else:
                minimal_flag = True

    return minimal_flag


def mode_frag(M, molecule, fc=[], keep_lone=[],
              keep_terminal=False, keep_rings=False, keep_types=True, force_linear=False, return_FC=False):
    """
    Returns the geometry and relevant property lists corresponding to the
    smallest fragment that is consistent with the mode being parametrized.
    If the 1 or 4 atom is a ring, then the entire ring is included.
    """
    # TODO Doesn't return the charge state of the new model compound 7/10/20
    # Initialize mass_dict (used for identifying the dihedral among a coincident set that will be explicitly scanned)
    if hasattr(mode_frag, "mass_dict") is False:
        mode_frag.mass_dict = Tables.MASSDICT

    # Check the consistency of the supplied arguments
    if molecule.gens < 1: raise FragGenAllException(
        "ERROR in mode_frag: gens variable must be an integer >= 1. Exiting...")
    if keep_terminal not in [True, False]: raise FragGenAllException("ERROR in mode_frag: keep_terminal \
                                                                     argument must be set to a bool. Exiting...")
    if keep_rings not in [True, False]: raise FragGenAllException("ERROR in mode_frag: keep_rings argument \
                                                                  must be set to a bool. Exiting...")

    # Generate the fragment geometry/attribute lists corresponding to the mode
    N_M, N_molecule = new_mode_geo(M, molecule, force_linear=force_linear)

    # determine loop list and fixed_bonds
    # loops over resonance structures and only fixes bonds if consistent bonding occurs across structures.
    # loop_list contains the atoms in the mode whose bond structure is fixed (if consistent across resonance structures)
    total_fixed_bonds = []
    for lb, bond_mat in enumerate(N_molecule.bond_mat):

        fixed_bonds = []
        if len(M) == 1:
            loop_list = N_M

        elif len(M) == 2:
            loop_list = N_M
            fixed_bonds += [(N_M[0], N_M[1], int(bond_mat[N_M[0], N_M[1]]))]

        elif len(N_M) == 3:
            loop_list = [N_M[1]]
            fixed_bonds += [(N_M[0], N_M[1], int(
                bond_mat[N_M[0], N_M[1]]))]  # XXX I think these should be removed, redundant with loop_list? 7/10/20
            fixed_bonds += [(N_M[1], N_M[2], int(bond_mat[N_M[1], N_M[2]]))]

        elif len(N_M) == 4:
            if sum([N_molecule.adj_mat[N_M[0], j] for j in N_M[1:]]) == 3:
                print("WARNING: USING IMPROPER DIHEDRAL CRITERIA ON {} CHECK THE RESULT".format(
                    [N_molecule.atom_types[i] for i in N_M]))
                loop_list = N_M[0]
            else:
                loop_list = N_M[1:3]
                fixed_bonds += [(N_M[1], N_M[2], int(bond_mat[N_M[1], N_M[2]]))]
        else:
            raise FragGenAllException(
                "ERROR in mode_frag: Protocol doesn't exist for modes involving more than 4 atoms.")

            # Include the atoms in the mode and connected atoms within the preserve list.
        for i in loop_list:
            fixed_bonds += [(i, j, int(k)) for j, k in enumerate(bond_mat[i]) if k > 1]

        total_fixed_bonds += [fixed_bonds]

    # only if all fixed_bonds are same, take it as fixed_bonds; else, fixed_bonds=[]
    if len(list(map(list, set(map(tuple, total_fixed_bonds))))) > 1:
        fixed_bonds = []
    else:
        fixed_bonds = total_fixed_bonds[0]

        # Include the atoms in the mode and connected atoms within the preserve list.
    preserve = []
    gs = N_molecule.AdjMat.graph_seps()
    for i in loop_list:
        preserve += [count_i for count_i, i in enumerate(gs[i]) if i < molecule.gens]
    preserve = set(preserve)

    # Perform Hydrogenation
    N_molecule, added_idx = add_hydrogens(N_molecule, preserve=preserve, fixed_bonds=fixed_bonds)

    # Update the link types
    tmp_atom_types = deepcopy(N_molecule.atom_types)
    N_molecule.AtomType.atom_types += ['link-' + tmp_atom_types[i] for i in added_idx]

    # Canonicalize by sorting the elements based on hashing (NOTE: range(len(N_atom_types)) is used here rather than "atoms" as in the keep_terminal is True option. 
    Masses = [mode_frag.mass_dict[N_molecule.elements[i]] for i in range(len(N_molecule.atom_types))]
    hash_list, atoms = [list(j) for j in zip(*sorted(
        [(N_molecule.AtomType.atom_hash(count_i, Masses), i) for count_i, i in
         enumerate(range(len(N_molecule.atom_types)))],
        reverse=True))]

    # Update lists/arrays based on the sorted atoms
    N_M = tuple([atoms.index(i) for i in N_M])
    N_Geometry = N_molecule.geometry[atoms]
    N_Adj_mat = N_molecule.adj_mat[atoms]
    N_Adj_mat = N_Adj_mat[:, atoms]
    N_Elements = [N_molecule.elements[i] for i in atoms]
    N_Atom_types = [N_molecule.atom_types[i] for i in atoms]

    # Note: the gens here is temporary
    N_molecule = Molecule(molecule.gens)
    N_molecule.parse_data(elements=N_Elements, geometry=N_Geometry, atom_types=N_Atom_types, adj_mat=N_Adj_mat,
                          hash_list=hash_list)

    # If keep_types is False then the atomtypes are recalculated here
    if keep_types is False:
        N_molecule.AtomType.UpdateAtomTypes()

    return N_M, N_molecule


# only used in frag_gen*.py
def add_hydrogens(molecule, preserve=[], saturate=True, retype=True, fixed_bonds=[]):
    """
    Add hydrogens based upon the supplied atom types.
    This function is only compatible with TAFFI atom types
    NOTE: Hydrogenation heuristics for geometry assume carbon behavior. This isn't usually a
    problem when the results are refined with transify, but more specific rules should be implemented in the future
    """
    # Initialize the saturation dictionary the first time this function is called
    if not hasattr(add_hydrogens, "sat_dict"):
        add_hydrogens.sat_dict = Tables.SATURATION

        add_hydrogens.lone_e = Tables.LONE_ELETRON

        add_hydrogens.frag = 0

    # Make copy of members to be mutable
    # TODO: eventually we want to change the data in molecule directly, but for refactoring this is less likely to have any bugs, should also think about whether there should be a mutable molecule class
    geo = deepcopy(molecule.geometry)
    adj_mat = deepcopy(molecule.adj_mat)
    atomtypes = deepcopy(molecule.atom_types)
    elements = deepcopy(molecule.elements)
    q_tot = deepcopy(molecule.q_tot)

    # Intermediate scalars
    H_length = 1.1
    N_atoms = len(geo)
    init_len = len(geo)

    # If the user specifies a set of atoms to preserve as is, then
    # then bonding_pref entry is set to full saturation.
    if preserve != []:
        bonding_pref = [(i, add_hydrogens.sat_dict[elements[i]]) for i in preserve]
    else:
        bonding_pref = []

    # Get the lewis structure
    # lone_electrons,bonding_electrons,core_electrons,bonding_pref = check_lewis(atomtypes,adj_mat,q_tot=q_tot,bonding_pref=bonding_pref,return_pref=True,fixed_bonds=fixed_bonds)
    # since the input doesn't contain fc infor yet, set fc all be zero
    fc = [0] * len(elements)
    lone_electrons, bonding_electrons, core_electrons, bonding_pref = molecule.AtomType.frag_find_lewis(q_tot=molecule.q_tot,
                                                                                      fc_0=fc, keep_lone=[],
                                                                                      fixed_bonds=fixed_bonds,
                                                                                      bonding_pref=bonding_pref,
                                                                                      return_pref=True,
                                                                                      check_lewis_flag=True)

    # Update the preserved atoms (check_lewis will extend this list if there are special groups (e.g., nitro moieties) that need to be conserved for the sake of the lewis structure)
    preserve = set([i[0] for i in bonding_pref])

    # Loop over the atoms in the geometry
    for count_i, i in enumerate(geo):

        # ID undercoordinated atoms
        if count_i in preserve:
            continue
        elif add_hydrogens.sat_dict[elements[count_i]] is not None:
            B_expected = add_hydrogens.sat_dict[elements[count_i]]
        else:
            raise FragGenAllException(
                "ERROR in add_hydrogens: could not determine the number of hydrogens to add to {}. Exiting...".format(
                    elements[count_i]))
        B_current = bonding_electrons[count_i]

        # Determine the number of nuclei that are attached and expected.
        N_current = np.sum(adj_mat[count_i])
        N_expected = N_current + (B_expected - B_current)

        # Add hydrogens to undercoordinated atoms
        if N_expected > N_current:

            old_inds = [count_j for count_j, j in enumerate(adj_mat[count_i]) if j == 1]

            # Protocols for 1 missing hydrogen
            if N_expected - N_current == 1:
                if N_expected == 1:
                    new = i + np.array([H_length, 0.0, 0.0])
                elif N_expected == 2:
                    new = -1.0 * normalize(geo[old_inds[0]] - i) * H_length + i + np.array(
                        [random.random(), random.random(),
                         random.random()]) * 0.01  # random factor added for non-carbon types to relax during FF-opt
                elif N_expected == 3:
                    new = -1.0 * normalize(
                        normalize(geo[old_inds[0]] - i) + normalize(geo[old_inds[1]] - i)) * H_length + i
                elif N_expected == 4:
                    print("WARNING in add_hydrogens: 1:4 (unchecked)")
                    new = -1.0 * normalize(
                        normalize(geo[old_inds[0]] - i) + normalize(geo[old_inds[1]] - i) + normalize(
                            geo[old_inds[2]] - i)) * H_length + i

                    # Update geometry, adj_mat, elements, and atomtypes with one new atoms
                geo = np.vstack([geo, new])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])]
                elements += ["H"]
                tmp = np.zeros([N_atoms + 1, N_atoms + 1])
                tmp[:N_atoms, :N_atoms] = adj_mat
                tmp[-1, count_i] = 1
                tmp[count_i, -1] = 1
                adj_mat = tmp
                N_atoms += 1

            # Protocols for 2 missing hydrogens
            # ISSUE, NEW ALGORITHM IS BASED ON BONDED ATOMS NOT BONDED CENTERS
            if N_expected - N_current == 2:
                if N_expected == 2:
                    new_1 = i + np.array([H_length, 0.0, 0.0])
                    new_2 = i - np.array([H_length, 0.0, 0.0])
                elif N_expected == 3:
                    rot_vec = normalize(
                        np.cross(geo[old_inds[0]] - i, np.array([random.random(), random.random(), random.random()])))
                    new_1 = normalize(axis_rot(geo[old_inds[0]], rot_vec, i, 120.0) - i) * H_length + i
                    new_2 = normalize(axis_rot(geo[old_inds[0]], rot_vec, i, 240.0) - i) * H_length + i
                elif N_expected == 4:
                    bisector = normalize(geo[old_inds[0]] - i + geo[old_inds[1]] - i)
                    new_1 = axis_rot(geo[old_inds[0]], bisector, i, 90.0)
                    new_2 = axis_rot(geo[old_inds[1]], bisector, i, 90.0)
                    rot_vec = normalize(np.cross(new_1 - i, new_2 - i))
                    angle = (109.5 - acos(np.dot(normalize(new_1 - i), normalize(new_2 - i))) * 180.0 / pi) / 2.0
                    new_1 = axis_rot(new_1, rot_vec, i, -angle)
                    new_2 = axis_rot(new_2, rot_vec, i, angle)
                    new_1 = -1 * H_length * normalize(new_1 - i) + i
                    new_2 = -1 * H_length * normalize(new_2 - i) + i

                # Update geometry, adj_mat, elements, and atomtypes with two new atoms
                geo = np.vstack([geo, new_1])
                geo = np.vstack([geo, new_2])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])] * 2
                elements += ["H", "H"]
                tmp = np.zeros([N_atoms + 2, N_atoms + 2])
                tmp[:N_atoms, :N_atoms] = adj_mat
                tmp[[-1, -2], count_i] = 1
                tmp[count_i, [-1, -2]] = 1
                adj_mat = tmp
                N_atoms += 2

            # Protocols for 3 missing hydrogens
            if N_expected - N_current == 3:
                if N_expected == 3:
                    rot_vec = np.array([0.0, 1.0, 0.0])
                    new_1 = i + np.array([H_length, 0.0, 0.0])
                    new_2 = axis_rot(new_1, rot_vec, i, 120.0)
                    new_3 = axis_rot(new_1, rot_vec, i, 240.0)
                if N_expected == 4:
                    rot_vec = normalize(
                        np.cross(geo[old_inds[0]] - i, np.array([random.random(), random.random(), random.random()])))
                    new_1 = H_length * normalize(axis_rot(geo[old_inds[0]], rot_vec, i, 109.5) - i) + i
                    new_2 = axis_rot(new_1, normalize(i - geo[old_inds[0]]), i, 120.0)
                    new_3 = axis_rot(new_2, normalize(i - geo[old_inds[0]]), i, 120.0)

                # Update geometry, adj_mat, elements, and atomtypes with three new atoms
                geo = np.vstack([geo, new_1])
                geo = np.vstack([geo, new_2])
                geo = np.vstack([geo, new_3])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])] * 3
                elements += ["H", "H", "H"]
                tmp = np.zeros([N_atoms + 3, N_atoms + 3])
                tmp[:N_atoms, :N_atoms] = adj_mat
                tmp[[-1, -2, -3], count_i] = 1
                tmp[count_i, [-1, -2, -3]] = 1
                adj_mat = tmp
                N_atoms += 3

            # Protocols for 4 missing hydrogens
            if N_expected - N_current == 4:
                if N_expected == 4:
                    new_1 = i + np.array([H_length, 0.0, 0.0])
                    rot_vec = normalize(
                        np.cross(new_1 - i, np.array([random.random(), random.random(), random.random()])))
                    new_2 = H_length * normalize(axis_rot(new_1, rot_vec, i, 109.5) - i) + i
                    new_3 = axis_rot(new_2, normalize(i - new_1), i, 120.0)
                    new_4 = axis_rot(new_3, normalize(i - new_1), i, 120.0)

                # Update geometry, adj_mat, elements, and atomtypes with three new atoms
                geo = np.vstack([geo, new_1])
                geo = np.vstack([geo, new_2])
                geo = np.vstack([geo, new_3])
                geo = np.vstack([geo, new_4])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])] * 4
                elements += ["H", "H", "H", "H"]
                tmp = np.zeros([N_atoms + 4, N_atoms + 4])
                tmp[:N_atoms, :N_atoms] = adj_mat
                tmp[[-1, -2, -3, -4], count_i] = 1
                tmp[count_i, [-1, -2, -3, -4]] = 1
                adj_mat = tmp
                N_atoms += 4
    # Note: gens here is temporary
    N_molecule = Molecule(molecule.gens)
    N_molecule.parse_data(geometry=geo, atom_types=atomtypes, elements=elements, adj_mat=adj_mat)

    if retype is True:
        N_molecule.AtomType.UpdateAtomTypes()
        return N_molecule, range(init_len, len(geo))
    else:
        return N_molecule, range(init_len, len(geo))


# only used by frag_gen*.py
def new_mode_geo(m_ind, molecule: Molecule, force_linear=False):
    """
    Description: Returns the canonical fragment corresponding to the mode defined associated with geo and atoms m_ind.

    Inputs:      m_ind:        list of indices involved in the mode
                 geo:          an Nx3 np.array holding the geometry of the molecule
                 adj_mat:      an NxN np.array holding the connectivity of the molecule
                 gens:         an integer specifying the number of generations involved in the geometry search
                                (Algorithm returns
                 force_linear: boolean, forces a non-cyclic structure.

    Returns:     m_ind:        list of indices involved in the mode (in terms of the new geometry)
                 N_geo:        new geometry for parameterizing the mode.
                 N_adj_mat:    new adjacency matrix
                 N_dup:        user supplied lists indexed to the original geometry, now indexed to the new geometry
    """
    m_ind_0 = deepcopy(m_ind)
    # Seed conditions for...
    # atoms: single atom
    # bonds: both atoms
    # angles: center atom
    # linear dihedral: center atoms
    # improper dihedral: center atom
    # b_top: holds the topology of the mode (which atoms are bonded to which), only used for remapping when using the force_linear algorithm
    if len(m_ind) == 1 or len(m_ind) == 2:
        m_ind = m_ind
        b_top = [[]]
    elif len(m_ind) == 2:  # XXX why is this here? 7/10/20
        m_ind = m_ind
        b_top = [[1], [0]]
    elif len(m_ind) == 3:
        m_ind = [m_ind[1]]
        b_top = [[1], [0, 2], [1]]
    elif len(m_ind) == 4:
        # Check for linear condition (a bond between the center atoms)
        if molecule.adj_mat[m_ind[1], m_ind[2]] == 1 and molecule.adj_mat[m_ind[2], m_ind[3]] == 1:
            m_ind = [m_ind[1], m_ind[2]]
            b_top = [[1], [0, 2], [1, 3], [2]]
            force_linear = True

        # Treat it as an improper
        else:
            m_ind = [m_ind[0]]
            b_top = [[1, 2, 3], [0], [0], [0]]
    #
    # Creates a non-linear but constitutionally identical topology. For linear structures there is no difference between the if/else statements
    #
    if force_linear is True:

        new_atoms = []
        bonds = [(m_ind.index(i), m_ind.index(count_j)) for i in m_ind for count_j, j in enumerate(molecule.adj_mat[i])
                 if
                 j == 1 and count_j in m_ind]  # holds the bonds in the new adj mat
        original_bonds = [(i, count_j) for i in m_ind for count_j, j in enumerate(molecule.adj_mat[i]) if
                          j == 1 and count_j in m_ind]  # holds the original bonds in the new adj mat
        parent_next = [[count_j for count_j, j in enumerate(molecule.adj_mat[i]) if j == 1 and count_j in m_ind] for i
                       in
                       m_ind]  # holds the "parents" of each added atom to avoid backtracking
        for i in range(molecule.gens + 1):

            # Used for indexing the bonds
            current_length = len(new_atoms)

            # Update the atoms that are included in the geometry (new_atoms)
            # as well as the connections that need to be iterated over (cons)
            # on the first generation the cons are seeded with m_ind
            if i == 0:
                new_atoms += m_ind
                cons = m_ind
            else:
                new_atoms += cons_next
                cons = cons_next

            # Exits after adding the final generation of connections
            if i == molecule.gens:
                break

            # Refresh the parent and "origin" lists. "origin" is different from parent because it
            # stores the index in new_atoms of each node in cons, this is needed for properly indexing the bonds.
            parent = parent_next
            cons_next = []
            parent_next = []
            origins = [current_length + j for j in range(len(cons))]

            # Loop over the connections at the current generation, find the next generation of
            # connections and avoid backtracking.
            for count_j, j in enumerate(cons):
                for count_k, k in enumerate(molecule.adj_mat[j]):
                    if k == 1 and count_k not in parent[count_j]:
                        cons_next += [count_k]
                        bonds += [(origins[count_j], len(new_atoms) + len(cons_next) - 1)]
                        original_bonds += [(j, count_k)]
                        parent_next += [[j]]

        # Generate the new adjacency matrix
        N_adj_mat = np.zeros([len(new_atoms), len(new_atoms)])

        # determine single or multiple bond_mats
        if len(molecule.bond_mat) > 0 and len(molecule.bond_mat) != len(molecule.geometry):  # multuple bond_mats
            N_bond_mats = []
            for count_i, i in enumerate(bonds):
                N_adj_mat[i[0], i[1]] = 1
                N_adj_mat[i[1], i[0]] = 1

            for lb in range(len(molecule.bond_mat)):
                N_bond_mat = np.zeros([len(new_atoms), len(new_atoms)])
                for count_i, i in enumerate(bonds):
                    N_bond_mat[i[0], i[1]] = molecule.bond_mat[lb][
                        original_bonds[count_i][0], original_bonds[count_i][1]]
                    N_bond_mat[i[1], i[0]] = molecule.bond_mat[lb][
                        original_bonds[count_i][0], original_bonds[count_i][1]]
                N_bond_mats += [N_bond_mat]
        else:
            N_bond_mat = np.zeros([len(new_atoms), len(new_atoms)])
            for count_i, i in enumerate(bonds):
                N_adj_mat[i[0], i[1]] = 1
                N_adj_mat[i[1], i[0]] = 1
                N_bond_mat[i[0], i[1]] = molecule.bond_mat[original_bonds[count_i][0], original_bonds[count_i][1]]
                N_bond_mat[i[1], i[0]] = molecule.bond_mat[original_bonds[count_i][0], original_bonds[count_i][1]]

    #
    # Use graph_seps algorithm: keeps all atoms and the topology out to gens bonds
    #
    else:

        # Graphical separations are used for determining which atoms and bonds to keep
        gs = molecule.AdjMat.graph_seps()
        # all atoms within "gens" of the m_ind atoms are kept
        # all bonds within "gens" of the m_ind atoms are kept (i.e. bonds between the "gens" separated atoms are NOT kept
        new_atoms = list(set([count_j for i in m_ind for count_j, j in enumerate(gs[i]) if j <= molecule.gens]))
        N_adj_mat = molecule.adj_mat[new_atoms, :][:, new_atoms]
        # determine single or multiple bond_mats
        if len(molecule.bond_mat) > 0 and len(molecule.bond_mat) != len(molecule.geometry):  # multuple bond_mats
            N_bond_mats = []
            for lb in range(len(molecule.bond_mat)):
                N_bond_mat = molecule.bond_mat[lb][new_atoms, :][:, new_atoms]
                N_bond_mats += [N_bond_mat]
        else:
            N_bond_mat = molecule.bond_mat[new_atoms, :][:, new_atoms]

        # remove the bonds between the "gens" separated atoms
        edge_ind = list(set([count_j for i in m_ind for count_j, j in enumerate(gs[i]) if j == molecule.gens]))
        edge_ind = [new_atoms.index(i) for i in edge_ind if min([gs[j, i] for j in m_ind]) == molecule.gens]
        for i in edge_ind:
            for j in edge_ind:
                N_adj_mat[i, j] = 0
                N_adj_mat[j, i] = 0
                N_bond_mat[i, j] = 0
                N_bond_mat[j, i] = 0

    # Create the new geometry and adj_mat
    N_geo = np.zeros([len(new_atoms), 3])
    for count_i, i in enumerate(new_atoms):
        N_geo[count_i, :] = molecule.geometry[i, :]

    N_elements = []
    N_atomtypes = []
    for i in new_atoms:
        N_elements.append(molecule.elements[i])
        N_atomtypes.append(molecule.atom_types[i])

    # Clean up the geometry
    N_geo = transify(N_geo, N_adj_mat, opt_terminals=False, opt_final=False)

    # Assign the mode ind
    # NOTE: the use of the list.index() method assumes that the mode indices are the first occuring in the geometry
    #       this should be a very good assumption for all conventional modes and seed scenarios (no exceptions have been found).
    if force_linear is True:

        # A more complicated remap is required for the force_linear algorithm because care needs to be
        # taken to ensure that the essential topology of m_ind_0 is preserved.

        # inds holds the mapped m_ind and it is initially seeded with the either
        # "None" values or the values from m_ind_0 that seeded the mode search
        inds = [None if i not in m_ind else new_atoms.index(i) for i in m_ind_0]

        # calculate the values for inds that preserve the original topology of m_ind
        for count_i, i in enumerate(inds):

            # Skip if already assigned
            if i is not None:
                continue

            # Find the indices in the old adjacency matrix for which m_ind_0[count_i] had connections with other m_ind_0 atoms.
            must_be_connected_to_old = [count_j for count_j, j in enumerate(molecule.adj_mat[m_ind_0[count_i]]) if
                                        j == 1 and count_j in m_ind_0]

            # Find connections in inds that have already been placed and must be presereved for assigning inds[count_i]
            must_be_connected_to_new = [inds[j] for j in b_top[count_i] if inds[j] is not None]

            # Find all of the atoms that have the right originating index
            candidates = [count_j for count_j, j in enumerate(new_atoms) if j == m_ind_0[count_i]]
            for j in candidates:

                # Check if the candidate j has the right connections in terms of
                # (i) the atom types it is connected to and
                # (ii) connections with respect to atoms in inds that have already been assigned
                connections = [new_atoms[count_k] for count_k, k in enumerate(N_adj_mat[j]) if k == 1]
                if False not in [k in connections for k in must_be_connected_to_old] and False not in [
                    N_adj_mat[j, k] == 1 for k in must_be_connected_to_new]:
                    inds[count_i] = j
                    break

        # Assign the mapping
        m_ind = inds

    else:
        m_ind = [new_atoms.index(i) for i in m_ind_0]

    # Note: gens here is temporary
    N_molecule = Molecule(molecule.gens)
    N_molecule.parse_data(elements=N_elements, geometry=N_geo, atom_types=N_atomtypes, adj_mat=N_adj_mat)
    if len(molecule.bond_mat) > 0 and len(molecule.bond_mat) != len(molecule.geometry):  # multuple bond_mats
        N_molecule.bond_mat = N_bond_mats
    else:
        N_molecule.bond_mat = N_bond_mat
    return m_ind, N_molecule


def map_frag(fixed_molecule, remapped_molecule):
    """
    This function handles remapping the atomtypes of one fragment on to another such that the modes (bonds
    angles, dihedrals, etc) in the remapped fragment are the same after remapping.
    """
    # The equivalency of the two topologies is based upon the mode types and number each topology possesses
    bond_types_0, angle_types_0, dihedral_types_0, one_five_types_0 = remapped_molecule.find_modes(return_all=1)[4:]
    comp_obj = (sorted(bond_types_0), sorted(angle_types_0), sorted(dihedral_types_0))

    # Seed the mapping by placing the first atom in the mappable set onto the first instance
    # of its hash type in the fixed set.
    # first = len(remapped_hash)-1
    first = 0
    mapping = {first: fixed_molecule.hash_list.index(
        remapped_molecule.hash_list[first])}  # fixed indices keyed to remapped indices
    R_mapping = {fixed_molecule.hash_list.index(
        remapped_molecule.hash_list[first]): first}  # remapped indices keyed to fixed indices
    remapped_place_me = set(range(len(remapped_molecule.atom_types)))  # remapped atom indices that need to be placed
    fixed_place_me = set(
        range(len(remapped_molecule.atom_types)))  # fixed atom indices that haven't found a match in the remapped set
    remapped_place_me.remove(first)  # remove the starting atom index
    fixed_place_me.remove(mapping[first])  # remove the starting atom index

    # Seed the connections list (idx in fixed, connection restraint in remapped) and iterate over the connections until no new connections are found
    cons = [(count_i, first) for count_i, i in enumerate(fixed_molecule.adj_mat[mapping[first]]) if
            i == 1 and count_i in fixed_place_me]
    for i in cons:
        if i[0] not in fixed_place_me: continue
        remapped_match = next(
            j for j in remapped_place_me if
            remapped_molecule.hash_list[j] == fixed_molecule.hash_list[i[0]] and remapped_molecule.adj_mat[
                i[1], j] == 1)
        mapping[remapped_match] = i[0]
        R_mapping[i[0]] = remapped_match
        fixed_place_me.remove(i[0])
        remapped_place_me.remove(remapped_match)
        cons += [(count_j, remapped_match) for count_j, j in enumerate(fixed_molecule.adj_mat[i[0]]) if
                 j == 1 and count_j in fixed_place_me]

        # Remap the atomtypes
    new_atomtypes = ["X"] * len(remapped_molecule.atom_types)
    for i in mapping: new_atomtypes[mapping[i]] = remapped_molecule.atom_types[i]

    # Check if the set of types and number of each mode match the original fragment
    # NOTE: the fixed_adj_mat is used for determining the connectivity
    # NOTE: the fixed_bond_mats is used for determining dihedral types
    tmp_mol = Molecule(fixed_molecule.gens)
    tmp_mol.parse_data(adj_mat=fixed_molecule.adj_mat, atom_types=new_atomtypes, bond_mat=fixed_molecule.bond_mat)
    bond_types_1, angle_types_1, dihedral_types_1, one_five_types_1 = tmp_mol.find_modes(return_all=1)[4:]
    if comp_obj != (sorted(bond_types_1), sorted(angle_types_1), sorted(dihedral_types_1)):
        raise FragGenAllException("ERROR in map_frag: No mapping was discovered that preserved the topology of the "
                                  "original fragment. Check that the two fragments are indeed identical. Exiting...")
        print("{:60s} {:60s}".format("Original", "New"))
        for count_i, i in enumerate(remapped_molecule.atom_types):
            print("{:<3d} {:60s} {:60s}".format(count_i, i, new_atomtypes[count_i]))
        return new_atomtypes, mapping, 1

    return new_atomtypes, mapping, 0


# Logger object redirects standard output to a file.
class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("frag_gen_all.log", "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        pass

if __name__ == "__main__":
    main(sys.argv[1:])
