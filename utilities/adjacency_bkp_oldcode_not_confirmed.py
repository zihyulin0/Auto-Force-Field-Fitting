import numpy as np
from copy import deepcopy
from utilities.PeriodictTable import Tables
from Excpetions import TAFFIException
from utilities.parse import read_alljson
from scipy.spatial.distance import cdist
import os, random
from model.structure import StructureBase

class AdjacencyException(TAFFIException):
    pass
class AdjacencyMatrix(StructureBase):
    """
    Operations related to adjacency matrix
    """
    def __init__(self):
        super().__init__()
        # adj_mat needs to be protected to avoid being accidentally changed when being inherited
        self._adj_mat = None

    @property
    def adj_mat(self):
        return self._adj_mat

    def BuildAdjMat(self):
        self._adj_mat = self.Table_generator()

    def graph_seps(self):
        """
        Returns a matrix of graphical separations for all nodes in a graph defined by the inputted adjacency matrix
        """

        # Create a new name for the object holding A**(N), initialized with A**(1)
        adj_mat = deepcopy(self.adj_mat)

        # Initialize an np.array to hold the graphical separations with -1 for all unassigned elements and 0 for the diagonal.
        seps = np.ones([len(adj_mat), len(adj_mat)]) * -1
        np.fill_diagonal(seps, 0)

        # Perform searches out to len(adj_mat) bonds (maximum distance for a graph with len(adj_mat) nodes
        for i in np.arange(len(adj_mat)):

            # All perform assignments to unassigned elements (seps==-1)
            # and all perform an assignment if the value in the adj_mat is > 0
            seps[np.where((seps == -1) & (adj_mat > 0))] = i + 1

            # Since we only care about the leading edge of the search and not the actual number of paths at higher orders, we can
            # set the larger than 1 values to 1. This ensures numerical stability for larger adjacency matrices.
            adj_mat[np.where(adj_mat > 1)] = 1

            # Break once all of the elements have been assigned
            if -1 not in seps:
                break

            # Take the inner product of the A**(i+1) with A**(1)
            adj_mat = np.dot(adj_mat, self.adj_mat)

        return seps

    def Table_generator(self, print_file=None, Radii_dict=False):
        """
        Generates the adjacency matrix based on UFF bond radii
        Inputs:     Elements: N-element List strings for each atom type
                    Geometry: Nx3 np.array holding the geometry of the molecule
                    File:  Optional. If Table_generator encounters a problem then it is often useful to have the name of the file the geometry came from printed.
        """
        Radii = Tables.Radii

        # Use Radii json file in Lib folder if sepcified
        if Radii_dict:
            # TODO check this module path
            module_path = '/'.join(os.path.abspath(__file__).split('/')[:-1])
            if os.path.isfile(module_path + '/Radii.json') is False:
                raise AdjacencyException("ERROR: {}/Radii.json doesn't exist, check for Radii.json file in the Library".format(module_path))
            Radii = read_alljson(module_path + '/Radii.json')


        Max_Bonds = Tables.MAXBONDS

        # Scale factor is used for determining the bonding threshold. 1.2 is a heuristic that give some lattitude in defining bonds since the UFF radii correspond to equilibrium lengths.
        scale_factor = 1.2

        # Print warning for uncoded elements.
        for i in self.elements:
            if i not in Radii.keys():
                raise AdjacencyException(
                    'ERROR in Table_generator: The geometry contains an element ({}) that the Table_generator function doesn\'t have bonding information for. This needs to be directly added to the Radii dictionary before proceeding. Exiting...'.format(i))

        # Generate distance matrix holding atom-atom separations (only save upper right)
        Dist_Mat = np.triu(cdist(self.geometry, self.geometry))

        # Find plausible connections
        x_ind, y_ind = np.where((Dist_Mat > 0.0) & (Dist_Mat < max([Radii[i] ** 2.0 for i in Radii.keys()])))

        # Initialize Adjacency Matrix
        self._adj_mat = np.zeros([len(self.geometry), len(self.geometry)])

        # Iterate over plausible connections and determine actual connections
        for count, i in enumerate(x_ind):

            # Assign connection if the ij separation is less than the UFF-sigma value times the scaling factor
            if Dist_Mat[i, y_ind[count]] < (Radii[self.elements[i]] + Radii[self.elements[y_ind[count]]]) * scale_factor:
                self._adj_mat[i, y_ind[count]] = 1

            if self.elements[i] == 'H' and self.elements[y_ind[count]] == 'H':
                if Dist_Mat[i, y_ind[count]] < (Radii[self.elements[i]] + Radii[self.elements[y_ind[count]]]) * 1.5:
                    self._adj_mat[i, y_ind[count]] = 1

        # Hermitize Adj_mat
        self._adj_mat = self._adj_mat + self._adj_mat.transpose()

        # Perform some simple checks on bonding to catch errors
        problem_dict = {i: 0 for i in Radii.keys()}
        for count_i, i in enumerate(self._adj_mat):

            if Max_Bonds[self.elements[count_i]] is not None and sum(i) > Max_Bonds[self.elements[count_i]]:
                problem_dict[self.elements[count_i]] += 1
                cons = sorted([(Dist_Mat[count_i, count_j], count_j) if count_j > count_i else (
                Dist_Mat[count_j, count_i], count_j) for count_j, j in enumerate(i) if j == 1])[::-1]
                while sum(self._adj_mat[count_i]) > Max_Bonds[self.elements[count_i]]:
                    sep, idx = cons.pop(0)
                    self._adj_mat[count_i, idx] = 0
                    self._adj_mat[idx, count_i] = 0
        #        if Elements[count_i] in conditions.keys():
        #            if sum(i) > conditions[Elements[count_i]]:

        # Print warning messages for obviously suspicious bonding motifs.
        element_to_name = {'H': 'hydogen', 'C': 'carbon', 'Si': 'silicon', 'F': 'fluorine', 'Cl': 'chlorine',
                           'Br': 'bromine', 'I': 'iodine', 'O': 'oxygen', 'N': 'nitrogen', 'B': 'boron'}
        conditions = {"H": 1, "C": 4, "Si": 4, "F": 1, "Cl": 1, "Br": 1, "I": 1, "O": 2, "N": 4, "B": 4}

        if sum([problem_dict[i] for i in problem_dict.keys()]) > 0:
            print("Table Generation Warnings:")
            if print_file is None:
                warning_msg = 'WARNING in Table_generator: {} {}(s) have more than {} bond.'
            else:
                warning_msg = 'WARNING in Table_generator: parsing {}'.format(
                    self.xyz) + ' {} {}(s) have more than {} bond.'
            for i in sorted(problem_dict.keys()):
                if problem_dict[i] > 0:
                    if i in conditions.keys():
                        print(warning_msg.format(problem_dict[i],element_to_name[i],conditions[i]))
            print("")

        return self._adj_mat


    def find_lewis(self, bonding_pref=[], fixed_bonds=[], fc_0=None, keep_lone=[],
                   return_pref=False, verbose=False, b_mat_only=False, return_FC=False, octet_opt=True,
                   check_lewis_flag=False):
        """
        # Returns a list with the number of electrons on each atom and a list with the number missing/surplus electrons on the atom
        #
        # Inputs:  elements:  a list of element labels indexed to the adj_mat
        #          adj_mat:   np.array of atomic connections
        #          bonding_pref: optional list of (index, bond_number) tuples that sets the target bond number of the indexed atoms
        #          q_tot:     total charge on the molecule
        #          fixed_bonds: optional list of (index_1,index_2,bond_number) tuples that creates fixed bonds between the index_1
        #                       and index_2 atoms. No further bonds will be added or subtracted between these atoms.
        #
        # Optional inputs for ion and radical cases:
        #          fc_0:      a list of formal charges on each atom
        #          keep_lone: a list of atom index for which contains a radical
        #
        # Returns: lone_electrons:
        #          bonding_electrons:
        #          core_electrons:
        #          bond_mat:  an NxN matrix holding the bond orders between all atoms in the adj_mat
        #          bonding_pref (optinal): optional list of (index, bond_number) tuples that sets the target bond number of the indexed atoms
        #
        """

        # Initialize the preferred lone electron dictionary the first time this function is called
        # TODO alternatively we can just acccess this from Tables, but I guess this way we sort of preserve the ability to change it
        # I don't the way the red pops out in pycharmm for this
        if not hasattr(find_lewis, "sat_dict"):
            find_lewis.lone_e = Tables.LONE_ELETRON
            # Initialize periodic table
            find_lewis.periodic = Tables.PERIODICT
            # Electronegativity ordering (for determining lewis structure)
            find_lewis.en = Tables.ELETRONEGATIVITY
            # Polarizability ordering (for determining lewis structure)
            find_lewis.pol = Tables.POLARIZABILITY
            # Bond Energy
            find_lewis.be = Tables.BOND_ENERGY
            # Initialize periodic table
            find_lewis.atomic_to_element = {find_lewis.periodic[i]: i for i in find_lewis.periodic.keys()}

        # somewhere here sets q_tot to zero which I don't like the member be
        # alternated especially when it's in base class
        q_tot = deepcopy(self.q_tot)

        # Consistency check on fc_0 argument, if supplied
        if fc_0 is not None:

            if len(fc_0) != len(self.elements):
                raise AdjacencyException("ERROR in find_lewis: the fc_0 and elements lists must have the same dimensions.")

            if int(sum(fc_0)) != int(q_tot):
                raise AdjacencyException("ERROR in find_lewis: the sum of formal charges does not equal q_tot.")

        # Initalize elementa and atomic_number lists for use by the function
        atomic_number = [find_lewis.periodic[i.lower()] for i in self.elements]
        adj_mat = deepcopy(self.adj_mat)

        # identify ring atoms and number of rings
        rings = []
        ring_size_list = range(11)[3:]  # at most 10 ring structure

        for ring_size in ring_size_list:
            for j, Ej in enumerate(self.elements):
                is_ring, ring_ind = self.return_ring_atom(j, ring_size=ring_size, alt_adj_mat=adj_mat)
                if is_ring and ring_ind not in rings:
                    rings += [ring_ind]
        rings = [list(ring_inds) for ring_inds in rings]

        # Initially assign all valence electrons as lone electrons
        lone_electrons = np.zeros(len(self.elements), dtype="int")
        bonding_electrons = np.zeros(len(self.elements), dtype="int")
        core_electrons = np.zeros(len(self.elements), dtype="int")
        valence = np.zeros(len(self.elements), dtype="int")
        bonding_target = np.zeros(len(self.elements), dtype="int")
        valence_list = np.zeros(len(self.elements), dtype="int")

        for count_i, i in enumerate(self.elements):

            # Grab the total number of (expected) electrons from the atomic number
            N_tot = atomic_number[count_i]

            # Determine the number of core/valence electrons based on row in the periodic table
            if N_tot > 54:
                raise AdjacencyException("ERROR in find_lewis: the algorithm isn't compatible with atomic numbers "
                                         "greater than 54 owing to a lack of rules for treating lanthanides. "
                                         "Exiting...")

            elif N_tot > 36:
                N_tot -= 36
                core_electrons[count_i] = 36
                valence[count_i] = 18

            elif N_tot > 18:
                N_tot -= 18
                core_electrons[count_i] = 18
                valence[count_i] = 18

            elif N_tot > 10:
                N_tot -= 10
                core_electrons[count_i] = 10
                valence[count_i] = 8

            elif N_tot > 2:
                N_tot -= 2
                core_electrons[count_i] = 2
                valence[count_i] = 8

            lone_electrons[count_i] = N_tot
            valence_list[count_i] = N_tot

            # Assign target number of bonds for this atom
            if count_i in [j[0] for j in bonding_pref]:
                bonding_target[count_i] = next(j[1] for j in bonding_pref if j[0] == count_i)
            else:
                bonding_target[count_i] = N_tot - find_lewis.lone_e[self.elements[count_i].lower()]

                # Loop over the adjmat and assign initial bonded electrons assuming single bonds (and adjust lone electrons accordingly)
        for count_i, i in enumerate(self.adj_mat):
            bonding_electrons[count_i] += sum(i)
            lone_electrons[count_i] -= sum(i)

        # Apply keep_lone: add one electron to such index
        for count_i in keep_lone:
            lone_electrons[count_i] += 1

        # Eliminate all radicals by forming higher order bonds
        change_list = range(len(lone_electrons))
        bonds_made = []
        loop_list = [(atomic_number[i], i) for i in range(len(lone_electrons))]
        loop_list = [i[1] for i in sorted(loop_list)]

        # Check for special chemical groups
        for i in range(len(self.elements)):
            # Handle nitro groups
            if self.is_nitro(i) is True:
                O_ind = [count_j for count_j, j in enumerate(self.adj_mat[i]) if
                         j == 1 and self.elements[count_j].lower() == "o" and sum(self.adj_mat[count_j]) == 1]
                bonding_pref = [j for j in bonding_pref if j[0] != i and j[0] not in O_ind]
                bonding_pref += [(i, 4)]
                bonding_pref += [(O_ind[0], 1)]
                bonding_pref += [(O_ind[1], 2)]
                bonding_electrons[O_ind[1]] += 1
                bonding_electrons[i] += 1
                lone_electrons[O_ind[1]] -= 1
                lone_electrons[i] -= 2
                lone_electrons[O_ind[0]] += 1
                adj_mat[i, O_ind[1]] += 1
                adj_mat[O_ind[1], i] += 1

            # Handle sulfoxide groups
            if self.is_sulfoxide(i) is True:
                O_ind = [count_j for count_j, j in enumerate(self.adj_mat[i]) if
                         j == 1 and self.elements[count_j].lower() == "o" and sum(self.adj_mat[count_j]) == 1]
                bonding_pref = [j for j in bonding_pref if j[0] != i and j[
                    0] not in O_ind]  # remove bonds involving the thioketone atoms from the bonding_pref list
                bonding_pref += [(i, 4)]
                bonding_pref += [(O_ind[0], 2)]
                bonding_electrons[O_ind[0]] += 1
                bonding_electrons[i] += 1
                lone_electrons[O_ind[0]] -= 1
                lone_electrons[i] -= 1
                bonds_made += [(i, O_ind[0])]
                adj_mat[i, O_ind[0]] += 1
                adj_mat[O_ind[0], i] += 1

            # Handle sulfonyl groups
            if self.is_sulfonyl(i) is True:
                O_ind = [count_j for count_j, j in enumerate(self.adj_mat[i]) if
                         j == 1 and self.elements[count_j].lower() == "o" and sum(self.adj_mat[count_j]) == 1]
                bonding_pref = [j for j in bonding_pref if j[0] != i and j[
                    0] not in O_ind]  # remove bonds involving the sulfoxide atoms from the bonding_pref list
                bonding_pref += [(i, 6)]
                bonding_pref += [(O_ind[0], 2)]
                bonding_pref += [(O_ind[1], 2)]
                bonding_electrons[O_ind[0]] += 1
                bonding_electrons[O_ind[1]] += 1
                bonding_electrons[i] += 2
                lone_electrons[O_ind[0]] -= 1
                lone_electrons[O_ind[1]] -= 1
                lone_electrons[i] -= 2
                bonds_made += [(i, O_ind[0])]
                bonds_made += [(i, O_ind[1])]
                adj_mat[i, O_ind[0]] += 1
                adj_mat[i, O_ind[1]] += 1
                adj_mat[O_ind[0], i] += 1
                adj_mat[O_ind[1], i] += 1

                # Handle phosphate groups
            if self.is_phosphate(i) is True:
                O_ind = [count_j for count_j, j in enumerate(self.adj_mat[i]) if
                         j == 1 and self.elements[count_j] in ["o", "O"]]  # Index of single bonded O-P oxygens
                O_ind_term = [j for j in O_ind if sum(self.adj_mat[j]) == 1]  # Index of double bonded O-P oxygens
                bonding_pref = [j for j in bonding_pref if j[0] != i and j[
                    0] not in O_ind]  # remove bonds involving the phosphate atoms from the bonding_pref list
                bonding_pref += [(i, 5)]
                bonding_pref += [(O_ind_term[0],
                                  2)]  # during testing it ended up being important to only add a bonding_pref tuple for one of the terminal oxygens
                bonding_electrons[O_ind_term[0]] += 1
                bonding_electrons[i] += 1
                lone_electrons[O_ind_term[0]] -= 1
                lone_electrons[i] -= 1
                bonds_made += [(i, O_ind_term[0])]
                adj_mat[i, O_ind_term[0]] += 1
                adj_mat[O_ind_term[0], i] += 1

            # Handle cyano groups
            if self.is_cyano(i, alt_adj_mat=adj_mat) is True:
                C_ind = [count_j for count_j, j in enumerate(self.adj_mat[i]) if
                         j == 1 and self.elements[count_j] in ["c", "C"] and sum(self.adj_mat[count_j]) == 2]
                bonding_pref = [j for j in bonding_pref if j[0] != i and j[
                    0] not in C_ind]  # remove bonds involving the cyano atoms from the bonding_pref list
                bonding_pref += [(i, 3)]
                bonding_pref += [(C_ind[0], 4)]
                bonding_electrons[C_ind[0]] += 2
                bonding_electrons[i] += 2
                lone_electrons[C_ind[0]] -= 2
                lone_electrons[i] -= 2
                bonds_made += [(i, C_ind[0])]
                bonds_made += [(i, C_ind[0])]
                adj_mat[i, C_ind[0]] += 2
                adj_mat[C_ind[0], i] += 2

            # Handle isocyano groups
            if self.is_isocyano(i) is True:
                C_ind = [count_j for count_j, j in enumerate(adj_mat[i]) if
                         j == 1 and self.elements[count_j] in ["c", "C"] and sum(adj_mat[count_j]) == 1]
                bonding_pref = [j for j in bonding_pref if j[0] != i and j[
                    0] not in C_ind]  # remove bonds involving the cyano atoms from the bonding_pref list
                bonding_pref += [(i, 4)]
                bonding_pref += [(C_ind[0], 3)]
                bonding_electrons[C_ind[0]] += 2
                bonding_electrons[i] += 2
                lone_electrons[C_ind[0]] -= 2
                lone_electrons[i] -= 2
                bonds_made += [(i, C_ind[0])]
                bonds_made += [(i, C_ind[0])]
                adj_mat[i, C_ind[0]] += 2
                adj_mat[C_ind[0], i] += 2

        # Apply fixed_bonds argument
        off_limits = []
        for i in fixed_bonds:

            # Initalize intermediate variables
            a = i[0]
            b = i[1]
            N = i[2]
            N_current = len([j for j in bonds_made if (a, b) == j or (b, a) == j]) + 1
            # Check that a bond exists between these atoms in the adjacency matrix
            if self.adj_mat[a, b] != 1:
                raise AdjacencyException("ERROR in find_lewis: fixed_bonds requests bond creation between \
                                         atoms {} and {} ({} bonds) but the adjacency matrix doesn't \
                                         reflect a bond. Exiting...".format(a, b, N))

            # Check that less than or an equal number of bonds exist between these atoms than is requested
            if N_current > N:
                raise AdjacencyException("ERROR in find_lewis: fixed_bonds requests bond creation between \
                                         atoms {} and {} ({} bonds)  but {} bonds already exist between these \
                                         atoms. There may be a conflict between the special groups handling \
                                         and the requested lewis_structure.".format(a, b, N, N_current))

            # Check that enough lone electrons exists on each atom to reach the target bond number
            if lone_electrons[a] < (N - N_current):
                print(
                    "Warning in find_lewis: fixed_bonds requests bond creation between atoms {} and {} ({} bonds) \
                    but atom {} only has {} lone electrons.".format(a, b, N, self.elements[a], lone_electrons[a]))

            # Check that enough lone electrons exists on each atom to reach the target bond number
            if lone_electrons[b] < (N - N_current):
                print(
                    "Warning in find_lewis: fixed_bonds requests bond creation between atoms {} and {} ({} bonds) \
                    but atom {} only has {} lone electrons.".format(a, b, N, self.elements[b], lone_electrons[b]))

            # Make the bonds between the atoms
            for j in range(N - N_current):
                bonding_electrons[a] += 1
                bonding_electrons[b] += 1
                lone_electrons[a] -= 1
                lone_electrons[b] -= 1
                bonds_made += [(a, b)]

            # Append bond to off_limits group so that further bond additions/breaks do not occur.
            off_limits += [(a, b), (b, a)]

        # Turn the off_limits list into a set for rapid lookup
        off_limits = set(off_limits)

        # Adjust formal charges (if supplied)
        if fc_0 is not None:
            for count_i, i in enumerate(fc_0):
                if i > 0:
                    # if lone_electrons[count_i] < i:
                    # print "ERROR in find_lewis: atom ({}, index {}) doesn't have enough lone electrons ({}) to be removed to satisfy the specified formal charge ({}).".format(elements[count_i],count_i,lone_electrons[count_i],i)
                    # quit()
                    lone_electrons[count_i] = lone_electrons[count_i] - i
                if i < 0:
                    lone_electrons[count_i] = lone_electrons[count_i] + int(abs(i))
            q_tot = 0

        # diagnostic print
        if verbose is True:
            print("Starting electronic structure:")
            print("\n{:40s} {:20} {:20} {:20} {:20} {}".format("elements", "lone_electrons", "bonding_electrons",
                                                               "core_electrons", "formal_charge", "bonded_atoms"))
            for count_i, i in enumerate(self.elements):
                print("{:40s} {:<20d} {:<20d} {:<20d} {:<20d} {}".format(self.elements[count_i], lone_electrons[count_i],
                                                                         bonding_electrons[count_i],
                                                                         core_electrons[count_i], \
                                                                         valence_list[count_i] - bonding_electrons[
                                                                             count_i] - lone_electrons[count_i], \
                                                                         ",".join(
                                                                             ["{}".format(count_j) for count_j, j in
                                                                              enumerate(adj_mat[count_i]) if j == 1])))

        # Initialize objects for use in the algorithm
        lewis_total = 1000
        lewis_lone_electrons = []
        lewis_bonding_electrons = []
        lewis_core_electrons = []
        lewis_valence = []
        lewis_bonding_target = []
        lewis_bonds_made = []
        lewis_adj_mat = []
        lewis_identical_mat = []

        # Determine the atoms with lone pairs that are unsatisfied as candidates for electron removal/addition to satisfy the total charge condition
        happy = [i[0] for i in bonding_pref if i[1] <= bonding_electrons[i[0]]]
        bonding_pref_ind = [i[0] for i in bonding_pref]

        # Determine is electrons need to be removed or added
        if q_tot > 0:
            adjust = -1
            octet_violate_e = []
            for count_j, j in enumerate(self.elements):
                if j.lower() in ["c", "n", "o", "f", "si", "p", "s", "cl"] and count_j not in bonding_pref_ind:
                    if bonding_electrons[count_j] * 2 + lone_electrons[count_j] > 8:
                        octet_violate_e += [count_j]
                elif j.lower() in ["br", "i"] and count_j not in bonding_pref_ind:
                    if bonding_electrons[count_j] * 2 + lone_electrons[count_j] > 18:
                        octet_violate_e += [count_j]

            normal_adjust = [count_i for count_i, i in enumerate(lone_electrons) if
                             i > 0 and count_i not in happy and count_i not in octet_violate_e]

        elif q_tot < 0:
            adjust = 1
            octet_violate_e = []
            for count_j, j in enumerate(self.elements):
                if j.lower() in ["c", "n", "o", "f", "si", "p", "s", "cl"] and count_j not in bonding_pref_ind:
                    if bonding_electrons[count_j] * 2 + lone_electrons[count_j] < 8:
                        octet_violate_e += [count_j]

                elif j.lower() in ["br", "i"] and count_j not in bonding_pref_ind:
                    if bonding_electrons[count_j] * 2 + lone_electrons[count_j] < 18:
                        octet_violate_e += [count_j]

            normal_adjust = [count_i for count_i, i in enumerate(lone_electrons) if
                             i > 0 and count_i not in happy and count_i not in octet_violate_e]

        else:
            adjust = 1
            octet_violate_e = []
            normal_adjust = [count_i for count_i, i in enumerate(lone_electrons) if i > 0 and count_i not in happy]

        # The outer loop checks each bonding structure produced by the inner loop for consistency with
        # the user specified "pref_bonding" and pref_argument with bonding electrons are
        for dummy_counter in range(lewis_total):

            lewis_loop_list = loop_list
            random.shuffle(lewis_loop_list)
            outer_counter = 0
            inner_max_cycles = 1000
            outer_max_cycles = 1000
            bond_sat = False

            lewis_lone_electrons.append(deepcopy(lone_electrons))
            lewis_bonding_electrons.append(deepcopy(bonding_electrons))
            lewis_core_electrons.append(deepcopy(core_electrons))
            lewis_valence.append(deepcopy(valence))
            lewis_bonding_target.append(deepcopy(bonding_target))
            lewis_bonds_made.append(deepcopy(bonds_made))
            lewis_adj_mat.append(deepcopy(adj_mat))
            lewis_counter = len(lewis_lone_electrons) - 1

            # Adjust the number of electrons by removing or adding to the available lone pairs
            # The algorithm simply adds/removes from the first N lone pairs that are discovered
            random.shuffle(octet_violate_e)
            random.shuffle(normal_adjust)
            adjust_ind = octet_violate_e + normal_adjust

            if len(adjust_ind) >= abs(q_tot):
                for i in range(abs(q_tot)):
                    lewis_lone_electrons[-1][adjust_ind[i]] += adjust
                    lewis_bonding_target[-1][adjust_ind[i]] += adjust
            else:
                for i in range(abs(q_tot)):
                    lewis_lone_electrons[-1][0] += adjust
                    lewis_bonding_target[-1][0] += adjust

            # Search for an optimal lewis structure
            while bond_sat is False:

                # Initialize necessary objects
                change_list = range(len(lewis_lone_electrons[lewis_counter]))
                inner_counter = 0
                bond_sat = True
                # Inner loop forms bonds to remove radicals or underbonded atoms until no further
                # changes in the bonding pattern are observed.
                while len(change_list) > 0:
                    change_list = []
                    for i in lewis_loop_list:

                        # List of atoms that already have a satisfactory binding configuration.
                        happy = [j[0] for j in bonding_pref if j[1] <= lewis_bonding_electrons[lewis_counter][j[0]]]

                        # If the current atom already has its target configuration then no further action is taken
                        if i in happy: continue

                        # If there are no lone electrons or too more bond formed then skip
                        if lewis_lone_electrons[lewis_counter][i] == 0: continue

                        # Take action if this atom has a radical or an unsatifisied bonding condition
                        if lewis_lone_electrons[lewis_counter][i] % 2 != 0 or lewis_bonding_electrons[lewis_counter][
                            i] != lewis_bonding_target[lewis_counter][i]:
                            # Try to form a bond with a neighboring radical (valence +1/-1 check ensures that no improper 5-bonded atoms are formed)
                            lewis_bonded_radicals = [(-find_lewis.en[self.elements[count_j].lower()], count_j) for count_j, j
                                                     in enumerate(self.adj_mat[i]) if
                                                     j == 1 and lewis_lone_electrons[lewis_counter][count_j] % 2 != 0 \
                                                     and 2 * (lewis_bonding_electrons[lewis_counter][count_j] + 1) + (
                                                                 lewis_lone_electrons[lewis_counter][count_j] - 1) <=
                                                     lewis_valence[lewis_counter][count_j] \
                                                     and lewis_lone_electrons[lewis_counter][
                                                         count_j] - 1 >= 0 and count_j not in happy]

                            lewis_bonded_lonepairs = [(-find_lewis.en[self.elements[count_j].lower()], count_j) for
                                                      count_j, j in enumerate(self.adj_mat[i]) if
                                                      j == 1 and lewis_lone_electrons[lewis_counter][count_j] > 0 \
                                                      and 2 * (lewis_bonding_electrons[lewis_counter][count_j] + 1) + (
                                                                  lewis_lone_electrons[lewis_counter][count_j] - 1) <=
                                                      lewis_valence[lewis_counter][count_j] and
                                                      lewis_lone_electrons[lewis_counter][count_j] - 1 >= 0 \
                                                      and count_j not in happy]

                            # Sort by atomic number (cheap way of sorting carbon before other atoms, should probably switch over to electronegativities)
                            lewis_bonded_radicals = [j[1] for j in sorted(lewis_bonded_radicals)]
                            lewis_bonded_lonepairs = [j[1] for j in sorted(lewis_bonded_lonepairs)]

                            # Correcting radicals is attempted first
                            if len(lewis_bonded_radicals) > 0:
                                lewis_bonding_electrons[lewis_counter][i] += 1
                                lewis_bonding_electrons[lewis_counter][lewis_bonded_radicals[0]] += 1
                                lewis_adj_mat[lewis_counter][i][lewis_bonded_radicals[0]] += 1
                                lewis_adj_mat[lewis_counter][lewis_bonded_radicals[0]][i] += 1
                                lewis_lone_electrons[lewis_counter][i] -= 1
                                lewis_lone_electrons[lewis_counter][lewis_bonded_radicals[0]] -= 1
                                change_list += [i, lewis_bonded_radicals[0]]
                                lewis_bonds_made[lewis_counter] += [(i, lewis_bonded_radicals[0])]

                            # Else try to form a bond with a neighboring atom with spare lone electrons (valence check ensures that no improper 5-bonded atoms are formed)
                            elif len(lewis_bonded_lonepairs) > 0:
                                lewis_bonding_electrons[lewis_counter][i] += 1
                                lewis_bonding_electrons[lewis_counter][lewis_bonded_lonepairs[0]] += 1
                                lewis_adj_mat[lewis_counter][i][lewis_bonded_lonepairs[0]] += 1
                                lewis_adj_mat[lewis_counter][lewis_bonded_lonepairs[0]][i] += 1
                                lewis_lone_electrons[lewis_counter][i] -= 1
                                lewis_lone_electrons[lewis_counter][lewis_bonded_lonepairs[0]] -= 1
                                change_list += [i, lewis_bonded_lonepairs[0]]
                                lewis_bonds_made[lewis_counter] += [(i, lewis_bonded_lonepairs[0])]
                                # lewis_bonds_en[lewis_counter] += 1.0/find_lewis.en[elements[i].lower()]/find_lewis.en[elements[lewis_bonded_lonepairs[0]].lower()]

                    # Increment the counter and break if the maximum number of attempts have been made
                    inner_counter += 1
                    if inner_counter >= inner_max_cycles:
                        print(
                            "WARNING: maximum attempts to establish a reasonable lewis-structure exceeded ({}).".format(
                                inner_max_cycles))

                # Check if the user specified preferred bond order has been achieved.
                if bonding_pref is not None:
                    unhappy = [i[0] for i in bonding_pref if i[1] != lewis_bonding_electrons[lewis_counter][i[0]]]
                    if len(unhappy) > 0:

                        # Break the first bond involving one of the atoms bonded to the under/over coordinated atoms
                        ind = set([unhappy[0]] + [count_i for count_i, i in enumerate(self.adj_mat[unhappy[0]]) if
                                                  i == 1 and (count_i, unhappy[0]) not in off_limits])

                        # Check if a rearrangment is possible, break if none are available
                        try:
                            break_bond = next(i for i in lewis_bonds_made[lewis_counter] if i[0] in ind or i[1] in ind)
                        except:
                            print(
                                "WARNING: no further bond rearrangments are possible and bonding_pref is still not satisfied.")
                            break

                        # Perform bond rearrangment
                        lewis_bonding_electrons[lewis_counter][break_bond[0]] -= 1
                        lewis_lone_electrons[lewis_counter][break_bond[0]] += 1
                        lewis_adj_mat[lewis_counter][break_bond[0]][break_bond[1]] -= 1
                        lewis_adj_mat[lewis_counter][break_bond[1]][break_bond[0]] -= 1
                        lewis_bonding_electrons[lewis_counter][break_bond[1]] -= 1
                        lewis_lone_electrons[lewis_counter][break_bond[1]] += 1

                        # Remove the bond from the list and reorder lewis_loop_list so that the indices involved in the bond are put last
                        lewis_bonds_made[lewis_counter].remove(break_bond)
                        lewis_loop_list.remove(break_bond[0])
                        lewis_loop_list.remove(break_bond[1])
                        lewis_loop_list += [break_bond[0], break_bond[1]]

                        # Update the bond_sat flag
                        bond_sat = False

                    # Increment the counter and break if the maximum number of attempts have been made
                    outer_counter += 1

                    # Periodically reorder the list to avoid some cyclical walks
                    if outer_counter % 100 == 0:
                        lewis_loop_list = self.reorder_list(lewis_loop_list, atomic_number)

                    # Print diagnostic upon failure
                    if outer_counter >= outer_max_cycles:
                        print("WARNING: maximum attempts to establish a lewis-structure consistent")
                        print("         with the user supplied bonding preference has been exceeded ({}).".format(
                            outer_max_cycles))
                        break

            # Re-apply keep_lone: remove one electron from such index
            for count_i in keep_lone:
                lewis_lone_electrons[lewis_counter][count_i] -= 1

            # Special cases, share pair of electrons
            total_electron = np.array(lewis_lone_electrons[lewis_counter]) + np.array(
                lewis_bonding_electrons[lewis_counter]) * 2

            # count for atom which doesn't satisfy
            # Notice: need systematical check for this part !!!
            unsatisfy = [count_t for count_t, te in enumerate(total_electron) if
                         te > 2 and te < 8 and te % 2 == 0 and self.elements[count_t].lower() in ["c", "n", "o", "f", "si",
                                                                                             "p", "s", "cl"]]
            for uns in unsatisfy:
                full_connect = [count_i for count_i, i in enumerate(self.adj_mat[uns]) if
                                i == 1 and total_electron[count_i] == 8 and lewis_lone_electrons[lewis_counter][
                                    count_i] >= 2]
                if len(full_connect) > 0:
                    lewis_lone_electrons[lewis_counter][full_connect[0]] -= 2
                    lewis_bonding_electrons[lewis_counter][uns] += 1
                    lewis_bonding_electrons[lewis_counter][full_connect[0]] += 1
                    lewis_adj_mat[lewis_counter][uns][full_connect[0]] += 1
                    lewis_adj_mat[lewis_counter][full_connect[0]][uns] += 1

                    # Delete last entry in the lewis np.arrays if the electronic structure is not unique: introduce identical_mat includes both info of bond_mats and formal_charges
            identical_mat = np.vstack([lewis_adj_mat[-1], np.array(
                [valence_list[k] - lewis_bonding_electrons[-1][k] - lewis_lone_electrons[-1][k] for k in
                 range(len(self.elements))])])
            lewis_identical_mat.append(identical_mat)

            if self.array_unique(lewis_identical_mat[-1], lewis_identical_mat[:-1]) is False:
                lewis_lone_electrons = lewis_lone_electrons[:-1]
                lewis_bonding_electrons = lewis_bonding_electrons[:-1]
                lewis_core_electrons = lewis_core_electrons[:-1]
                lewis_valence = lewis_valence[:-1]
                lewis_bonding_target = lewis_bonding_target[:-1]
                lewis_bonds_made = lewis_bonds_made[:-1]
                lewis_adj_mat = lewis_adj_mat[:-1]
                lewis_identical_mat = lewis_identical_mat[:-1]

        # Find the total number of lone electrons in each structure
        lone_electrons_sums = []
        for i in range(len(lewis_lone_electrons)):
            lone_electrons_sums.append(sum(lewis_lone_electrons[i]))

        # Find octet violations in each structure
        octet_violations = []
        for i in range(len(lewis_lone_electrons)):
            ov = 0
            if octet_opt is True:
                for count_j, j in enumerate(self.elements):
                    if j.lower() in ["c", "n", "o", "f", "si", "p", "s", "cl", "br",
                                     "i"] and count_j not in bonding_pref_ind:
                        if lewis_bonding_electrons[i][count_j] * 2 + lewis_lone_electrons[i][count_j] != 8 and \
                                lewis_bonding_electrons[i][count_j] * 2 + lewis_lone_electrons[i][count_j] != 18:
                            ov += 1
            octet_violations.append(ov)

        ## Calculate bonding energy
        lewis_bonds_energy = []
        for bonds_made in lewis_bonds_made:
            for lb, bond_made in enumerate(bonds_made): bonds_made[lb] = tuple(sorted(bond_made))
            count_bonds_made = ["{}-{}-{}".format(min(atomic_number[bm[0]], atomic_number[bm[1]]),
                                                  max(atomic_number[bm[0]], atomic_number[bm[1]]), bonds_made.count(bm))
                                for bm in set(bonds_made)]
            lewis_bonds_energy += [
                sum([find_lewis.be[cbm] if cbm in find_lewis.be.keys() else -10000.0 for cbm in count_bonds_made])]

        # normalize the effect
        lewis_bonds_energy = [-be / max(1, max(lewis_bonds_energy)) for be in lewis_bonds_energy]

        ## Find the total formal charge for each structure
        formal_charges_sums = []
        for i in range(len(lewis_lone_electrons)):
            fc = 0
            for j in range(len(self.elements)):
                fc += abs(valence_list[j] - lewis_bonding_electrons[i][j] - lewis_lone_electrons[i][j])
            formal_charges_sums.append(fc)

        ## Find formal charge eletronegativity contribution
        lewis_formal_charge = [
            [valence_list[i] - lewis_bonding_electrons[_][i] - lewis_lone_electrons[_][i] for i in range(len(self.elements))]
            for _ in range(len(lewis_lone_electrons))]
        lewis_keep_lone = [[count_i for count_i, i in enumerate(lone) if i % 2 != 0] for lone in lewis_lone_electrons]
        lewis_fc_aro = []  # form aromatic rings to stablize charge/radical
        lewis_fc_en = []  # Electronegativity for stabling charge/radical
        lewis_fc_pol = []  # Polarizability for stabling charge/radical
        lewis_fc_hc = []  # Hyper-conjugation contribution
        lewis_fc_all = []  # Avoid forming allene

        for i in range(len(lewis_lone_electrons)):

            formal_charge = lewis_formal_charge[i]
            radical_atom = lewis_keep_lone[i]
            fc_ind = [(count_j, j) for count_j, j in enumerate(formal_charge) if j != 0]

            for R_ind in radical_atom:  # assign +0.5 for radical
                fc_ind += [(R_ind, 0.5)]

            # initialize en,pol and hc
            fc_aro, fc_en, fc_pol, fc_hc, fc_all = 0, 0, 0, 0, 0
            bond_madesi = lewis_bonds_made[i]

            # consider fc_aro if a ring exists
            if len(rings) > 0 and len(fc_ind) > 0:

                # loop over all rings
                for ring in rings:

                    # if radical/charged atoms in this ring, and this ring consists of even number atoms
                    # if len(ring) % 2 == 0 and (True in [ind[0] in ring for ind in fc_ind]):
                    if len(ring) % 2 == 0:

                        bond_madesi_ring = [bmade for bmade in bond_madesi if (bmade[0] in ring and bmade[1] in ring)]

                        # check whether generate an aromatic ring
                        if len(bond_madesi_ring) == len(ring) / 2 and len(
                            set(sum([list(bond) for bond in bond_madesi_ring], []))) == len(ring): fc_aro += -1

            # consider whether allene will form
            # check whether two bonds form on (a,b) and (a,c)
            if len(bond_madesi) >= 2:
                # find whether an atom forms two bonds
                atoms = list(sum(bond_madesi, ()))
                common_atom = [ind for ind in set(atoms) if atoms.count(ind) > 1]
                if len(common_atom) > 0:
                    for atom_c in common_atom:
                        # find which bond formed involves this common atom
                        which_bond = [ind for ind in bond_madesi if atom_c in ind]
                        bonded_atoms = list(sum(which_bond, ()))
                        # if a triple bond forms, there will be two atoms, if not, if will be allene
                        if len(set(bonded_atoms)) > 2: fc_all += 1

            # Loop over formal charges and radicals
            for count_fc in fc_ind:

                ind = count_fc[0]
                charge = count_fc[1]
                # Count the self contribution: (-) on the most electronegative atom and (+) on the least electronegative atom
                fc_en += 10 * charge * find_lewis.en[self.elements[ind].lower()]

                # Find the nearest and next-nearest atoms for each formal_charge/radical contained atom
                gs = self.graph_seps()
                nearest_atoms = [count_k for count_k, k in enumerate(lewis_adj_mat[i][ind]) if k >= 1]
                NN_atoms = list(set([count_j for count_j, j in enumerate(gs[ind]) if j == 2]))

                # only count when en > en(C)
                fc_en += charge * (sum([find_lewis.en[self.elements[count_k].lower()] for count_k in nearest_atoms if
                                        find_lewis.en[self.elements[count_k].lower()] > 2.54]) + \
                                   sum([find_lewis.en[self.elements[count_k].lower()] for count_k in NN_atoms if
                                        find_lewis.en[self.elements[count_k].lower()] > 2.54]) * 0.1)

                if charge < 0:  # Polarizability only affects negative charge
                    fc_pol += charge * sum([find_lewis.pol[self.elements[count_k].lower()] for count_k in nearest_atoms])

                # find hyper-conjugation strcuture
                nearby_carbon = [nind for nind in nearest_atoms if self.elements[nind].lower() == 'c']
                for carbon_ind in nearby_carbon:
                    carbon_nearby = [nind for nind in NN_atoms if
                                     lewis_adj_mat[i][carbon_ind][nind] >= 1 and self.elements[nind].lower() in ['c', 'h']]
                    if len(carbon_nearby) == 3:
                        fc_hc -= charge * (
                                    len([nind for nind in carbon_nearby if self.elements[nind].lower() == 'c']) * 2 + len(
                                [nind for nind in carbon_nearby if self.elements[nind].lower() == 'h']))

            lewis_fc_aro.append(fc_aro)
            lewis_fc_all.append(fc_all)
            lewis_fc_en.append(fc_en)
            lewis_fc_pol.append(fc_pol)
            lewis_fc_hc.append(fc_hc)

            # normalize the effect
        lewis_fc_en = [lfc / max(1, max(abs(np.array(lewis_fc_en)))) for lfc in lewis_fc_en]
        lewis_fc_pol = [lfp / max(1, max(abs(np.array(lewis_fc_pol)))) for lfp in lewis_fc_pol]

        # Add the total number of radicals to the total formal charge to determine the criteria.
        # The radical count is scaled by 0.01 and the lone pair count is scaled by 0.001. This results
        # in the structure with the lowest formal charge always being returned, and the radical count
        # only being considered if structures with equivalent formal charges are found, and likewise with
        # the lone pair count. The structure(s) with the lowest score will be returned.
        lewis_criteria = []
        for i in range(len(lewis_lone_electrons)):
            lewis_criteria.append(10.0 * octet_violations[i] + abs(formal_charges_sums[i]) + 0.1 * sum(
                [1 for j in lewis_lone_electrons[i] if j % 2 != 0]) + \
                                  0.05 * lewis_fc_aro[i] + 0.01 * lewis_fc_all[i] + 0.005 * lewis_fc_en[i] + 0.001 *
                                  lewis_fc_pol[i] + 0.0001 * lewis_fc_hc[i] + \
                                  0.0001 * lewis_bonds_energy[i])
        # print(lewis_criteria)
        # print(lewis_bonds_made)
        # print(octet_violations,formal_charges_sums)
        # print(lewis_fc_aro,lewis_fc_all,lewis_fc_pol,lewis_fc_hc,lewis_bonds_energy)
        # exit()
        # sort from least to most and return a list containing the origial list's indices in the correct order
        best_lewis = [i[0] for i in sorted(enumerate(lewis_criteria), key=lambda x: x[1])]
        best_lewis = [i for i in best_lewis if lewis_criteria[i] == lewis_criteria[best_lewis[0]]]

        # Finally check formal charge to keep those with
        lewis_re_fc = [lewis_formal_charge[_] + lewis_keep_lone[_] for _ in best_lewis]
        appear_times = [lewis_re_fc.count(i) for i in lewis_re_fc]
        best_lewis = [best_lewis[i] for i in range(len(lewis_re_fc)) if appear_times[i] == max(appear_times)]

        # Apply keep_lone information, remove the electron to form lone electron
        for i in best_lewis:
            for j in keep_lone:
                lewis_lone_electrons[i][j] -= 1

        # Print diagnostics
        if verbose is True:
            for i in best_lewis:
                print("Bonding Matrix  {}".format(i))
                print("Formal_charge:  {}".format(formal_charges_sums[i]))
                print("Lewis_criteria: {}\n".format(lewis_criteria[i]))
                print("{:<40s} {:<40s} {:<15s} {:<15s}".format("Elements", "Bond_Mat", "Lone_Electrons", "FC"))
                for j in range(len(self.elements)):
                    print("{:<40s} {}    {} {}".format(self.elements[j], " ".join([str(k) for k in lewis_adj_mat[i][j]]),
                                                       lewis_lone_electrons[i][j],
                                                       valence_list[j] - lewis_bonding_electrons[i][j] -
                                                       lewis_lone_electrons[i][j]))
                print(" ")

        # If only the bonding matrix is requested, then only that is returned
        if b_mat_only is True:
            if return_FC is False:
                return [lewis_adj_mat[_] for _ in best_lewis]
            else:
                return [lewis_adj_mat[_] for _ in best_lewis], [lewis_formal_charge[_] for _ in best_lewis]

        # return like check_lewis function
        if check_lewis_flag is True:
            if return_pref is True:
                return lewis_lone_electrons[best_lewis[0]], lewis_bonding_electrons[best_lewis[0]], \
                lewis_core_electrons[best_lewis[0]], bonding_pref
            else:
                return lewis_lone_electrons[best_lewis[0]], lewis_bonding_electrons[best_lewis[0]], \
                lewis_core_electrons[best_lewis[0]]

        # Optional bonding pref return to handle cases with special groups
        if return_pref is True:
            if return_FC is False:
                return [lewis_lone_electrons[_] for _ in best_lewis], [lewis_bonding_electrons[_] for _ in best_lewis], \
                    [lewis_core_electrons[_] for _ in best_lewis], [lewis_adj_mat[_] for _ in best_lewis], bonding_pref
            else:
                return [lewis_lone_electrons[_] for _ in best_lewis], [lewis_bonding_electrons[_] for _ in best_lewis], \
                    [lewis_core_electrons[_] for _ in best_lewis], [lewis_adj_mat[_] for _ in best_lewis], [
                    lewis_formal_charge[_] for _ in best_lewis], bonding_pref

        else:
            if return_FC is False:
                return [lewis_lone_electrons[_] for _ in best_lewis], [lewis_bonding_electrons[_] for _ in best_lewis], \
                    [lewis_core_electrons[_] for _ in best_lewis], [lewis_adj_mat[_] for _ in best_lewis]
            else:
                return [lewis_lone_electrons[_] for _ in best_lewis], [lewis_bonding_electrons[_] for _ in best_lewis], \
                    [lewis_core_electrons[_] for _ in best_lewis], [lewis_adj_mat[_] for _ in best_lewis], [
                    lewis_formal_charge[_] for _ in best_lewis]
    @staticmethod
    def reorder_list(loop_list, atomic_number):
        """
        Helper function to check_lewis and get_bonds that rolls the loop_list carbon elements
        """
        c_types = [count_i for count_i, i in enumerate(loop_list) if atomic_number[i] == 6]
        others = [count_i for count_i, i in enumerate(loop_list) if atomic_number[i] != 6]
        if len(c_types) > 1:
            c_types = c_types + [c_types.pop(0)]
        return [loop_list[i] for i in c_types + others]

    @staticmethod
    def array_unique(a, a_list):
        """
        Description: Checks if an np.array "a" is unique compared with a list of np.arrays "a_list"
        at the first match False is returned.
        """
        for i in a_list:
            if np.array_equal(a, i):
                return False
        return True

    # Return bool depending on if the atom is a nitro nitrogen atom
    def is_nitro(self, i, alt_adj_mat=None):
        if alt_adj_mat is not None:
            adj_mat = alt_adj_mat
        else:
            adj_mat = self.adj_mat

        status = False
        if self.elements[i] not in ["N", "n"]:
            return False
        O_ind = [count_j for count_j, j in enumerate(adj_mat[i]) if
                 j == 1 and self.elements[count_j] in ["o", "O"] and sum(adj_mat[count_j]) == 1]
        if len(O_ind) == 2:
            return True
        else:
            return False

    # Return bool depending on if the atom is a sulfoxide sulfur atom
    def is_sulfoxide(self, i, alt_adj_mat=None):
        if alt_adj_mat is not None:
            adj_mat = alt_adj_mat
        else:
            adj_mat = self.adj_mat

        status = False
        if self.elements[i] not in ["S", "s"]:
            return False
        O_ind = [count_j for count_j, j in enumerate(adj_mat[i]) if
                 j == 1 and self.elements[count_j] in ["o", "O"] and sum(adj_mat[count_j]) == 1]
        connect = sum(adj_mat[i])

        if len(O_ind) == 1 and int(connect) == 3:
            return True
        else:
            return False

    # Return bool depending on if the atom is a sulfonyl sulfur atom
    def is_sulfonyl(self, i, alt_adj_mat=None):
        if alt_adj_mat is not None:
            adj_mat = alt_adj_mat
        else:
            adj_mat = self.adj_mat

        status = False
        if self.elements[i] not in ["S", "s"]:
            return False
        O_ind = [count_j for count_j, j in enumerate(adj_mat[i]) if
                 j == 1 and self.elements[count_j] in ["o", "O"] and sum(adj_mat[count_j]) == 1]
        connect = sum(adj_mat[i])

        if len(O_ind) == 2 and int(connect) == 4:
            return True
        else:
            return False

    # Return bool depending on if the atom is a phosphate phosphorus atom
    def is_phosphate(self, i, alt_adj_mat=None):
        if alt_adj_mat is not None:
            adj_mat = alt_adj_mat
        else:
            adj_mat = self.adj_mat

        status = False
        if self.elements[i] not in ["P", "p"]:
            return False
        O_ind = [count_j for count_j, j in enumerate(adj_mat[i]) if j == 1 and self.elements[count_j] in ["o", "O"]]
        O_ind_term = [j for j in O_ind if sum(adj_mat[j]) == 1]
        if len(O_ind) == 4 and sum(adj_mat[i]) == 4 and len(O_ind_term) > 0:
            return True
        else:
            return False

    # Return bool depending on if the atom is a cyano nitrogen atom
    def is_cyano(self, i, alt_adj_mat=None):
        if alt_adj_mat is not None:
            adj_mat = alt_adj_mat
        else:
            adj_mat = self.adj_mat

        status = False
        if self.elements[i] not in ["N", "n"] or sum(adj_mat[i]) > 1:
            return False

        C_ind = [count_j for count_j, j in enumerate(adj_mat[i]) if
                 j == 1 and self.elements[count_j] in ["c", "C"] and sum(adj_mat[count_j]) == 2]

        if len(C_ind) == 1:
            return True

        else:
            return False

    # Return bool depending on if the atom is a cyano nitrogen atom
    def is_isocyano(self, i, alt_adj_mat=None):
        if alt_adj_mat is not None:
            adj_mat = alt_adj_mat
        else:
            adj_mat = self.adj_mat

        if len(self.elements) <= 2:
            return False

        status = False
        if self.elements[i] not in ["N", "n"] or sum(adj_mat[i]) > 1:
            return False

        C_ind = [count_j for count_j, j in enumerate(adj_mat[i]) if
                 j == 1 and self.elements[count_j] in ["c", "C"] and sum(adj_mat[count_j]) == 1]
        if len(C_ind) == 1:
            return True
        else:
            return False

    # is function for fragments
    # Return bool depending on if the atom is a sulfoxide sulfur atom
    def is_frag_sulfoxide(self, i, alt_adj_mat=None):
        if alt_adj_mat is not None:
            adj_mat = alt_adj_mat
        else:
            adj_mat = self.adj_mat

        status = False
        if self.elements[i] not in ["S", "s"]:
            return False
        O_ind = [count_j for count_j, j in enumerate(adj_mat[i]) if
                 j == 1 and self.elements[count_j] in ["o", "O"] and sum(adj_mat[count_j]) == 1]
        connect = sum(adj_mat[i])

        if len(O_ind) >= 1 and int(connect) == 3:
            return True
        else:
            return False

    # Return bool depending on if the atom is a sulfonyl sulfur atom
    def is_frag_sulfonyl(self, i, alt_adj_mat=None):
        if alt_adj_mat is not None:
            adj_mat = alt_adj_mat
        else:
            adj_mat = self.adj_mat

        status = False
        if self.elements[i] not in ["S", "s"]:
            return False
        O_ind = [count_j for count_j, j in enumerate(adj_mat[i]) if
                 j == 1 and self.elements[count_j] in ["o", "O"] and sum(adj_mat[count_j]) == 1]
        connect = sum(adj_mat[i])

        if len(O_ind) >= 2 and int(connect) == 4:
            return True
        else:
            return False

    # Return bool depending on if the atom is a sulfonyl sulfur atom
    def is_frag_ethenone(self, i, alt_adj_mat=None):
        if alt_adj_mat is not None:
            adj_mat = alt_adj_mat
        else:
            adj_mat = self.adj_mat

        status = False
        if self.elements[i] not in ["C", "c"]:
            return False

        OS_ind = [count_j for count_j, j in enumerate(adj_mat[i]) if
                  j == 1 and self.elements[count_j] in ["o", "O", "s", "S"] and sum(adj_mat[count_j]) == 1]
        CN_ind = [count_j for count_j, j in enumerate(adj_mat[i]) if
                  j == 1 and self.elements[count_j] in ["c", "C", "n", "N"]]
        connect = sum(adj_mat[i])

        if len(OS_ind) == 1 and len(CN_ind) == 1 and int(connect) == 2:
            return True
        else:
            return False



    def return_ring_atom(self, idx, start=None, ring_size=10, counter=0, avoid_set=None, in_ring=None, alt_adj_mat=None):
        """
        Return true if idx is a ring atom
        sometime you'd want 
        """
        if alt_adj_mat is not None:
            adj_mat = alt_adj_mat
        else:
            adj_mat = self.adj_mat

        # Consistency/Termination checks
        if ring_size < 3:
            raise AdjacencyException("ERROR in ring_atom: ring_size variable must be set to an integer greater than 2!")

        if counter == ring_size:
            return False, []

        # Automatically assign start to the supplied idx value. For recursive calls this is set manually
        if start is None:
            start = idx

        if avoid_set is None:
            avoid_set = set([])

        if in_ring is None:
            in_ring = set([idx])

        # Trick: The fact that the smallest possible ring has three nodes can be used to simplify
        #        the algorithm by including the origin in avoid_set until after the second step
        if counter >= 2 and start in avoid_set:
            avoid_set.remove(start)

        elif counter < 2 and start not in avoid_set:
            avoid_set.add(start)

        # Update the avoid_set with the current idx value
        avoid_set.add(idx)

        # Loop over connections and recursively search for idx
        status = 0
        cons = [count_i for count_i, i in enumerate(adj_mat[idx]) if i == 1 and count_i not in avoid_set]
        # print cons,counter,start,avoid_set
        if len(cons) == 0:
            return False, []

        elif start in cons:
            return True, in_ring

        else:
            for i in cons:
                if \
                self.return_ring_atom(i, start=start, ring_size=ring_size, counter=counter + 1, avoid_set=avoid_set,
                                 in_ring=in_ring, alt_adj_mat=alt_adj_mat)[0] == True:
                    in_ring.add(i)
                    return True, in_ring
            return False, []

    def transify(self, start=-1, end=-1, elements=None, opt_terminals=False, opt_final=False):
        """
        This function takes a geometry and adjacency matrix and returns an "all-trans" aligned geometric conformer
        """

        # Sanity checks
        if start < -1 or start >= len(self.geometry): raise AdjacencyException(
            "ERROR in transify: start must be set to an index within geo or -1. Exiting...")
        if end < -1 or end >= len(self.geometry): raise AdjacencyException(
            "ERROR in transify: end must be set to an index within geo or -1. Exiting...")
        if len(self.geometry) != len(self.adj_mat): raise AdjacencyException(
            "ERROR in transify: geo and adj_mat must have the same dimensions. Exiting...")
        if opt_terminals == True and elements is None: raise AdjacencyException(
            "ERROR in transify: element information must be supplied for terminals to be optimized. Exiting...")
        if opt_terminals == True and elements is not None and len(elements) != len(
            self.geometry): raise AdjacencyException(
            "ERROR in transify: length of elements is not equal to the length of geo. Exiting...")
        if opt_final == True and elements is None: raise AdjacencyException(
            "ERROR in transify: element information must be supplied for the final structure to be optimized. Exiting...")
        if opt_final == True and elements is not None and len(elements) != len(self.geometry): raise AdjacencyException(
            "ERROR in transify: length of elements is not equal to the length of geo. Exiting...")

        # Return the geometry if it is empty
        if len(self.geometry) == 0:
            return self.geometry

        # Calculate the pair-wise graphical separations between all atoms
        seps = self.graph_seps()

        # Find the most graphically separated pair of atoms.
        if start == -1 and end == -1:
            max_ind = np.where(seps == seps.max())
            start, end = max_ind[0][0], max_ind[1][0]

        # Find the most graphically separated atom from the start atom
        elif end == -1:
            end = np.argmax(seps[start])

            # Find the shortest pathway between these points
        pathway = Dijkstra(self.adj_mat, start, end)

        # If the molecule doesn't have any dihedrals then return
        if len(pathway) < 4:
            return self.geometry

        # Initialize the list of terminal atoms and ring atoms (used in a couple places so initialized here relatively early)
        terminals = set([count_i for count_i, i in enumerate(self.adj_mat) if np.sum(i) == 1])
        rings = set([count_i for count_i, i in enumerate(self.geometry) if self.ring_atom(self.adj_mat, count_i) == True])

        # Work through the backbone dihedrals and straighten them out
        for i in range(1, len(pathway) - 2):

            if pathway[i] in rings: continue

            # Collect the atoms that are connected to the 2 atom of the dihedral but not to the 3 atom of the dihedral (and vice versus)
            group_1 = return_connected(adj_mat, start=pathway[i], avoid=[
                pathway[i + 1]])  # not used since only the forward portion of the pathway gets rotated
            group_2 = return_connected(adj_mat, start=pathway[i + 1], avoid=[pathway[i]])

            # Skip if the two groups are equal (happens in the case of rings)
            if group_1 == group_2: continue

            # Calculate the rotation vector and angle
            rot_vec = geo[pathway[i + 1]] - geo[pathway[i]]
            stationary = geo[pathway[i]]
            theta = np.pi - dihedral_calc(geo[pathway[i - 1:i + 3]])

            # Perform the rotation
            for j in group_2:
                geo[j] = axis_rot(geo[j], rot_vec, stationary, theta, mode='radian')

            # Check for overlaps
            theta = np.pi - dihedral_calc(geo[pathway[i - 1:i + 3]])

            # Consider including a check for overlaps and performing a gauche rotation as an alternative.
            # if len(np.where(cdist(geo,geo) < 1.0 & seps > 3)) > 0:

        # Identify the branch points
        branches = [count_i for count_i, i in enumerate(adj_mat) if len([j for count_j, j in enumerate(i) if
                                                                         j == 1 and count_j not in terminals and count_j not in rings]) > 2]
        avoid_list = terminals

        # Loop over the branch points and correct their dihedrals
        for count_i, i in enumerate(branches):

            branching_cons = [count_j for count_j, j in enumerate(adj_mat[i]) if j == 1 and count_j not in pathway]

            for c in branching_cons:

                # Find the connections to this branch point that are terminal atoms
                conn_terminals = [j for j in branching_cons if j != c]

                # If a terminal connection exists then it is used to orient the branch by prepending it to the branch index list
                if len(conn_terminals) > 0:
                    branch_ind = [conn_terminals[0]] + list(return_connected(adj_mat, start=i,
                                                                             avoid=[count_j for count_j, j in
                                                                                    enumerate(adj_mat[i]) if
                                                                                    j == 1 and count_j in pathway] + conn_terminals))
                    geo[branch_ind[1:]] = transify(geo[branch_ind], adj_mat[branch_ind, :][:, branch_ind], start=0)[1:]

                # If no terminal connections exists then the first dihedral of the branch is not adjusted. The logic is that the sp2/sp3 alignment is better judged by the initial guess.
                else:
                    branch_ind = list(return_connected(adj_mat, start=i,
                                                       avoid=[count_j for count_j, j in enumerate(adj_mat[i]) if
                                                              j == 1 and count_j in pathway] + conn_terminals))
                    geo[branch_ind] = transify(geo[branch_ind], adj_mat[branch_ind, :][:, branch_ind],
                                               start=next(count_j for count_j, j in enumerate(branch_ind) if j == i))

        # If optimize structure is set to True then the transified structure is relaxed using obminimize
        if opt_final == True:
            opt_geo(geo, adj_mat, elements, ff='mmff94')

        # If optimize terminals is set to true then a conformer search is performed over the terminal groups
        if opt_terminals == True:
            geo = opt_terminal_centers(geo, adj_mat, elements, ff='mmff94')

        return geo


