from model.adjacency import AdjacencyMatrix
from model.structure import StructureBase, StructureException
from utilities.PeriodictTable import Tables
from utilities.Excpetions import TAFFIException
import numpy as np
import random, sys
from copy import deepcopy
from itertools import chain, combinations
from types import SimpleNamespace


class AtomTypesException(TAFFIException):
    pass


def type_adjmat(label):
    """
    This function generates a local adjacency matrix from an atomtype label

    :param label: one taffi atomtype
    :type label: str

    :returns:
        *adj_mat* - local adjacency matrix
        *labels* - elements from that atom type
    """

    # Initialize breaks (indices for brackets in "label"), atoms (start:end index tuples for the atom labels),
    # labels ( basis of atom labels indexed to the adjmat), and adj_mat (the local adjacency matrix with connectivity relationships)
    breaks = [count_i for count_i, i in enumerate(label) if i in ['[', ']', '=', '#']]
    atoms = [(i + 1, i + breaks[count_i + 1] - i) for count_i, i in enumerate(breaks[:-1]) if
             breaks[count_i + 1] - i > 1]
    labels = [label[i[0]:i[1]] for i in atoms]
    adj_mat = np.zeros([len(atoms), len(atoms)])

    # Loop over atoms
    for count_i, i in enumerate(atoms):

        # Initialize variables
        starting_index = breaks.index(i[1])  # index of the nearest bracket forward from the current atom label
        break_count = 0  # counter for keeping track of who needs parsing and who is connected to who

        # Loop over brackets
        for count_j, j in enumerate(breaks[starting_index:-1]):

            # Increment break_count + 1 for "open" brackets
            if label[j] == "[": break_count += 1

            # Increment break_count - 1 for "closed" brackets
            if label[j] == "]": break_count -= 1

            # If the break_count variable is -1 then all connections have been found
            if break_count == -1:
                break

            # When break_count == 1 and the parser resides at an open bracket, the next atom past the bracket is connected
            if break_count == 1 and label[j] == "[":
                idx = next(count_k for count_k, k in enumerate(atoms) if k[0] == j + 1)
                adj_mat[count_i, idx] = 1
                adj_mat[idx, count_i] = 1
    return adj_mat, labels


class AtomTypesBaseOperation(StructureBase):
    """
    Atom types' base operation methods
    a lot of recursive functions, most of them can't be used without higher level wrapper
    """

    def __init__(self, structure_data={}, AdjMat=None, adj_mat=None):
        super().__init__()
        self._AdjMat = AdjacencyMatrix()

        if AdjMat is not None:
            self._AdjMat = AdjMat
            # this also parse the elements, geometry, q_tot from AdjMat to AtomTypes
            self.parse_data_from_AdjMat(AdjMat)
        if adj_mat is not None:
            self._AdjMat.adj_mat = adj_mat
        if structure_data != {}:
            # The sequence of parsing AdjMat first then structure data is important
            # the structure parsed here will also CHANGE the structure in self._AdjMat
            self.parse_data(**structure_data)

    @property
    def AdjMat(self):
        # this make sure it won't be empty when being called
        if self._AdjMat is None:
            raise AtomTypesException('AdjMat not set')
        return self._AdjMat

    @property
    def adj_mat(self):
        if self._AdjMat is None:
            raise AtomTypesException('AdjMat not set, cannot retrieve its adj_mat')
        return self._AdjMat.adj_mat

    @property
    def graph_seps(self):
        return self._AdjMat.graph_seps()

    @StructureBase.elements.setter
    def elements(self, value):
        self._AdjMat.elements = value
        super(AtomTypesBaseOperation, AtomTypesBaseOperation).elements.__set__(self, value)

    @StructureBase.geometry.setter
    def geometry(self, value):
        self._AdjMat.geometry = value
        super(AtomTypesBaseOperation, AtomTypesBaseOperation).geometry.__set__(self, value)

    @StructureBase.q_tot.setter
    def q_tot(self, value):
        self._AdjMat.q_tot = value
        super(AtomTypesBaseOperation, AtomTypesBaseOperation).q_tot.__set__(self, value)

    def parse_data_from_AdjMat(self, AdjMat):
        self._AdjMat = AdjMat
        self.elements = AdjMat.elements
        self.geometry = AdjMat.geometry
        self.q_tot = AdjMat.q_tot

    def build_adj_mat(self):
        """
        build only adj_mat, which means AdjMat should already be setup
        """
        if self._AdjMat is None:
            raise AtomTypesException("AdjMatrix obj within Atomtype not yet build")
        self._AdjMat.build_adj_mat()

    def taffi_type(self, ind, masses, gens=2, avoid=[], fc=[]):
        """
        recursive function
        adjacency matrix based algorithm for identifying the taffi atom type
        the gens change within, so need to be an argument instead of a class member

        :param ind: atom index
        :type ind: int

        :param masses: masses of the atoms
        :type masses: list

        :return: taffi_atom type
        """

        periodic = Tables.PERIODICT

        # Check fc condition
        if len(fc) == 0:
            fc = [''] * len(self.elements)
        if len(fc) != len(self.elements):
            raise AtomTypesException("ERROR in taffi_type: fc must have the same dimensions as elements and A")

        # Find connections, avoid is used to avoid backtracking
        cons = [count_i for count_i, i in enumerate(self.adj_mat[ind]) if i == 1 and count_i not in avoid]

        # Sort the connections based on the hash function
        if len(cons) > 0:
            cons = list(zip(*sorted([(self.atom_hash(i, masses, gens=gens - 1), i) for i in cons])[::-1]))[1]

            # Calculate the subbranches
        # NOTE: recursive call with the avoid list results
        if gens == 0:
            subs = []
        else:
            subs = [self.taffi_type(i, masses, gens=gens - 1, avoid=[ind], fc=fc) for i in cons]

        # Calculate formal charge terms
        return "{}".format(periodic[self.elements[ind].lower()]) + fc[ind] + "".join(
            ["[" + i + "]" for i in subs])

    def atom_hash(self, ind, M, alpha=100.0, beta=0.1, gens=10):
        """
        hashing function for canonicalizing geometries on the basis of their adjacency matrices and elements
        :type ind: int
        :param ind: index of the atom being hashed

        :type M: list
        :param M: masses of the atoms in the molecule

        """
        if gens <= 0:
            return self.rec_sum(ind, M, beta, gens=0)
        else:
            return alpha * sum(self.adj_mat[ind]) + self.rec_sum(ind, M, beta, gens)

    def rec_sum(self, ind, M, beta, gens, avoid_list=[]):
        """
        recursive function for summing up the masses at each generation of connections.

        :type ind: int
        :param ind: index of the atom being hashed

        :type M: list
        :param M: masses of the atoms in the molecule

        :type gens: int
        :param gens: taffi_type depth
        """
        if gens != 0:
            tmp = M[ind] * beta
            new = [count_j for count_j, j in enumerate(self.adj_mat[ind]) if j == 1 and count_j not in avoid_list]
            if len(new) > 0:
                for i in new:
                    tmp += self.rec_sum(i, M, beta * 0.1, gens - 1, avoid_list=avoid_list + [ind])
                return tmp
            else:
                return tmp
        else:
            return M[ind] * beta

    # list taffi type, not used anymore
    def taffi_type_list(self, ind, adj_list, masses, gens=2, avoid=[]):

        periodic = Tables.PERIODICT

        # Find connections, avoid is used to avoid backtracking
        cons = [i for i in adj_list[ind] if i not in avoid]

        # Sort the connections based on the hash function
        if len(cons) > 0:
            cons = \
                list(zip(*sorted([(self.atom_hash_list(i, adj_list, masses, gens=gens - 1), i) for i in cons])[::-1]))[
                    1]

        # Calculate the subbranches
        # NOTE: recursive call with the avoid list results
        if gens == 0:
            subs = []
        else:
            subs = [self.taffi_type_list(i, adj_list, masses, gens=gens - 1, avoid=[ind]) for i in cons]

        return "{}".format(periodic[self.elements[ind].lower()]) + "".join(["[" + i + "]" for i in subs])

    def atom_hash_list(self, ind, A, M, alpha=100.0, beta=0.1, gens=10):
        """
        hashing function for canonicalizing geometries on the basis of their adjacency lists and elements
        ind  : index of the atom being hashed
        M    : masses of the atoms in the molecule
        gens : depth of the search used for the hash
        """
        if gens <= 0:
            return self.rec_sum_list(ind, A, M, beta, gens=0)
        else:
            return alpha * len(A[ind]) + self.rec_sum_list(ind, A, M, beta, gens)

    def rec_sum_list(self, ind, A, M, beta, gens, avoid_list=[]):
        """
        recursive function for summing up the masses at each generation of connections.

        :param ind: atom index
        :type ind: int

        :param A: adjacency matrix
        :type A: np matrix

        :param M: masses
        :type M: list
        """
        if gens != 0:
            tmp = M[ind] * beta
            new = [j for j in A[ind] if j not in avoid_list]
            if len(new) > 0:
                for i in new:
                    tmp += self.rec_sum_list(i, A, M, beta * 0.1, gens - 1, avoid_list=avoid_list + [ind])
                return tmp
            else:
                return tmp
        else:
            return M[ind] * beta

    def ring_atom_new(self, idx, start=None, ring_size=10, counter=0, avoid=[]):
        """
        Return true if idx is a ring atom
        """

        # Consistency/Termination checks
        if ring_size < 3:
            raise AtomTypesException("ERROR in ring_atom: ring_size variable must be set to an integer greater than 2!")
        if counter == ring_size:
            return False

        # Automatically assign start to the supplied idx value. For recursive calls this is set manually
        if start is None:
            start = idx

        # Loop over connections and recursively search for idx
        cons = [count_i for count_i, i in enumerate(self.adj_mat[idx]) if i == 1 and count_i not in avoid]
        if len(cons) == 0:
            return False
        elif start in cons:
            return True
        else:
            for i in cons:
                if self.ring_atom_new(i, start=start, ring_size=ring_size, counter=counter + 1,
                                      avoid=[idx]): return True
            return False

    def is_nitro(self, i, alt_adj_mat=None):
        """
        Return bool depending on if the atom is a nitro nitrogen atom
        """
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

    def is_sulfoxide(self, i, alt_adj_mat=None):
        """
        Return bool depending on if the atom is a sulfoxide sulfur atom
        """
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

    def is_sulfonyl(self, i, alt_adj_mat=None):
        """
        Return bool depending on if the atom is a sulfonyl sulfur atom
        """
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

    def is_phosphate(self, i, alt_adj_mat=None):
        """
        Return bool depending on if the atom is a phosphate phosphorus atom
        """
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

    def is_cyano(self, i, alt_adj_mat=None):
        """
         # Return bool depending on if the atom is a cyano nitrogen atom
        """
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

    def is_frag_sulfoxide(self, i, alt_adj_mat=None):
        """
    # is function for fragments
    # Return bool depending on if the atom is a sulfoxide sulfur atom
        """
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

    def is_frag_sulfonyl(self, i, alt_adj_mat=None):
        """
        # Return bool depending on if the atom is a sulfonyl sulfur atom
        """
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
    #              at the first match False is returned.
        """
        for i in a_list:
            if np.array_equal(a, i):
                return False
        return True


class AtomType(AtomTypesBaseOperation):
    """
    Atom type related operations
    Unlike base, methods here can be used independently
    TODO: have to check adj_mat changes is always consitent
    """

    def __init__(self, gens: int, structure_data={}, AdjMat=None, adj_mat=None):
        super().__init__(structure_data=structure_data, AdjMat=AdjMat, adj_mat=adj_mat)
        self.gens = gens

        try:
            self.atom_types = [0] * len(self.elements)  # place holder
        except StructureException:
            self.atom_types = [0]

    def build_from_typeadj_fun(self, atomtype):
        """
        build elements and adj_mat from type_adj function
        this would only set elements, adj_mat, so be very careful that if we try to regenerate adj_mat from AdjMat
        it won't work as it doesn't have geometry info
        """
        adj_mat, self.elements = type_adjmat(atomtype)
        self._AdjMat.adj_mat = adj_mat

    def UpdateAtomTypes(self, **kwargs):
        self.atom_types = self.id_types(gens=self.gens, **kwargs)

    # TODO I want to make id_types from being directly called, the custom gens is very dangerous
    def id_types(self, gens=2, avoid=[], fc=None, keep_lone=None, return_index=False, algorithm="matrix"):
        """
        identifies the taffi atom types for all elements from an adjacency matrix/list (A) and element identify.
        keep gens as option for potential individual method call
        """
        mass_dict = Tables.MASSDICT
        e_dict = Tables.ELEMENTDICT

        # If atomic numbers are supplied in place of elements
        try:
            self.elements = [e_dict[int(_)] for _ in self.elements]
        except:
            pass

        # Initialize fc if undefined by user
        if fc is None and keep_lone is None:
            fc = [[]]
            fc[0] = [0] * len(self.elements)
            keep_lone = [[]]

        elif keep_lone is None:
            keep_lone = [[] for i in range(len(fc))]

        elif fc is None:
            fc = [[0] * len(self.elements) for i in range(len(keep_lone))]

        if len(fc[0]) != len(self.elements):
            raise AtomTypesException("ERROR in id_types: fc must have the same dimensions as elements and A")

        # bonding index: refer to which bonding/fc it uses to determine atomtypes
        bond_index = [range(len(fc))] * len(self.elements)

        # check duplication:
        set_fc = list(map(list, set(map(tuple, fc))))
        set_keep_lone = list(map(list, set(map(tuple, keep_lone))))

        total_fc = [fc[i] + keep_lone[i] for i in range(len(fc))]
        set_total = list(map(list, set(map(tuple, total_fc))))
        keep_ind = sorted([next(count_m for count_m, m in enumerate(total_fc) if m == j) for j in set_total])

        fc_0 = deepcopy(fc)
        keeplone_0 = deepcopy(keep_lone)

        if max(len(set_fc), len(set_keep_lone)) == 1:
            fc_0 = fc_0[0]
            keep_lone = keeplone_0[0]

            # Calculate formal charge terms and radical terms
            fc_s = ['' for i in range(len(self.elements))]
            for i in range(len(self.elements)):
                if i in keep_lone: fc_s[i] += '*'
                if fc_0[i] > 0: fc_s[i] += abs(fc_0[i]) * "+"
                if fc_0[i] < 0: fc_s[i] += abs(fc_0[i]) * "-"
            fc_0 = fc_s

            # Assemble prerequisite masses and Loop over the inidices that need to be id'ed
            masses = [mass_dict[i] for i in self.elements]
            N_masses = deepcopy(masses)
            for i in range(len(self.elements)):
                N_masses[i] += (fc_0[i].count('+') * 100.0 + fc_0[i].count('-') * 90.0 + fc_0[i].count('*') * 80.0)

            if algorithm == "matrix":
                atom_types = ["[" + self.taffi_type(i, N_masses, gens, fc=fc_0) + "]" for i in
                              range(len(self.elements))]
            elif algorithm == "list":
                # not really used anymore
                atom_types = ["[" + self.taffi_type_list(i, self.adj_mat, N_masses, gens) + "]" for i in
                              range(len(self.elements))]

        # resonance structure appear, identify special atoms and keep both formal charge information (now only support matrix input)
        else:
            # Assemble prerequisite masses
            masses = [mass_dict[i] for i in self.elements]
            charge_atoms = [[index for index, charge in enumerate(fc_i) if charge != 0] for fc_i in
                            fc_0]  # find charge contained atoms
            CR_atoms = [charge_atoms[i] + keeplone_0[i] for i in range(len(fc_0))]
            keep_CR_atoms = [charge_atoms[i] + keeplone_0[i] for i in keep_ind]  # equal to set of CR_atom
            special_atoms = [index for index in list(set(chain.from_iterable(keep_CR_atoms))) if
                             list(chain.from_iterable(keep_CR_atoms)).count(index) < len(set_fc) * len(
                                 set_keep_lone)]  # find resonance atoms
            normal_atoms = [ind for ind in range(len(self.elements)) if ind not in special_atoms]
            atom_types = []

            # Calculate formal charge terms
            for _ in range(len(fc_0)):
                fc_s = ['' for i in range(len(self.elements))]
                for i in range(len(self.elements)):
                    if i in keeplone_0[_]: fc_s[i] += '*'
                    if fc_0[_][i] > 0: fc_s[i] += abs(fc_0[_][i]) * "+"
                    if fc_0[_][i] < 0: fc_s[i] += abs(fc_0[_][i]) * "-"
                fc_0[_] = fc_s

            for ind in range(len(self.elements)):
                if ind in normal_atoms:
                    # Graphical separations are used for determining which atoms and bonds to keep
                    gs = self.AdjMat.graph_seps()

                    # all atoms within "gens" of the ind atoms are kept
                    keep_atoms = list(set([count_j for count_j, j in enumerate(gs[ind]) if j <= gens]))
                    contain_special = [N_s for N_s in keep_atoms if N_s in special_atoms]
                    N_special = len(contain_special)

                    # if N_special = 0 select first formal charge
                    if N_special == 0:
                        fc = fc_0[0]

                    # if N_special = 1, select formal charge for that special atom
                    elif N_special == 1:
                        bond_ind = [l for l in range(len(fc_0)) if contain_special[0] in CR_atoms[l]]
                        fc = fc_0[bond_ind[0]]
                        bond_index[ind] = sorted(bond_ind)

                    # if N_special >= 2, introduce additional Criteria to determine the bond matrix
                    else:
                        fc_criteria = [0] * len(fc_0)
                        # find the nearest special atom
                        nearest_special_atom = [N_s for N_s in special_atoms if self.adj_mat[ind][N_s] == 1]
                        for l in range(len(fc_0)):
                            fc_criteria[l] = -len(
                                [index for index in nearest_special_atom if index in CR_atoms[l]]) - 0.1 * len(
                                [index for index in contain_special if
                                 index not in nearest_special_atom and index in CR_atoms[l]])

                        bond_ind = [bind for bind, cr in enumerate(fc_criteria) if cr == min(fc_criteria)]
                        fc = fc_0[bond_ind[0]]
                        bond_index[ind] = sorted(bond_ind)

                else:
                    bond_ind = [l for l in range(len(fc_0)) if ind in CR_atoms[l]]
                    fc = fc_0[bond_ind[0]]
                    bond_index[ind] = bond_ind

                # add charge to atom_type sorting
                N_masses = deepcopy(masses)
                for i in range(len(self.elements)):
                    N_masses[i] += (fc[i].count('+') * 100.0 + fc[i].count('-') * 90.0 + fc[i].count('*') * 80.0)

                atom_types += ["[" + self.taffi_type(ind, N_masses, gens, fc=fc) + "]"]

        # Add ring atom designation for atom types that belong are intrinsic to rings
        # (depdends on the value of gens)
        for count_i, i in enumerate(atom_types):
            if self.ring_atom_new(count_i, ring_size=(gens + 2)) == True:
                atom_types[count_i] = "R" + atom_types[count_i]

        if return_index:
            return atom_types, bond_index
        else:
            return atom_types

    def find_lewis(self, bonding_pref=[], q_tot=0, fixed_bonds=[], fc_0=None, keep_lone=[], return_pref=False,
                   verbose=False, b_mat_only=False, return_FC=False, octet_opt=True, check_lewis_flag=False):
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
        # keep the ability to place with these tables manually
        find_lewis_dict = SimpleNamespace()
        find_lewis_dict.lone_e = Tables.LONE_ELECTRON
        find_lewis_dict.periodic = Tables.PERIODICT
        find_lewis_dict.en = Tables.ELETRONEGATIVITY
        find_lewis_dict.pol = Tables.POLARIZABILITY
        find_lewis_dict.be = Tables.BOND_ENERGY
        find_lewis_dict.atomic_to_element = {find_lewis_dict.periodic[i]: i for i in find_lewis_dict.periodic.keys()}

        # Consistency check on fc_0 argument, if supplied
        if fc_0 is not None:
            if len(fc_0) != len(self.elements):
                raise AtomTypesException(
                    "ERROR in find_lewis: the fc_0 and elements lists must have the same dimensions.")
            if int(sum(fc_0)) != int(q_tot):
                raise AtomTypesException("ERROR in find_lewis: the sum of formal charges does not equal q_tot.")

        # Initalize elementa and atomic_number lists for use by the function
        atomic_number = [find_lewis_dict.periodic[i.lower()] for i in self.elements]
        adj_mat = deepcopy(self.adj_mat)

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
                raise AtomTypesException("ERROR in find_lewis: the algorithm isn't compatible with atomic numbers "
                                         "greater than 54 owing to a lack of rules for treating lanthanides. Exiting...")
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
                bonding_target[count_i] = N_tot - find_lewis_dict.lone_e[self.elements[count_i].lower()]

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
                O_ind = [count_j for count_j, j in enumerate(self.adj_mat[i]) if j == 1 and
                         self.elements[count_j].lower() == "o" and sum(self.adj_mat[count_j]) == 1]
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
            if self.is_cyano(i) is True:
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
            if self.is_isocyano(i, alt_adj_mat=adj_mat) is True:
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
                raise AtomTypesException("ERROR in find_lewis: fixed_bonds requests bond creation between \
                                         atoms {} and {} ({} bonds) but the adjacency matrix doesn't \
                                         reflect a bond. Exiting...".format(a, b, N))

            # Check that less than or an equal number of bonds exist between these atoms than is requested
            if N_current > N:
                raise AtomTypesException("ERROR in find_lewis: fixed_bonds requests bond creation between \
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
                print(
                    "{:40s} {:<20d} {:<20d} {:<20d} {:<20d} {}".format(self.elements[count_i], lone_electrons[count_i],
                                                                       bonding_electrons[count_i],
                                                                       core_electrons[count_i],
                                                                       valence_list[count_i] - bonding_electrons[
                                                                           count_i] - lone_electrons[count_i],
                                                                       ",".join(["{}".format(count_j) for count_j, j in
                                                                                 enumerate(adj_mat[count_i]) if
                                                                                 j == 1])))

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
                            lewis_bonded_radicals = [(-find_lewis_dict.en[self.elements[count_j].lower()], count_j) for
                                                     count_j, j in enumerate(self.adj_mat[i]) if
                                                     j == 1 and lewis_lone_electrons[lewis_counter][count_j] % 2 != 0 \
                                                     and 2 * (lewis_bonding_electrons[lewis_counter][count_j] + 1) + (
                                                                 lewis_lone_electrons[lewis_counter][count_j] - 1) <=
                                                     lewis_valence[lewis_counter][count_j] \
                                                     and lewis_lone_electrons[lewis_counter][
                                                         count_j] - 1 >= 0 and count_j not in happy]

                            lewis_bonded_lonepairs = [(-find_lewis_dict.en[self.elements[count_j].lower()], count_j) for
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
            unsatisfy = [count_t for count_t, te in enumerate(total_electron) if te > 2 and te < 8 and te % 2 == 0]
            for uns in unsatisfy:
                full_connect = [count_i for count_i, i in enumerate(self.adj_mat[uns]) if i == 1
                                and total_electron[count_i] == 8 and lewis_lone_electrons[lewis_counter][count_i] >= 2]
                if len(full_connect) > 0:
                    lewis_lone_electrons[lewis_counter][full_connect[0]] -= 2
                    lewis_bonding_electrons[lewis_counter][uns] += 1
                    lewis_bonding_electrons[lewis_counter][full_connect[0]] += 1
                    lewis_adj_mat[lewis_counter][uns][full_connect[0]] += 1
                    lewis_adj_mat[lewis_counter][full_connect[0]][uns] += 1

            # Delete last entry in the lewis np.arrays if the electronic structure is not unique: introduce identical_mat includes both info of bond_mats and formal_charges
            identical_mat = np.vstack([lewis_adj_mat[-1], np.array([valence_list[k] - lewis_bonding_electrons[-1][k]
                                                                    - lewis_lone_electrons[-1][k] for k in
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
                        if (lewis_bonding_electrons[i][count_j] * 2 + lewis_lone_electrons[i][count_j] != 8
                                and lewis_bonding_electrons[i][count_j] * 2 + lewis_lone_electrons[i][count_j] != 18):
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
                sum([find_lewis_dict.be[cbm] if cbm in find_lewis_dict.be.keys() else -10000.0 for cbm in
                     count_bonds_made])]
        # normalize the effect
        lewis_bonds_energy = [-be / max(1, max(lewis_bonds_energy)) for be in lewis_bonds_energy]

        ## Find the total formal charge for each structure
        formal_charges_sums = []
        for i in range(len(lewis_lone_electrons)):
            fc = 0
            for j in range(len(self.elements)):
                fc += valence_list[j] - lewis_bonding_electrons[i][j] - lewis_lone_electrons[i][j]
            formal_charges_sums.append(fc)

        ## Find formal charge eletronegativity contribution
        lewis_formal_charge = [[valence_list[i] - lewis_bonding_electrons[_][i] - lewis_lone_electrons[_][i] for i in
                                range(len(self.elements))] for _ in range(len(lewis_lone_electrons))]
        lewis_keep_lone = [[count_i for count_i, i in enumerate(lone) if i % 2 != 0] for lone in lewis_lone_electrons]
        lewis_fc_en = []  # Electronegativity for stabling charge/radical
        lewis_fc_pol = []  # Polarizability for stabling charge/radical
        lewis_fc_hc = []  # Hyper-conjugation contribution
        for i in range(len(lewis_lone_electrons)):
            formal_charge = lewis_formal_charge[i]
            radical_atom = lewis_keep_lone[i]
            fc_ind = [(count_j, j) for count_j, j in enumerate(formal_charge) if j != 0]
            for R_ind in radical_atom:  # assign +0.5 for radical
                fc_ind += [(R_ind, 0.5)]

            # initialize en,pol and hc
            fc_en, fc_pol, fc_hc = 0, 0, 0

            # Loop over formal charges and radicals
            for count_fc in fc_ind:
                ind = count_fc[0]
                charge = count_fc[1]
                # Count the self contribution: (-) on the most electronegative atom and (+) on the least electronegative atom
                fc_en += 10 * charge * find_lewis_dict.en[self.elements[ind].lower()]

                # Find the nearest and next-nearest atoms for each formal_charge/radical contained atom
                gs = self.AdjMat.graph_seps()
                nearest_atoms = [count_k for count_k, k in enumerate(lewis_adj_mat[i][ind]) if k >= 1]
                NN_atoms = list(set([count_j for count_j, j in enumerate(gs[ind]) if j == 2]))

                # only count when en > en(C)
                fc_en += charge * (
                            sum([find_lewis_dict.en[self.elements[count_k].lower()] for count_k in nearest_atoms if
                                 find_lewis_dict.en[self.elements[count_k].lower()] > 2.54]) +
                            sum([find_lewis_dict.en[self.elements[count_k].lower()] for count_k in NN_atoms if
                                 find_lewis_dict.en[self.elements[count_k].lower()] > 2.54]) * 0.1)

                if charge < 0:  # Polarizability only affects negative charge
                    fc_pol += charge * sum(
                        [find_lewis_dict.pol[self.elements[count_k].lower()] for count_k in nearest_atoms])

                # find hyper-conjugation strcuture
                nearby_carbon = [nind for nind in nearest_atoms if self.elements[nind].lower() == 'c']
                for carbon_ind in nearby_carbon:
                    carbon_nearby = [nind for nind in NN_atoms if
                                     lewis_adj_mat[i][carbon_ind][nind] >= 1 and self.elements[nind].lower() in ['c',
                                                                                                                 'h']]
                    if len(carbon_nearby) == 3: fc_hc -= charge * (len([nind for nind in carbon_nearby
                                                                        if self.elements[nind].lower() == 'c']) * 2
                                                                   + len(
                                [nind for nind in carbon_nearby if self.elements[nind].lower() == 'h']))

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
            # lewis_criteria.append( 10.0*octet_violations[i] + abs(formal_charges_sums[i]) +
            # 0.1*sum([ 1 for j in lewis_lone_electrons[i] if j % 2 != 0 ]) + 0.001*lewis_bonds_energy[i]
            # + 0.00001*lewis_fc_en[i] + 0.000001*lewis_fc_pol[i] + 0.0000001*lewis_fc_hc[i])
            lewis_criteria.append(10.0 * octet_violations[i] + abs(formal_charges_sums[i]) +
                                  0.1 * sum([1 for j in lewis_lone_electrons[i] if j % 2 != 0]) + 0.01 * lewis_fc_en[
                                      i] + 0.005 * lewis_fc_pol[i] + \
                                  0.0001 * lewis_fc_hc[i] + 0.0001 * lewis_bonds_energy[i])
        best_lewis = [i[0] for i in sorted(enumerate(lewis_criteria), key=lambda x: x[
            1])]  # sort from least to most and return a list containing the origial list's indices in the correct order
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
                    print(
                        "{:<40s} {}    {} {}".format(self.elements[j], " ".join([str(k) for k in lewis_adj_mat[i][j]]),
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

    def frag_find_lewis(self, bonding_pref=[], fixed_bonds=[], q_tot=0, fc_0=None, keep_lone=[],
                        return_pref=False, return_FC=False, octet_opt=True, check_lewis_flag=False):
        """
        This function is used to deal with fragments. frag_find_lewis will just form double/triple bonds based on given fc/keep_lone information rather than generate fc/keep_lone.

        :param bonding_pref: optional list of (index, bond_number) tuples that sets the target bond number of the indexed atoms
        :param q_tot:     total charge on the molecule
        :param fixed_bonds: optional list of (index_1,index_2,bond_number) tuples that creates fixed bonds between the index_1
                               and index_2 atoms. No further bonds will be added or subtracted between these atoms.
        :param fc_0:      a list of formal charges on each atom
        :param keep_lone: a list of atom index for which contains a radical

        :returns: lone_electrons:
                  bonding_electrons:
                  core_electrons:
                  bond_mat:  an NxN matrix holding the bond orders between all atoms in the adj_mat
                  bonding_pref (optinal): optional list of (index, bond_number) tuples that sets the target bond number of the indexed atoms

        """
        frag_find_lewis_dict = SimpleNamespace()
        frag_find_lewis_dict.lone_e = Tables.LONE_ELECTRON
        frag_find_lewis_dict.periodic = Tables.PERIODICT
        frag_find_lewis_dict.en = Tables.ELETRONEGATIVITY
        frag_find_lewis_dict.pol = Tables.POLARIZABILITY
        frag_find_lewis_dict.be = Tables.BOND_ENERGY
        frag_find_lewis_dict.atomic_to_element = {frag_find_lewis_dict.periodic[i]: i for i in
                                                  frag_find_lewis_dict.periodic.keys()}

        # Consistency check on fc_0 argument, if supplied
        if fc_0 is not None:
            if len(fc_0) != len(self.elements):
                raise AtomTypesException(
                    "ERROR in frag_find_lewis: the fc_0 and elements lists must have the same dimensions.")
            if int(sum(fc_0)) != int(q_tot):
                raise AtomTypesException("ERROR in frag_find_lewis: the sum of formal charges does not equal q_tot.")

        # Initalize elementa and atomic_number lists for use by the function
        atomic_number = [frag_find_lewis_dict.periodic[i.lower()] for i in self.elements]
        adj_mat = deepcopy(self.adj_mat)

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
                print(
                    "ERROR in frag_find_lewis: the algorithm isn't compatible with atomic numbers greater than 54 owing to a lack of rules for treating lanthanides. Exiting...")
                quit()
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
                bonding_target[count_i] = N_tot - frag_find_lewis_dict.lone_e[self.elements[count_i].lower()]

                # Loop over the adjmat and assign initial bonded electrons assuming single bonds (and adjust lone electrons accordingly)
        for count_i, i in enumerate(self.adj_mat):
            bonding_electrons[count_i] += sum(i)
            lone_electrons[count_i] -= sum(i)

        # Apply keep_lone: add one electron to such index
        for count_i in keep_lone:
            lone_electrons[count_i] += 1
            bonding_target[count_i] -= 1

        # Eliminate all radicals by forming higher order bonds
        change_list = range(len(lone_electrons))
        bonds_made = []
        loop_list = [(atomic_number[i], i) for i in range(len(lone_electrons))]
        loop_list = [i[1] for i in sorted(loop_list)]

        # Loop over bonding_pref, find whether exist two of them can form bonds
        if len(bonding_pref) > 1:
            bonding_pref_ind = [i[0] for i in bonding_pref]
            comb = combinations(bonding_pref_ind, 2)
            pref_pair = [sorted(pair) for pair in comb if self.adj_mat[pair[0]][pair[1]] == 1]
        else:
            pref_pair = []

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
            if self.is_frag_sulfoxide(i) is True:
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
            if self.is_frag_sulfonyl(i) is True:
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
            if self.is_cyano(i) is True:
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
            if self.is_isocyano(i, alt_adj_mat=adj_mat) is True:
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
                raise AtomTypesException(
                    "ERROR in frag_find_lewis: fixed_bonds requests bond creation between atoms {} and {} ({} bonds) but the adjacency matrix doesn't reflect a bond. Exiting...".format(
                        a, b, N))

            # Check that less than or an equal number of bonds exist between these atoms than is requested
            if N_current > N:
                raise AtomTypesException(
                    "ERROR in frag_find_lewis: fixed_bonds requests bond creation between atoms {} and {} ({} bonds) but {} bonds already exist between these atoms. There may be a conflict between the special groups handling and the requested lewis_structure.".format(
                        a, b, N, N_current))

            # Check that enough lone electrons exists on each atom to reach the target bond number
            if lone_electrons[a] < (N - N_current):
                print(
                    "Warning in frag_find_lewis: fixed_bonds requests bond creation between atoms {} and {} ({} bonds)  but atom {} only has {} lone electrons.".format(
                        a, b, N, self.elements[a], lone_electrons[a]))

            # Check that enough lone electrons exists on each atom to reach the target bond number
            if lone_electrons[b] < (N - N_current):
                print(
                    "Warning in frag_find_lewis: fixed_bonds requests bond creation between atoms {} and {} ({} bonds) but atom {} only has {} lone electrons.".format(
                        a, b, N, self.elements[b], lone_electrons[b]))

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
                    #    print "ERROR in find_lewis: atom ({}, index {}) doesn't have enough lone electrons ({}) to be removed to satisfy the specified formal charge ({}).".format(self.elements[count_i],count_i,lone_electrons[count_i],i)
                    #    quit()
                    lone_electrons[count_i] = lone_electrons[count_i] - i
                if i < 0:
                    lone_electrons[count_i] = lone_electrons[count_i] + int(abs(i))
            q_tot = 0

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
                            lewis_bonded_lonepairs = [
                                (-frag_find_lewis_dict.en[self.elements[count_j].lower()], count_j) for
                                count_j, j in enumerate(self.adj_mat[i]) if
                                j == 1 and lewis_lone_electrons[lewis_counter][count_j] > 0 \
                                and 2 * (lewis_bonding_electrons[lewis_counter][count_j] + 1) + (
                                        lewis_lone_electrons[lewis_counter][count_j] - 1) <=
                                lewis_valence[lewis_counter][count_j] and \
                                lewis_lone_electrons[lewis_counter][
                                    count_j] - 1 >= 0 and count_j not in happy]

                            # Sort by atomic number (cheap way of sorting carbon before other atoms, should probably switch over to electronegativities)
                            lewis_bonded_lonepairs = [j[1] for j in sorted(lewis_bonded_lonepairs)]

                            # Try to form a bond with a neighboring atom with spare lone electrons (valence check ensures that no improper 5-bonded atoms are formed)
                            if len(lewis_bonded_lonepairs) > 0:
                                lewis_bonding_electrons[lewis_counter][i] += 1
                                lewis_bonding_electrons[lewis_counter][lewis_bonded_lonepairs[0]] += 1
                                lewis_adj_mat[lewis_counter][i][lewis_bonded_lonepairs[0]] += 1
                                lewis_adj_mat[lewis_counter][lewis_bonded_lonepairs[0]][i] += 1
                                lewis_lone_electrons[lewis_counter][i] -= 1
                                lewis_lone_electrons[lewis_counter][lewis_bonded_lonepairs[0]] -= 1
                                change_list += [i, lewis_bonded_lonepairs[0]]
                                lewis_bonds_made[lewis_counter] += [(i, lewis_bonded_lonepairs[0])]

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

                        potential_bond = [i for i in lewis_bonds_made[lewis_counter] if (
                                (i[0] in ind or i[1] in ind) and (
                                i[0] not in bonding_pref_ind or i[1] not in bonding_pref_ind))]
                        if len(potential_bond) == 0:
                            potential_bond = [i for i in lewis_bonds_made[lewis_counter] if i[0] in ind or i[1] in ind]

                            # Check if a rearrangment is possible, break if none are available
                        try:
                            break_bond = next(i for i in potential_bond)
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
                        print(
                            "WARNING: maximum attempts to establish a lewis-structure consistent with the user supplied bonding preference has been exceeded ({}).".format(
                                outer_max_cycles))
                        break

            # Re-apply keep_lone: remove one electron from such index
            for count_i in keep_lone:
                lewis_lone_electrons[lewis_counter][count_i] -= 1

            # Delete last entry in the lewis np.arrays if the electronic structure is not unique
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
                    if j.lower() in ["c", "n", "o", "f", "si", "p", "s", "cl"] and count_j not in bonding_pref_ind:
                        if lewis_bonding_electrons[i][count_j] * 2 + lewis_lone_electrons[i][count_j] < 8:
                            ov += (8 - lewis_bonding_electrons[i][count_j] * 2 - lewis_lone_electrons[i][count_j])
                        else:
                            ov += 2 * (lewis_bonding_electrons[i][count_j] * 2 + lewis_lone_electrons[i][count_j] - 8)
                    if j.lower() in ["br", 'i'] and count_j not in bonding_pref_ind:
                        if lewis_bonding_electrons[i][count_j] * 2 + lewis_lone_electrons[i][count_j] < 18:
                            ov += (18 - lewis_bonding_electrons[i][count_j] * 2 - lewis_lone_electrons[i][count_j])
                        else:
                            ov += 2 * (lewis_bonding_electrons[i][count_j] * 2 + lewis_lone_electrons[i][count_j] - 18)
            octet_violations.append(ov)

        # Find the total formal charge for each structure
        formal_charges_sums = []
        for i in range(len(lewis_lone_electrons)):
            fc = 0
            for j in range(len(self.elements)):
                fc += valence_list[j] - lewis_bonding_electrons[i][j] - lewis_lone_electrons[i][j]
            formal_charges_sums.append(fc)

        # Find formal charge eletronegativity contribution
        lewis_fc_en = []  # Electronegativity for stabling charge/radical
        lewis_fc_pol = []  # Polarizability for stabling charge/radical
        lewis_fc_hc = []  # Hyper-conjugation contribution
        formal_charge = deepcopy(fc_0)
        radical_atom = deepcopy(keep_lone)

        # Consider the facts that will stabilize ions/radicals
        for i in range(len(lewis_lone_electrons)):
            fc_ind = [(count_j, j) for count_j, j in enumerate(formal_charge) if j != 0]
            for R_ind in radical_atom:  # assign +0.5 for radical
                fc_ind += [(R_ind, 0.5)]

            # initialize en,pol and hc
            fc_en, fc_pol, fc_hc = 0, 0, 0

            # Loop over formal charges and radicals
            for count_fc in fc_ind:
                ind = count_fc[0]
                charge = count_fc[1]

                # Count the self contribution: (-) on the most electronegative atom and (+) on the least electronegative atom
                fc_en += 10 * charge * frag_find_lewis_dict.en[self.elements[ind].lower()]

                # Find the nearest and next-nearest atoms for each formal_charge/radical contained atom
                gs = self.graph_seps()
                nearest_atoms = [count_k for count_k, k in enumerate(lewis_adj_mat[i][ind]) if k >= 1]
                NN_atoms = list(set([count_j for count_j, j in enumerate(gs[ind]) if j == 2]))

                # only count en > en(C)
                fc_en += charge * (
                            sum([frag_find_lewis_dict.en[self.elements[count_k].lower()] for count_k in nearest_atoms if
                                 frag_find_lewis_dict.en[self.elements[count_k].lower()] > 2.54]) +
                            sum([frag_find_lewis_dict.en[self.elements[count_k].lower()] for count_k in NN_atoms if
                                 frag_find_lewis_dict.en[self.elements[count_k].lower()] > 2.54]) * 0.1)

                if charge < 0:  # Polarizability only affects negative charge ?
                    fc_pol += charge * sum(
                        [frag_find_lewis_dict.pol[self.elements[count_k].lower()] for count_k in nearest_atoms])

                # find hyper-conjugation strcuture
                nearby_carbon = [nind for nind in nearest_atoms if self.elements[nind].lower() == 'c']
                for carbon_ind in nearby_carbon:
                    carbon_nearby = [nind for nind in NN_atoms if
                                     lewis_adj_mat[i][carbon_ind][nind] >= 1 and self.elements[nind].lower() in ['c',
                                                                                                                 'h']]
                    radical_on_carbon = lewis_lone_electrons[i][carbon_ind]
                    if (len(carbon_nearby) + radical_on_carbon) == 3: fc_hc -= charge * (
                            len([nind for nind in carbon_nearby if self.elements[nind].lower() == 'c']) * 2 +
                            len([nind for nind in carbon_nearby if
                                 self.elements[nind].lower() == 'h']) + radical_on_carbon)

            lewis_fc_en.append(fc_en)
            lewis_fc_pol.append(fc_pol)
            lewis_fc_hc.append(fc_hc)
            # normalize the effect
        lewis_fc_en = [lfc / max(1, max(abs(np.array(lewis_fc_en)))) for lfc in lewis_fc_en]
        lewis_fc_pol = [lfp / max(1, max(abs(np.array(lewis_fc_pol)))) for lfp in lewis_fc_pol]

        # The bond formation will try to keep resonance structures
        # First find identical atoms; then determine whether the distance between them is 2, if so, find the connecting atom to this pair to form two "resonance preferance pairs". After
        # finding all such pairs, count it in lewis critria.
        mass_dict = Tables.MASSDICT
        masses = [mass_dict[i] for i in self.elements]
        hash_list = [self.atom_hash(ind, masses) for ind in range(len(self.elements))]

        # Find identical atoms
        same_atoms = [(i, [count_j for count_j, j in enumerate(hash_list) if j == i]) for i in set(hash_list) if
                      hash_list.count(i) > 1]
        gs = self.AdjMat.graph_seps()
        res_atoms = []

        # Keep identical atoms whose distance = 2
        for i in same_atoms:
            if len(i[1]) == 2 and gs[i[1][0]][i[1][1]] == 2: res_atoms += [[tuple(i[1]), 2]]
            if len(i[1]) > 2:
                comb = combinations(i[1], 2)
                res_atoms += [[j, len(i[1])] for j in comb if gs[j[0]][j[1]] == 2]

        # Find connecting atom to form resonance pairs
        res_pair = []
        res_pair_score = []
        for pair_list in res_atoms:
            pair = pair_list[0]
            center_atom = \
            [ind for ind, i in enumerate(self.adj_mat[pair[0]]) if i == 1 and self.adj_mat[pair[1]][ind] == 1][
                0]
            res_pair += [tuple(sorted((pair[0], center_atom))), tuple(sorted((pair[1], center_atom)))]
            res_pair_score += [pair_list[1], pair_list[1]]

        # loop over lewis_bonds_made to determine the lewis criteria
        lewis_res_bonds = []
        lewis_bonds_energy = []
        for bonds_made in lewis_bonds_made:
            res_bonds = 0
            bonds_energy = 0
            for lb, bond_made in enumerate(bonds_made): bonds_made[lb] = tuple(sorted(bond_made))
            for bond_made in bonds_made:
                # in the bonds in res_pair, enlarge by the number of symmetry atoms
                if bond_made in res_pair:
                    res_bonds += res_pair_score[res_pair.index(bond_made)]
                # if the bonds in pref_bonds, enlarge the effect by 2; for cation, prefer to form bonds while for anion/radical prefer to keep single bond
                factor = 0.1  # make bonds_energy comparable with originally used bonds_en
                if sorted([bond_made[0], bond_made[1]]) in pref_pair: factor *= 2
                if fc_0[bond_made[0]] > 0 or fc_0[bond_made[1]] > 0: factor *= 2
                if fc_0[bond_made[0]] < 0 or fc_0[bond_made[1]] < 0: factor *= 0.5
                if bond_made[0] in keep_lone or bond_made[1] in keep_lone: factor *= 0.5
                bond_type = "{}-{}-{}".format(min(atomic_number[bond_made[0]], atomic_number[bond_made[1]]),
                                              max(atomic_number[bond_made[0]], atomic_number[bond_made[1]]),
                                              bonds_made.count(bond_made))
                if bond_type in frag_find_lewis_dict.be.keys():
                    bonds_energy += factor * frag_find_lewis_dict.be[bond_type]
                else:
                    bonds_energy -= factor * (-10000.0)
            lewis_bonds_energy += [bonds_energy]
            lewis_res_bonds += [res_bonds]
        # normalize the effect
        lewis_bonds_energy = [-be / max(1, max(lewis_bonds_energy)) for be in lewis_bonds_energy]
        lewis_res_bonds = [-re / max(1, max(lewis_res_bonds)) for re in lewis_res_bonds]

        # Add the total number of radicals to the total formal charge to determine the criteria.
        # The radical count is scaled by 0.01 and the lone pair count is scaled by 0.001. This results
        # in the structure with the lowest formal charge always being returned, and the radical count
        # only being considered if structures with equivalent formal charges are found, and likewise with
        # the lone pair count. The structure(s) with the lowest score will be returned.
        lewis_criteria = []
        for i in range(len(lewis_lone_electrons)):
            lewis_criteria.append(
                10.0 * octet_violations[i] + lewis_fc_en[i] + lewis_fc_pol[i] + 0.1 * lewis_fc_hc[i] + 0.01 *
                lewis_bonds_energy[i] + 0.0001 * lewis_res_bonds[i])

        best_lewis = [i[0] for i in sorted(enumerate(lewis_criteria), key=lambda x: x[
            1])]  # sort from least to most and return a list containing the origial list's indices in the correct order
        best_lewis = [i for i in best_lewis if lewis_criteria[i] == lewis_criteria[best_lewis[0]]]

        # Apply keep_lone information, remove the electron to form lone electron
        for i in best_lewis:
            for j in keep_lone:
                lewis_lone_electrons[i][j] -= 1

        # return check_lewis function
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
                    [lewis_core_electrons[_] for _ in best_lewis], [lewis_adj_mat[_] for _ in best_lewis], \
                    [[valence_list[i] - lewis_bonding_electrons[_][i] - lewis_lone_electrons[_][i] for i in
                      range(len(self.elements))] for _ in best_lewis], bonding_pref

        else:
            if return_FC is False:
                return [lewis_lone_electrons[_] for _ in best_lewis], [lewis_bonding_electrons[_] for _ in best_lewis], \
                    [lewis_core_electrons[_] for _ in best_lewis], [lewis_adj_mat[_] for _ in best_lewis]
            else:
                return [lewis_lone_electrons[_] for _ in best_lewis], [lewis_bonding_electrons[_] for _ in best_lewis], \
                    [lewis_core_electrons[_] for _ in best_lewis], [lewis_adj_mat[_] for _ in best_lewis], \
                    [[valence_list[i] - lewis_bonding_electrons[_][i] - lewis_lone_electrons[_][i] for i in
                      range(len(self.elements))] for _ in best_lewis]

def main(argv):
    test = AtomType(2)
    test.elements = ['O','S']
    print(test.elements)
    print(test._AdjMat.elements)

if __name__ == '__main__':
    main(sys.argv[1:])