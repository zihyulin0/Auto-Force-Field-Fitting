from model.structure import StructureBase
from model.adjacency import AdjacencyMatrix
from model.atomtypes import AtomType
from utilities.PeriodictTable import Tables
from utilities.Excpetions import TAFFIException
from copy import deepcopy
from utilities.writers import write_xyz

class MoleculeException(TAFFIException):
    pass

class Molecule(StructureBase):
    """
    whenever elements, geometry, q_tot being set, it will change that in Molecule, AdjMat, Atomtypes
    but when adj_mat or AdjMat changes, need to manually update that in Atomtypes
    """
    def __init__(self, gens):
        super().__init__()
        self.gens = gens
        # just the place holder here, since atomtypes would require a build adjmat
        self.AdjMat = AdjacencyMatrix()
        self.AtomType = AtomType(self.gens)
        self.bond_mat = None
        self.hash_list = None
        self.fc = None
        # model compound use only,
        # TODO we might need to think about whether to have an class for frac or model molecule
        self.mc_prop = {}
        # internal use only
        self._lone_electrons = None
        self._core_electrons = None
        self._bonding_electrons = None

    @StructureBase.elements.setter
    def elements(self, value):
        self.AdjMat.elements = value
        self.AtomType.elements = value
        super(Molecule, Molecule).elements.__set__(self, value)
    @StructureBase.geometry.setter
    def geometry(self, value):
        self.AdjMat.geometry = value
        self.AtomType.geometry = value
        super(Molecule, Molecule).geometry.__set__(self, value)

    @StructureBase.q_tot.setter
    def q_tot(self, value):
        self.AdjMat.q_tot = value
        self.AtomType.q_tot = value
        super(Molecule, Molecule).q_tot.__set__(self, value)

    @property
    def adj_mat(self):
        return self.AdjMat.adj_mat

    @property
    def atom_types(self):
        return self.AtomType.atom_types

    def parse_data(self, **kwargs):
        elements = kwargs.get('elements', None)
        geometry = kwargs.get('geometry', None)
        q_tot = kwargs.get('q_tot', None)
        bond_mat = kwargs.get('bond_mat', None)
        atom_types = kwargs.get('atom_types', None)
        adj_mat = kwargs.get('adj_mat', None)
        hash_list = kwargs.get('hash_list', None)

        if elements is not None:
            self.elements = elements
        if geometry is not None:
            self.geometry = geometry
        if q_tot is not None:
            self.q_tot = q_tot
        if bond_mat is not None:
            self.bond_mat = bond_mat
        if hash_list is not None:
            self.hash_list = hash_list
        if atom_types is not None:
            self.AtomType.atom_types = atom_types
        if adj_mat is not None:
            self.AdjMat.adj_mat = adj_mat
            self.AtomType.parse_data_from_AdjMat(self.AdjMat)


    def get_atom_types_pipe(self, xyz, method):
        """
        A sophisticate pipeline to get atom types starting from reading xyz and does canonicalization
        """
        # parse geometry to molecule
        self.parse_from_xyz(xyz)
        # build adj_mat
        self.AdjMat.build_adj_mat()
        self.AtomType.parse_data_from_AdjMat(self.AdjMat)
        if method == 'complicate':
            # involves canoicalize and formal charge
            # canonicalize the geometry
            self.canon_geo()
            self.update_atomtypes_extra()
        elif method == 'simple':
            self.AtomType.UpdateAtomTypes()
        else:
            raise MoleculeException('specified pipeline {} to generate atomtypes does not exist'.format(method))

    def update_bondmat_fc(self):
        """
        get bonding matrix and formal charge
        """
        self._lone_electrons, self._bonding_electrons, self._core_electrons, self.bond_mat, self.fc = \
            self.AdjMat.find_lewis(return_pref=False, verbose=False, return_FC=True)

    def update_atomtypes_extra(self):
        """
        update atom types with extra procedure with the radicals
        """
        if self._core_electrons is None:
            self.update_bondmat_fc()

        # find radicals for atom type assignment
        keep_lones = [[count_j for count_j, j in enumerate(self._lone_electron) if j % 2 != 0] for self._lone_electron in
                      self._lone_electrons]
        # finally determine atom types
        self.AtomType.UpdateAtomTypes(fc=self.fc, keep_lone=keep_lones, return_index=False)

    def write_xyz(self, file_prefix):
        write_xyz(file_prefix, self.elements, self.geometry, additional_col=self.atom_types)

    def canon_geo(self):
        """
        Canonicalizes the ordering of atoms in a geometry based on a hash function. Atoms that hash to equivalent
        values retain their relative order from the input geometry.

        This will update geo, adj_mat, elements, atomtypes (but that's probably still empty) since canon_geo will be
        called before that
        """

        # Canonicalize by sorting the elements based on hashing
        masses = [Tables.MASSDICT[i] for i in self.elements]
        hash_list, atoms = [list(j) for j in
                            zip(*sorted([(self.AtomType.atom_hash(i, masses), i) for i in range(len(self.geometry))], reverse=True))]

        # Update lists/arrays based on atoms
        self.geometry = self.geometry[atoms]
        self.elements = [self.elements[i] for i in atoms]
        # update adj_mat
        self.AdjMat._adj_mat = self.AdjMat._adj_mat[atoms]
        self.AdjMat._adj_mat = self.AdjMat._adj_mat[:, atoms]
        self.AtomType.parse_data_from_AdjMat(self.AdjMat)
        # update atom_types (at this point, it's probably still just an empty list
        self.AtomType.atom_types = [self.AtomType.atom_types[i] for i in atoms]


    def find_modes(self, return_all=0, d1_opt=False):
        # A wrapper for the commands to parse the dihedrals from the adjacency matrix and geometry.
        #           Atom_types isn't necessary here, this section of code just hasn't been cleaned up.
        # Returns:  list of (dihedral_type,angle) tuples.


        # Initialize lists of each instance and type of FF object.
        # instances are stored as tuples of the atoms involved
        # (e.g., bonds between atoms 1 and 13 and 17 and 5 would be stored as [(1,13),(17,5)]
        # Similarly, types are stored as tuples of atom types.
        Atom_types = deepcopy(self.atom_types)
        Atom_types = [next(j for j in i.split('link-') if j != '') for i in
                      Atom_types]  # split('-link') call is necessary for handling fragment atoms
        Bonds = []
        Angles = []
        Dihedrals = []
        Dihedral_types = []
        One_fives = []

        # Gen = 1 types are generated for angle and dihedral terms
        if d1_opt:
            g1_types = []
            for i in Atom_types:
                tmp = AtomType()
                tmp.build_from_typeadj_fun(i)
                g1_types += [tmp.id_types(gens=1)[0]]

        # Find bonds #
        for count_i, i in enumerate(self.adj_mat):
            Bonds += [canon_bond((Atom_types[count_i], Atom_types[count_j]), (count_i, count_j)) for count_j, j in
                      enumerate(i) if j == 1 and count_j > count_i]
        Bond_types, Bonds = map(list, zip(*Bonds))

        # Remove -UA tag from Bond_types (united-atom has no meaning for bonds)
        Bond_types = [(i[0].split('-UA')[0], i[1].split('-UA')[0]) for i in Bond_types]

        # Find angles #
        # When d1_opt is supplied, the 1 and 3 atoms are typed based on depth=1 types
        for i in Bonds:
            if d1_opt:
                Angles += [canon_angle((g1_types[count_j], Atom_types[i[0]], g1_types[i[1]]), (count_j, i[0], i[1])) for
                           count_j, j in enumerate(self.adj_mat[i[0]]) if j == 1 and count_j != i[1]]
                Angles += [canon_angle((g1_types[i[0]], Atom_types[i[1]], g1_types[count_j]), (i[0], i[1], count_j)) for
                           count_j, j in enumerate(self.adj_mat[i[1]]) if j == 1 and count_j != i[0]]
            else:
                Angles += [canon_angle((Atom_types[count_j], Atom_types[i[0]], Atom_types[i[1]]), (count_j, i[0], i[1]))
                           for count_j, j in enumerate(self.adj_mat[i[0]]) if j == 1 and count_j != i[1]]
                Angles += [canon_angle((Atom_types[i[0]], Atom_types[i[1]], Atom_types[count_j]), (i[0], i[1], count_j))
                           for count_j, j in enumerate(self.adj_mat[i[1]]) if j == 1 and count_j != i[0]]

        Angle_types, Angles = remove_duplicate_modes(*map(list, zip(*Angles)))

        # Remove -UA tag from Angle_types (united-atom has no meaning for angles)
        Angle_types = [(i[0].split('-UA')[0], i[1].split('-UA')[0], i[2].split('-UA')[0]) for i in Angle_types]

        # Find dihedrals #
        # When d1_opt is supplied, the 1 and 4 atoms are typed based on depth=1 types
        for i in Angles:
            if d1_opt:
                Dihedrals += [canon_dihedral((g1_types[count_j], Atom_types[i[0]], Atom_types[i[1]], g1_types[i[2]]),
                                             (count_j, i[0], i[1], i[2])) for count_j, j in enumerate(self.adj_mat[i[0]]) if
                              j == 1 and count_j not in [i[1], i[2]]]
                Dihedrals += [canon_dihedral((g1_types[i[0]], Atom_types[i[1]], Atom_types[i[2]], g1_types[count_j]),
                                             (i[0], i[1], i[2], count_j)) for count_j, j in enumerate(self.adj_mat[i[2]]) if
                              j == 1 and count_j not in [i[0], i[1]]]
            else:
                Dihedrals += [
                    canon_dihedral((Atom_types[count_j], Atom_types[i[0]], Atom_types[i[1]], Atom_types[i[2]]),
                                   (count_j, i[0], i[1], i[2])) for count_j, j in enumerate(self.adj_mat[i[0]]) if
                    j == 1 and count_j not in [i[1], i[2]]]
                Dihedrals += [
                    canon_dihedral((Atom_types[i[0]], Atom_types[i[1]], Atom_types[i[2]], Atom_types[count_j]),
                                   (i[0], i[1], i[2], count_j)) for count_j, j in enumerate(self.adj_mat[i[2]]) if
                    j == 1 and count_j not in [i[0], i[1]]]
        if len(Dihedrals) != 0:
            Dihedral_types, Dihedrals = remove_duplicate_modes(*map(list, zip(*Dihedrals)))

        # Add Dihedral_type to dihedrals
        for count_i, i in enumerate(Dihedrals):
            if 2 in [j[i[1], i[2]] for j in self.bond_mat]:
                Dihedral_types[count_i] = tuple(list(Dihedral_types[count_i]) + ["harmonic"])
            else:
                Dihedral_types[count_i] = tuple(list(Dihedral_types[count_i]) + ["opls"])

                # Find 1-5s
        # NOTE: no effort is made to sort based on types because these are only used for coul and lj corrections
        for i in Dihedrals:
            # Find atoms attached to first atom of each dihedral
            One_fives += [(count_j, i[0], i[1], i[2], i[3]) for count_j, j in enumerate(self.adj_mat[i[0]]) if
                          j == 1 and count_j not in [i[1], i[2], i[3]]]

            # Find atoms attached to the fourth atom of each dihedral
            One_fives += [(i[0], i[1], i[2], i[3], count_j) for count_j, j in enumerate(self.adj_mat[i[3]]) if
                          j == 1 and count_j not in [i[0], i[1], i[2]]]

        One_five_types = [(Atom_types[i[0]], Atom_types[i[1]], Atom_types[i[2]], Atom_types[i[3]], Atom_types[i[4]]) for
                          i in One_fives]

        if return_all == 1:
            return Bonds, Angles, Dihedrals, One_fives, Bond_types, Angle_types, Dihedral_types, One_five_types
        else:
            return Bonds, Angles, Dihedrals, One_fives

# TODO need to think about where these canon functions should live in
def canon_bond(types, ind=None):
    """
    # Description: returns a canonicallized TAFFI bond. TAFFI bonds are written so that the lesser *atom_type* between 1 and 2 is first.
    #
    # inputs:      types: a list of taffi atom types defining the bond
    #              ind:   a list of indices corresponding to the bond
    #
    # returns:     a canonically ordered bond (and list of indices if ind was supplied)
    """


    # consistency checks
    if len(types) != 2:
        raise MoleculeException("ERROR in canon_bond: the supplied dihedral doesn't have two elements. Exiting...")
    if ind != None and len(ind) != 2:
        raise MoleculeException("ERROR in canon_bond: the iterable supplied to ind doesn't have two elements. Exiting...")

    # bond types are written so that the lesser *atom_type* between 1 and 2 is first.
    if types[0] <= types[1]:
        if ind == None:
            return types
        else:
            return types, ind
    else:
        if ind == None:
            return types[::-1]
        else:
            return types[::-1], ind[::-1]

def canon_angle(types, ind=None):
    """
     Description: returns a canonicallized TAFFI angle. TAFFI angles are written so that the lesser *atom_type* between 1 and 3 is first.
    #
    # inputs:      types: a list of taffi atom types defining the angle
    #              ind:   a list of indices corresponding to the angle
    #
    # returns:     a canonically ordered angle (and list of indices if ind was supplied)
    """

    # consistency checks
    if len(types) != 3:
        raise MoleculeException("ERROR in canon_angle: the supplied dihedral doesn't have three elements. Exiting...")
    if ind != None and len(ind) != 3:
        raise MoleculeException("ERROR in canon_angle: the iterable supplied to ind doesn't have three elements. Exiting...")

    # angle types are written so that the lesser *atom_type* between 1 and 3 is first.
    if types[0] <= types[2]:
        if ind == None:
            return types
        else:
            return types, ind
    else:
        if ind == None:
            return types[::-1]
        else:
            return types[::-1], ind[::-1]

def canon_dihedral(types_0, ind=None):
    """
    # Description: returns a canonicallized TAFFI dihedral. TAFFI dihedrals are written so that the lesser *atom_type* between 1 and 4 is first.
    #              In the event that 1 and 4 are of the same type, then the lesser of 2 and 3 goes first.
    #
    # inputs:      types: a list of taffi atom types defining the dihedral
    #              ind:   a list of indices corresponding to the dihedral
    #
    # returns:     a canonically ordered dihedral (and list of indices if ind was supplied)
    """

    # consistency checks
    if len(types_0) < 4:
        raise MoleculeException("ERROR in canon_dihedral: the supplied dihedral has less than four elements. Exiting...")
    if ind != None and len(ind) != 4:
        raise MoleculeException("ERROR in canon_dihedral: the iterable supplied to ind doesn't have four elements. Exiting...")

    # Grab the types and style component (the fifth element if available)
    types = list(types_0[:4])
    if len(types_0) > 4:
        style = [types_0[4]]
    else:
        style = []

    # dihedral types are written so that the lesser *atom_type* between 1 and 4 is first.
    # In the event that 1 and 4 are of the same type, then the lesser of 2 and 3 goes first
    if types[0] == types[3]:
        if types[1] <= types[2]:
            if ind == None:
                return tuple(types + style)
            else:
                return tuple(types + style), ind
        else:
            if ind == None:
                return tuple(types[::-1] + style)
            else:
                return tuple(types[::-1] + style), ind[::-1]
    elif types[0] < types[3]:
        if ind == None:
            return tuple(types + style)
        else:
            return tuple(types + style), ind
    else:
        if ind == None:
            return tuple(types[::-1] + style)
        else:
            return tuple(types[::-1] + style), ind[::-1]

def canon_improper(types, ind=None):
    """
    # Description: returns a canonicallized TAFFI improper. TAFFI impropers are written so that
    #              the three peripheral *atom_types* are written in increasing order.
    #
    # inputs:      types: a list of taffi atom types defining the improper
    #              ind:   a list of indices corresponding to the improper
    #
    # returns:     a canonically ordered improper (and list of indices if ind was supplied)
    """

    # consistency checks
    if len(types) != 4:
        raise MoleculeException("ERROR in canon_improper: the supplied improper doesn't have four elements. Exiting...")
    if ind != None and len(ind) != 4:
        raise MoleculeException("ERROR in canon_improper: the iterable supplied to ind doesn't have four elements. Exiting...")

    # improper types are written so that the lesser *atom_type* between 1 and 4 is first.
    # In the event that 1 and 4 are of the same type, then the lesser of 2 and 3 goes first
    if ind == None:
        return tuple([types[0]] + sorted(types[1:]))
    else:
        tmp_types, tmp_ind = zip(*sorted(zip(types[1:], ind[1:])))
        return tuple([types[0]] + list(tmp_types[:])), tuple([ind[0]] + list(tmp_ind[:]))

def remove_duplicate_modes(types, modes):
    """
    # Helper function for Find_modes_tmp that removes duplicate modes
    """
    a = [tuple(sorted(_)) for _ in modes]
    ind = [a.index(_) for _ in set(a)]
    return [types[_] for _ in ind], [modes[_] for _ in ind]


def main():
    """
    testing if the structure is passed down to all
    """
    molecule = Molecule(2)
    molecule.parse_from_xyz('test.xyz')
    print(molecule.AdjMat.elements)
    print(molecule.AtomType.elements)


if __name__ == '__main__':
    main()