#!/bin/env python
# Author: Brett Savoie (brettsavoie@gmail.com)

# Add TAFFY Lib to path
import sys, os
from subprocess import PIPE, Popen
from model.adjacency import AdjacencyMatrix
from model.molecule import Molecule
from model.atomtypes import AtomTypesBaseOperation
from utilities.Excpetions import TAFFIException
import numpy as np
from copy import deepcopy


class TransifyException(TAFFIException):
    pass


# TODO I think this might be fine as it is for now, but we might want to change that

def transify(geo, adj_mat, start=-1, end=-1, elements=None, opt_terminals=False, opt_final=False):
    # TODO how opt_geo depends on a lot of thin makes me want to move this
    """
    This function takes a geometry and adjacency matrix and returns an "all-trans" aligned geometric conformer
    have a recursive call in here
    """

    # Sanity checks
    if start < -1 or start >= len(geo):
        raise TransifyException("ERROR in transify: start must be set to an index "
                                "within geo or -1. Exiting...")
    if end < -1 or end >= len(geo):
        raise TransifyException("ERROR in transify: end must be set to an index "
                                "within geo or -1. Exiting...")
    if len(geo) != len(adj_mat):
        raise TransifyException("ERROR in transify: geo and adj_mat must have the same "
                                "dimensions. Exiting...")
    if opt_terminals and elements is None:
        raise TransifyException("ERROR in transify: element information "
                                "must be supplied for terminals to be "
                                "optimized. Exiting...")
    if opt_terminals and elements is not None and len(elements) != len(geo):
        raise TransifyException("ERROR in transify: length of elements is not equal to the length of geo. Exiting...")
    if opt_final and elements is None:
        raise TransifyException("ERROR in transify: element information must be supplied for the final structure to "
                                "be optimized. Exiting...")
    if opt_final and elements is not None and len(elements) != len(geo):
        raise TransifyException("ERROR in transify: length of elements is not equal to the length of geo. Exiting...")

    # Return the geometry if it is empty
    if len(geo) == 0:
        return geo

    AdjMat = AdjacencyMatrix(adjmat=adj_mat)
    AdjMat.parse_data(geometry=geo)
    atype_op = AtomTypesBaseOperation(AdjMat=AdjMat)

    # Calculate the pair-wise graphical separations between all atoms
    seps = AdjMat.graph_seps()

    # Find the most graphically separated pair of atoms. 
    if start == -1 and end == -1:
        max_ind = np.where(seps == seps.max())
        start, end = max_ind[0][0], max_ind[1][0]

    # Find the most graphically separated atom from the start atom
    elif end == -1:
        end = np.argmax(seps[start])

        # Find the shortest pathway between these points
    pathway = AdjMat.Dijkstra(start, end)

    # If the molecule doesn't have any dihedrals then return
    if len(pathway) < 4:
        return geo

    # Initialize the list of terminal atoms and ring atoms (used in a couple places so initialized here relatively early)
    terminals = set([count_i for count_i, i in enumerate(adj_mat) if np.sum(i) == 1])
    rings = set([count_i for count_i, i in enumerate(geo) if atype_op.ring_atom_new(count_i) == True])

    # Work through the backbone dihedrals and straighten them out
    for i in range(1, len(pathway) - 2):

        if pathway[i] in rings: continue

        # Collect the atoms that are connected to the 2 atom of the dihedral but not to the 3 atom of the dihedral (and vice versus)
        group_1 = AdjMat.return_connected(start=pathway[i], avoid=[
            pathway[i + 1]])  # not used since only the forward portion of the pathway gets rotated
        group_2 = AdjMat.return_connected(start=pathway[i + 1], avoid=[pathway[i]])

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
                branch_ind = [conn_terminals[0]] + list(AdjMat.return_connected(start=i,
                                                                                avoid=[count_j for count_j, j in
                                                                                       enumerate(adj_mat[i]) if
                                                                                       j == 1 and count_j in pathway] + conn_terminals))
                geo[branch_ind[1:]] = transify(geo[branch_ind], adj_mat[branch_ind, :][:, branch_ind], start=0)[1:]

            # If no terminal connections exists then the first dihedral of the branch is not adjusted. The logic is that the sp2/sp3 alignment is better judged by the initial guess. 
            else:
                branch_ind = list(AdjMat.return_connected(start=i,
                                                          avoid=[count_j for count_j, j in enumerate(adj_mat[i]) if
                                                                 j == 1 and count_j in pathway] + conn_terminals))
                geo[branch_ind] = transify(geo[branch_ind], adj_mat[branch_ind, :][:, branch_ind],
                                           start=next(count_j for count_j, j in enumerate(branch_ind) if j == i))

    # If optimize structure is set to True then the transified structure is relaxed using obminimize
    if opt_final :
        opt_geo(geo, adj_mat, elements, ff='uff')

    # If optimize terminals is set to true then a conformer search is performed over the terminal groups
    if opt_terminals:
        geo = opt_terminal_centers(geo, adj_mat, elements, ff='uff')

    return geo



def opt_geo(geo, adj_mat, elements, ff='uff'):
    """
     Description: This function calls obminimize (open babel geometry optimizer function) to optimize the current geometry

    :param geo:      Nx3 array of atomic coordinates
    :param adj_mat:  NxN array of connections
    :param elements: N list of element labels
    :param ff:       force-field specification passed to obminimize (uff, gaff)

     :returns:     geo:      Nx3 array of optimized atomic coordinates

    """
    # Write a temporary molfile for obminimize to use
    tmp_filename = '.tmp.mol'
    count = 0
    while os.path.isfile(tmp_filename):
        count += 1
        if count == 10:
            raise TransifyException("ERROR in opt_geo: could not find a suitable filename for the tmp geometry. Exiting...")
        else:
            tmp_filename = ".tmp" + tmp_filename

            # Use the mol_write function imported from the write_functions.py
    # to write the current geometry and topology to file
    tmp_mol = Molecule(2)
    tmp_mol.parse_data(elements=elements, geometry=geo, adj_mat= adj_mat)
    tmp_mol.mol_write(tmp_filename, append_opt=False)

    # Popen(stdout=PIPE).communicate() returns stdout,stderr as strings
    output = \
    Popen('obminimize -sd -c 1e-20 -n 100000 -ff {} {}'.format(ff, tmp_filename).split(), stdin=PIPE, stdout=PIPE,
          stderr=PIPE, bufsize=-1).communicate()[0]
    output = str(output)

    # Parse the optimized *.xyz geometry
    counter = 0
    for count_i, i in enumerate(output.split(r'\n')[2:]):
        fields = i.split()
        if len(fields) != 4:
            continue
        else:
            geo[counter] = [float(j) for j in fields[1:]]
            counter += 1

    # Remove the tmp file that was read by obminimize
    try:
        os.remove(tmp_filename)
    except:
        pass
    return geo

def opt_terminal_centers(geo, adj_mat, elements, ff='uff'):
    """
    This function instantiates all conformers involving adjustments of terminal centers (non-terminal atoms with
    only one non-terminal connection, such as a methyl carbon). The function performs repeated calls to obminimize
    (open babel geometry optimizer function) to minimize each conformer and returns the lowest energy geometry.

    :param geo:      Nx3 array of atomic coordinates
    :param adj_mat:  NxN array of connections
    :param elements: N list of element labels
    :param ff:       force-field specification passed to obminimize (uff, gaff)

    returns:     geo:      Nx3 array of optimized atomic coordinates
#
    """
    # Generate lists of rotation centers, the terminal and non-terminal atoms connected to each center and the rotation vectors for each center
    centers = terminal_centers(adj_mat)
    terminal_cons = [[count_j for count_j, j in enumerate(adj_mat[i]) if j == 1 and np.sum(adj_mat[count_j]) == 1] for i
                     in centers]
    nonterminal_cons = [
        next(count_j for count_j, j in enumerate(adj_mat[i]) if j == 1 and np.sum(adj_mat[count_j]) != 1) for i in
        centers]
    rot_vecs = [geo[nonterminal_cons[count_i]] - geo[i] for count_i, i in enumerate(centers)]

    # Print a warning if there are a large number of centers
    if len(centers) > 5:
        print(
            "WARNING in opt_terminal_centers: Iterating through all terminal center conformers. There are {} centers, so this might take a while.".format(
                len(centers)))

        # Collect the rotation operations
    ops = []
    for i in terminal_cons:

        if len(i) == 1:
            ops += [[0.0]]
        elif len(i) == 2 and False in [elements[i[0]] == elements[j] for j in i]:
            ops += [[j * 90.0 for j in range(4)]]
        elif len(i) == 2:
            ops += [[j * 90.0 for j in range(2)]]
        elif len(i) == 3 and False in [elements[i[0]] == elements[j] for j in i]:
            ops += [[j * 60.0 for j in range(6)]]
        elif len(i) == 3:
            ops += [[j * 60.0 for j in range(2)]]

    # Loop over all combinations of rotation operations and retain the lowest energy conformer
    E_min = None
    geo_min = None
    min_num = None
    for count_i, i in enumerate(combos(ops)):
        new_geo = deepcopy(geo)
        for count_j, j in enumerate(i):
            for k in terminal_cons[count_j]:
                new_geo[k] = axis_rot(new_geo[k], rot_vecs[count_j], v2=new_geo[centers[count_j]], theta=j)

        # UNCOMMENT TO KEEP ALL CONFORMERS: For debugging it can be useful to see the individual conformers.
        # mol_write("appended.mol",elements,new_geo,adj_mat,append_opt=True)

        # Optimize the current geometry
        new_geo = opt_geo(new_geo, adj_mat, elements, ff='uff')

        # Generate a temporary filename for opt_terminal_centers to use
        tmp_filename = '.tmp.mol'
        count = 0
        while os.path.isfile(tmp_filename):
            count += 1
            if count == 10:
                print(
                    "ERROR in opt_terminal_centers: could not find a suitable filename for the tmp geometry. Exiting...")
                return geo
            else:
                tmp_filename = ".tmp" + tmp_filename

                # Calculate the energy of the optimized conformer
        # use os.devnull to suppress stderr from the obminimize call
        tmp_mol = Molecule()
        tmp_mol.parse_data(elements=elements, geometry=new_geo, adj_mat=adj_mat)
        tmp_mol.mol_write(tmp_filename, append_opt=True)

        # Popen(stdout=PIPE).communicate() returns stdout,stderr as strings
        output = Popen('obenergy -ff {} {}'.format(ff, tmp_filename).split(), stdin=PIPE, stdout=PIPE, stderr=PIPE,
                       bufsize=-1).communicate()[0]
        output = str(output)

        # Grab the energy from the output
        Energy = next(float(j.split()[3]) for j in output.split(r'\n') if
                      len(j.split()) > 4 and j.split()[0] == "TOTAL" and j.split()[1] == "ENERGY" and j.split()[
                          2] == "=")
        if E_min is None or Energy < E_min:
            E_min = Energy
            geo_min = deepcopy(new_geo)
            min_num = count_i
        os.remove(tmp_filename)

    # UNCOMMENT TO OPTIMIZE ALL CONFORMERS: For debugging it can be useful to see the individual optimized conformers.
    # with open('appended.opt.xyz','w') as f:
    #     with open(os.devnull,'w') as devnull:
    #         call('obminimize -sd -c 1e-10 -ff uff appended.mol'.split(),stdout=f,stderr=devnull)

    return geo_min


def terminal_centers(adj_mat):
    """
    # Returns the list of terminal centers, where a terminal center is an atom with only one
    # bond to a non-terminal atom.
    """
    return [ count_i for count_i,i in enumerate(adj_mat) if sum(i) > 1 and len([ j for count_j,j in enumerate(i) if j == 1 and sum(adj_mat[count_j]) > 1 ]) == 1 ]

def combos(x):
    """
    Generates all unique combinations (independent of order) of objects in a list of lists (x)

    :param x: a list of lists

    :returns: an iterable list of objects that yields a unique combination of objects (one from each sublist in x) from the x.
    """
    if len(x) == 0:
        yield []
    elif hasattr(x[0], '__iter__'):
        for i in x[0]:
            for j in combos(x[1:]):
                yield [i]+j
    else:
        for i in x:
            yield [i]


def dihedral_calc(xyz):
    """
    # Description: Calculates the dihedral angle (in radians) for
    #              a quadruplet of atoms
    #
    # Inputs:      xyz: a 4x3 numpy array, where each row holds the
    #                   cartesian position of each atom and row indices
    #                   determine placement in the dihedral ( i.e, which
    #                   atoms correspond to 1-2-3-4 in the dihedral )
    #
    # Returns:    angle: dihedral angle in radians
    """
    # Calculate the 2-1 vector
    v1 = (xyz[1] - xyz[0])

    # Calculate the 3-2 vector
    v2 = (xyz[2] - xyz[1])

    # Calculate the 4-3 bond vector
    v3 = (xyz[3] - xyz[2])

    # Calculate dihedral (funny np.arctan2 formula avoids the use of arccos and retains the sign of the angle)
    angle = np.arctan2(np.dot(v1, np.cross(v2, v3)) * (np.dot(v2, v2)) ** (0.5),
                       np.dot(np.cross(v1, v2), np.cross(v2, v3)))

    return angle


def axis_rot(Point, v1, v2=[0.0, 0.0, 0.0], theta=0.0, mode='angle'):
    """
    # Description:
    # Rotate Point by an angle, theta, about the vector with an orientation of v1 passing through v2.
    # Performs counter-clockwise rotations (i.e., if the direction vector were pointing
    # at the spectator, the rotations would appear counter-clockwise)
    # For example, a 90 degree rotation of a 0,0,1 about the canonical
    # y-axis results in 1,0,0.
    #
    # Point: 1x3 array, coordinates to be rotated
    # v1: 1x3 array, rotation direction vector
    # v2: 1x3 array, point the rotation passes through (default is the origin)
    # theta: scalar, magnitude of the rotation (defined by default in degrees) (default performs no rotation)
    """

    # Temporary variable for performing the transformation
    rotated = np.array([Point[0], Point[1], Point[2]])

    # If mode is set to 'angle' then theta needs to be converted to radians to be compatible with the
    # definition of the rotation vectors
    if mode == 'angle':
        theta = theta * np.pi / 180.0

    # Rotation carried out using formulae defined here (11/22/13) http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/)
    # Adapted for the assumption that v1 is the direction vector and v2 is a point that v1 passes through
    a = v2[0]
    b = v2[1]
    c = v2[2]
    u = v1[0]
    v = v1[1]
    w = v1[2]
    L = u ** 2 + v ** 2 + w ** 2

    # Rotate Point
    x = rotated[0]
    y = rotated[1]
    z = rotated[2]

    # x-transformation
    rotated[0] = (a * (v ** 2 + w ** 2) - u * (b * v + c * w - u * x - v * y - w * z)) \
                 * (1.0 - np.cos(theta)) + L * x * np.cos(theta) + L ** (0.5) * (
                         -c * v + b * w - w * y + v * z) * np.sin(theta)

    # y-transformation
    rotated[1] = (b * (u ** 2 + w ** 2) - v * (a * u + c * w - u * x - v * y - w * z)) \
                 * (1.0 - np.cos(theta)) + L * y * np.cos(theta) + L ** (0.5) * (
                         c * u - a * w + w * x - u * z) * np.sin(theta)

    # z-transformation
    rotated[2] = (c * (u ** 2 + v ** 2) - w * (a * u + b * v - u * x - v * y - w * z)) \
                 * (1.0 - np.cos(theta)) + L * z * np.cos(theta) + L ** (0.5) * (
                         -b * u + a * v - v * x + u * y) * np.sin(theta)

    rotated = rotated / L
    return rotated


def main(argv):
    return


if __name__ == "__main__":
    main(sys.argv[1:])
