import numpy as np
from copy import deepcopy
from utilities.PeriodictTable import Tables
from utilities.Excpetions import TAFFIException
from utilities.parse import read_alljson
from scipy.spatial.distance import cdist
import os, random
from model.structure import StructureBase

class AdjacencyException(TAFFIException):
    pass

class AdjacencyMatrix(StructureBase):
    """
    Operations related to adjacency matrix
    Main problem: you cannot guarantee the structure and the existing adj_mat is consistent
    It might generate one adj_mat with a set of geometries and later modified the geometry wihtout generate adj_mat again
    """
    def __init__(self, adjmat=None):
        super().__init__()
        # adj_mat needs to be protected to avoid being accidentally changed when being inherited
        self._adj_mat = adjmat

    @property
    def adj_mat(self):
        return self._adj_mat

    @adj_mat.setter
    def adj_mat(self, value):
        self._adj_mat = value

    def build_adj_mat(self):
        self.Table_generator()

    def graph_seps(self):
        """
        Returns a matrix of graphical separations for all nodes in a graph defined by the inputted adjacency matrix
        """
        if self.adj_mat is None:
            raise AdjacencyException("graph_seps called before adj_mat is generated")

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
        :returns: adjacency matrix
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

        return

    def return_connected(self, start=0, avoid=[]):
        """
        Returns the set of connected nodes to the start node, while avoiding any connections through nodes in the avoid list.
        """

        # Initialize the avoid list with the starting index
        avoid = set(avoid + [start])

        # new_0 holds the most recently encountered nodes, beginning with start
        # new_1 is a set holding all of the encountered nodes
        new_0 = [start]
        new_1 = set([start])

        # keep looping until no new nodes are encountered
        while len(new_0) > 0:
            # reinitialize new_0 with new connections
            new_0 = [count_j for i in new_0 for count_j, j in enumerate(self.adj_mat[i]) if j == 1 and count_j not in avoid]

            # update the new_1 set and avoid list with the most recently encountered new nodes
            new_1.update(new_0)
            avoid.update(new_0)

        # return the set of encountered nodes
        return new_1

    def Dijkstra(self, start=0, end=-1):
        """
        Description: This is a simple implementation of the Dijkstra algorithm for
                  finding the backbone of a polymer
        """

        # Default to the last node in Adj_mat if end is unassigned or less than 0
        if end < 0:
            end = len(self.adj_mat) + end

        # Remove terminal sites (sites with only a single length 2
        # self walk). Continue until all of the terminal structure
        # has been removed from the topology.
        Adj_trimmed = np.copy(self.adj_mat)

        # Initialize Distances, Previous, and Visited lists
        Distances = np.array([100000] * (len(self.adj_mat)))  # Holds shortest distance to origin from each site
        Distances[start] = 0  # Sets the separation of the initial node from the initial node to zero
        Previous = np.array([-1] * len(self.adj_mat))  # Holds the previous site on the short distance to origin
        Visited = [0] * len(self.adj_mat)  # Holds which sites have been visited

        # Initialize current site (i) and neighbors list
        i = start  # current site
        neighbors = []

        # Iterate through sites. At each step find the shortest distance between all of hte
        # current sites neighbors and the START. Update the shortest distances of all sites
        # and choose the next site to iterate on based on which has the shortest distance of
        # among the UNVISITED set of sites. Once the terminal site is identified the shortest
        # path has been found
        while (0 in Visited):

            # If the current site is the terminal site, then the algorithm is finished
            if i == end:
                break

            # Add new neighbors to the list
            neighbors = [count_j for count_j, j in enumerate(Adj_trimmed[i]) if j == 1]

            # Remove the current site from the list of unvisited
            Visited[i] = 1

            # Iterate over neighbors and update shortest paths
            for j in neighbors:

                # Update distances for current generation of connections
                if Distances[i] + Adj_trimmed[j, i] < Distances[j]:
                    Distances[j] = Distances[i] + Adj_trimmed[j, i]
                    Previous[j] = i

            # Find new site based on the minimum separation (only go to unvisited sites!)
            tmp = min([j for count_j, j in enumerate(Distances) if Visited[count_j] == 0])
            i = [count_j for count_j, j in enumerate(Distances) if j == tmp and Visited[count_j] == 0][0]

        # Find shortest path by iterating backwards
        # starting with the end site.
        Shortest_path = [end]
        i = end
        while (i != start):
            Shortest_path = Shortest_path + [Previous[i]]
            i = Previous[i]

        # Reverse order of the list to go from start to finish
        Shortest_path = Shortest_path[::-1]
        return Shortest_path
