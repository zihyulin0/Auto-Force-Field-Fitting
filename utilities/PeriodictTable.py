class Tables():
    @property
    def Radii(self):
        # Initialize UFF bond radii (Rappe et al. JACS 1992)
        # NOTE: Units of angstroms
        # NOTE: These radii neglect the bond-order and electronegativity corrections in the original paper. Where several values exist for the same atom, the largest was used.
        Radii = {'H': 0.354, 'He': 0.849, \
                 'Li': 1.336, 'Be': 1.074, 'B': 0.838, 'C': 0.757, 'N': 0.700, 'O': 0.658, 'F': 0.668, 'Ne': 0.920, \
                 'Na': 1.539, 'Mg': 1.421, 'Al': 1.244, 'Si': 1.117, 'P': 1.117, 'S': 1.064, 'Cl': 1.044, 'Ar': 1.032, \
                 'K': 1.953, 'Ca': 1.761, 'Sc': 1.513, 'Ti': 1.412, 'V': 1.402, 'Cr': 1.345, 'Mn': 1.382, 'Fe': 1.335,
                 'Co': 1.241, 'Ni': 1.164, 'Cu': 1.302, 'Zn': 1.193, 'Ga': 1.260, 'Ge': 1.197, 'As': 1.211, 'Se': 1.190,
                 'Br': 1.192, 'Kr': 1.147, \
                 'Rb': 2.260, 'Sr': 2.052, 'Y': 1.698, 'Zr': 1.564, 'Nb': 1.473, 'Mo': 1.484, 'Tc': 1.322, 'Ru': 1.478,
                 'Rh': 1.332, 'Pd': 1.338, 'Ag': 1.386, 'Cd': 1.403, 'In': 1.459, 'Sn': 1.398, 'Sb': 1.407, 'Te': 1.386,
                 'I': 1.382, 'Xe': 1.267, \
                 'Cs': 2.570, 'Ba': 2.277, 'La': 1.943, 'Hf': 1.611, 'Ta': 1.511, 'W': 1.526, 'Re': 1.372, 'Os': 1.372,
                 'Ir': 1.371, 'Pt': 1.364, 'Au': 1.262, 'Hg': 1.340, 'Tl': 1.518, 'Pb': 1.459, 'Bi': 1.512, 'Po': 1.500,
                 'At': 1.545, 'Rn': 1.42, \
                 'default': 0.7}

        # SAME AS ABOVE BUT WITH A SMALLER VALUE FOR THE Al RADIUS ( I think that it tends to predict a bond where none are expected
        Radii = {'H': 0.39, 'He': 0.849, \
                 'Li': 1.336, 'Be': 1.074, 'B': 0.838, 'C': 0.757, 'N': 0.700, 'O': 0.658, 'F': 0.668, 'Ne': 0.920, \
                 'Na': 1.539, 'Mg': 1.421, 'Al': 1.15, 'Si': 1.050, 'P': 1.117, 'S': 1.064, 'Cl': 1.044, 'Ar': 1.032, \
                 'K': 1.953, 'Ca': 1.761, 'Sc': 1.513, 'Ti': 1.412, 'V': 1.402, 'Cr': 1.345, 'Mn': 1.382, 'Fe': 1.335,
                 'Co': 1.241, 'Ni': 1.164, 'Cu': 1.302, 'Zn': 1.193, 'Ga': 1.260, 'Ge': 1.197, 'As': 1.211, 'Se': 1.190,
                 'Br': 1.192, 'Kr': 1.147, \
                 'Rb': 2.260, 'Sr': 2.052, 'Y': 1.698, 'Zr': 1.564, 'Nb': 1.473, 'Mo': 1.484, 'Tc': 1.322, 'Ru': 1.478,
                 'Rh': 1.332, 'Pd': 1.338, 'Ag': 1.386, 'Cd': 1.403, 'In': 1.459, 'Sn': 1.398, 'Sb': 1.407, 'Te': 1.386,
                 'I': 1.382, 'Xe': 1.267, \
                 'Cs': 2.570, 'Ba': 2.277, 'La': 1.943, 'Hf': 1.611, 'Ta': 1.511, 'W': 1.526, 'Re': 1.372, 'Os': 1.372,
                 'Ir': 1.371, 'Pt': 1.364, 'Au': 1.262, 'Hg': 1.340, 'Tl': 1.518, 'Pb': 1.459, 'Bi': 1.512, 'Po': 1.500,
                 'At': 1.545, 'Rn': 1.42, \
                 'default': 0.7}
        return Radii

    @property
    def MAXBONDS(self):
        return {'H': 2, 'He': 1, \
                     'Li': None, 'Be': None, 'B': 4, 'C': 4, 'N': 4, 'O': 2, 'F': 1, 'Ne': 1, \
                     'Na': None, 'Mg': None, 'Al': 4, 'Si': 4, 'P': None, 'S': None, 'Cl': 1, 'Ar': 1, \
                     'K': None, 'Ca': None, 'Sc': None, 'Ti': None, 'V': None, 'Cr': None, 'Mn': None, 'Fe': None,
                     'Co': None, 'Ni': None, 'Cu': None, 'Zn': None, 'Ga': 3, 'Ge': None, 'As': None, 'Se': None,
                     'Br': 1, 'Kr': None, \
                     'Rb': None, 'Sr': None, 'Y': None, 'Zr': None, 'Nb': None, 'Mo': None, 'Tc': None, 'Ru': None,
                     'Rh': None, 'Pd': None, 'Ag': None, 'Cd': None, 'In': None, 'Sn': None, 'Sb': None, 'Te': None,
                     'I': 1, 'Xe': None, \
                     'Cs': None, 'Ba': None, 'La': None, 'Hf': None, 'Ta': None, 'W': None, 'Re': None, 'Os': None,
                     'Ir': None, 'Pt': None, 'Au': None, 'Hg': None, 'Tl': None, 'Pb': None, 'Bi': None, 'Po': None,
                     'At': None, 'Rn': None}
    @property
    def MASSDICT(self):
        return {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 'C': 12.011,
                     'N': 14.00674, 'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, \
                     'Na': 22.989768, 'Mg': 24.3050, 'Al': 26.981539, 'Si': 28.0855, 'P': 30.973762,
                     'S': 32.066, 'Cl': 35.4527, 'Ar': 39.948, \
                     'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955910, 'Ti': 47.867, 'V': 50.9415,
                     'Cr': 51.9961, 'Mn': 54.938049, 'Fe': 55.845, 'Co': 58.933200, 'Ni': 58.6934,
                     'Cu': 63.546, 'Zn': 65.39, \
                     'Ga': 69.723, 'Ge': 72.61, 'As': 74.92159, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.80, \
                     'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585, 'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94,
                     'Tc': 98.0, 'Ru': 101.07, 'Rh': 102.90550, 'Pd': 106.42, 'Ag': 107.8682,
                     'Cd': 112.411, \
                     'In': 114.818, 'Sn': 118.710, 'Sb': 121.760, 'Te': 127.60, 'I': 126.90447,
                     'Xe': 131.29, \
                     'Cs': 132.90545, 'Ba': 137.327, 'La': 138.9055, 'Hf': 178.49, 'Ta': 180.9479,
                     'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217, 'Pt': 195.078,
                     'Au': 196.96655, 'Hg': 200.59, \
                     'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.98038, 'Po': 209.0, 'At': 210.0, 'Rn': 222.0}
    @property
    def ELEMENTDICT(self):
        return {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', \
                               11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', \
                               19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co',
                               28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', \
                               37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh',
                               46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', \
                               55: 'Cs', 56: 'Ba', 57: 'La', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir',
                               78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn'}

    @property
    def PERIODICT(self):
        return {"h": 1, "he": 2, \
                     "li": 3, "be": 4, "b": 5, "c": 6, "n": 7, "o": 8, "f": 9, "ne": 10, \
                     "na": 11, "mg": 12, "al": 13, "si": 14, "p": 15, "s": 16, "cl": 17, "ar": 18, \
                     "k": 19, "ca": 20, "sc": 21, "ti": 22, "v": 23, "cr": 24, "mn": 25, "fe": 26,
                     "co": 27, "ni": 28, "cu": 29, "zn": 30, "ga": 31, "ge": 32, "as": 33, "se": 34,
                     "br": 35, "kr": 36, \
                     "rb": 37, "sr": 38, "y": 39, "zr": 40, "nb": 41, "mo": 42, "tc": 43, "ru": 44,
                     "rh": 45, "pd": 46, "ag": 47, "cd": 48, "in": 49, "sn": 50, "sb": 51, "te": 52,
                     "i": 53, "xe": 54, \
                     "cs": 55, "ba": 56, "hf": 72, "ta": 73, "w": 74, "re": 75, "os": 76, "ir": 77,
                     "pt": 78, "au": 79, "hg": 80, "tl": 81, "pb": 82, "bi": 83, "po": 84, "at": 85,
                     "rn": 86}
    @property
    def ELETRONEGATIVITY(self):
        return {"h": 2.3, "he": 4.16, \
                             "li": 0.91, "be": 1.58, "b": 2.05, "c": 2.54, "n": 3.07, "o": 3.61, "f": 4.19, "ne": 4.79, \
                             "na": 0.87, "mg": 1.29, "al": 1.61, "si": 1.91, "p": 2.25, "s": 2.59, "cl": 2.87,
                             "ar": 3.24, \
                             "k": 0.73, "ca": 1.03, "sc": 1.19, "ti": 1.38, "v": 1.53, "cr": 1.65, "mn": 1.75,
                             "fe": 1.80, "co": 1.84, "ni": 1.88, "cu": 1.85, "zn": 1.59, "ga": 1.76, "ge": 1.99,
                             "as": 2.21, "se": 2.42, "br": 2.69, "kr": 2.97, \
                             "rb": 0.71, "sr": 0.96, "y": 1.12, "zr": 1.32, "nb": 1.41, "mo": 1.47, "tc": 1.51,
                             "ru": 1.54, "rh": 1.56, "pd": 1.58, "ag": 1.87, "cd": 1.52, "in": 1.66, "sn": 1.82,
                             "sb": 1.98, "te": 2.16, "i": 2.36, "xe": 2.58, \
                             "cs": 0.66, "ba": 0.88, "la": 1.09, "hf": 1.16, "ta": 1.34, "w": 1.47, "re": 1.60,
                             "os": 1.65, "ir": 1.68, "pt": 1.72, "au": 1.92, "hg": 1.76, "tl": 1.79, "pb": 1.85,
                             "bi": 2.01, "po": 2.19, "at": 2.39, "rn": 2.60}
    @property
    def LONE_ELETRON(self):
        return {'h': 0, 'he': 2, \
                                 'li': 0, 'be': 2, 'b': 0, 'c': 0, 'n': 2, 'o': 4, 'f': 6, 'ne': 8, \
                                 'na': 0, 'mg': 2, 'al': 0, 'si': 0, 'p': 2, 's': 4, 'cl': 6, 'ar': 8, \
                                 'k': 0, 'ca': 2, 'sc': None, 'ti': None, 'v': None, 'cr': None, 'mn': None, 'fe': None,
                                 'co': None, 'ni': None, 'cu': None, 'zn': None, 'ga': 10, 'ge': 0, 'as': 3, 'se': 4,
                                 'br': 6, 'kr': None, \
                                 'rb': 0, 'sr': 2, 'y': None, 'zr': None, 'nb': None, 'mo': None, 'tc': None,
                                 'ru': None, 'rh': None, 'pd': None, 'ag': None, 'cd': None, 'in': None, 'sn': None,
                                 'sb': None, 'te': None, 'i': 6, 'xe': None, \
                                 'cs': 0, 'ba': 2, 'la': None, 'hf': None, 'ta': None, 'w': None, 're': None,
                                 'os': None, 'ir': None, 'pt': None, 'au': None, 'hg': None, 'tl': None, 'pb': None,
                                 'bi': None, 'po': None, 'at': None, 'rn': None}
    @property
    def POLARIZABILITY(self):
        return {"h": 4.5, "he": 1.38, \
                              "li": 164.0, "be": 377, "b": 20.5, "c": 11.3, "n": 7.4, "o": 5.3, "f": 3.74, "ne": 2.66, \
                              "na": 163.0, "mg": 71.2, "al": 57.8, "si": 37.3, "p": 25.0, "s": 19.4, "cl": 14.6,
                              "ar": 11.1, \
                              "k": 290.0, "ca": 161.0, "sc": 97.0, "ti": 100.0, "v": 87.0, "cr": 83.0, "mn": 68.0,
                              "fe": 62.0, "co": 55, "ni": 49, "cu": 47.0, "zn": 38.7, "ga": 50.0, "ge": 40.0,
                              "as": 30.0, "se": 29.0, "br": 21.0, "kr": 16.8, \
                              "rb": 320.0, "sr": 197.0, "y": 162, "zr": 112.0, "nb": 98.0, "mo": 87.0, "tc": 79.0,
                              "ru": 72.0, "rh": 66, "pd": 26.1, "ag": 55, "cd": 46.0, "in": 65.0, "sn": 53.0,
                              "sb": 43.0, "te": 28.0, "i": 32.9, "xe": 27.3, }
    @property
    def BOND_ENERGY(self):
        """
        Bond energy dictionary {}-{}-{} refers to atom1, atom2 additional bonds number (1 refers to double bonds)
        If energy for multiple bonds is missing, it means it's unusual to form multiple bonds, such value will be -10000.0, if energy for single bonds if missing, directly take multiple bonds energy as the difference
        From https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Chemical_Bonding/Fundamentals_of_Chemical_Bonding/Bond_Energies
        find_lewis.be = { "6-6-1": 267, "6-6-2":492, "6-7-1":310, "6-7-2":586, "6-8-1":387, "6-8-2":714, "7-8-1":406, "7-7-1":258, "7-7-2":781, "8-8-1":349, "8-16-1":523, "16-16-1":152}
        Or from https://www2.chemistry.msu.edu/faculty/reusch/OrgPage/bndenrgy.htm ("6-15-1" is missing)
        Remove 6-16-1:73
        """
        return {"6-6-1": 63, "6-6-2": 117, "6-7-1": 74, "6-7-2": 140, "6-8-1": 92.5, "6-8-2": 172.5,
                             "7-7-1": 70.6, "7-7-2": 187.6, "7-8-1": 88, "8-8-1": 84, "8-15-1": 20, "8-16-1": 6,
                             "15-15-1": 84, "15-15-2": 117, "15-16-1": 70}

    @property
    def SATURATION(self):
        return {'H': 1, 'He': 1, \
                                  'Li': 1, 'Be': 2, 'B': 3, 'C': 4, 'N': 3, 'O': 2, 'F': 1, 'Ne': 1, \
                                  'Na': 1, 'Mg': 2, 'Al': 3, 'Si': 4, 'P': 3, 'S': 2, 'Cl': 1, 'Ar': 1, \
                                  'K': 1, 'Ca': 2, 'Sc': None, 'Ti': None, 'V': None, 'Cr': None, 'Mn': None,
                                  'Fe': None, 'Co': None, 'Ni': None, 'Cu': None, 'Zn': None, 'Ga': None, 'Ge': None,
                                  'As': None, 'Se': None, 'Br': 1, 'Kr': None, \
                                  'Rb': 1, 'Sr': 2, 'Y': None, 'Zr': None, 'Nb': None, 'Mo': None, 'Tc': None,
                                  'Ru': None, 'Rh': None, 'Pd': None, 'Ag': None, 'Cd': None, 'In': None, 'Sn': None,
                                  'Sb': None, 'Te': None, 'I': 1, 'Xe': None, \
                                  'Cs': 1, 'Ba': 2, 'La': None, 'Hf': None, 'Ta': None, 'W': None, 'Re': None,
                                  'Os': None, 'Ir': None, 'Pt': None, 'Au': None, 'Hg': None, 'Tl': None, 'Pb': None,
                                  'Bi': None, 'Po': None, 'At': None, 'Rn': None}
