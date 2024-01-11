import os, json, codecs
from Excpetions import TAFFIException
import numpy as np

class ReadingException(TAFFIException):
    pass

def read_alljson(jsonfile):
    if os.path.isfile(jsonfile) is False:
        raise ReadingException("Error: json file: {} not found".format(jsonfile))
    obj_text = codecs.open(jsonfile, 'r').read()
    jload = json.loads(obj_text)
    return jload
def read_n_column(infile, n, skip_n=0):
    """
    read a file with n columns, this will skip any empty lines or line start with #

    :type infile: str
    :param infile: input file name

    :type n: int
    :param n: number of columns

    :type skip_n: int
    :param skip_n: skip n lines

    :rtype: numpy array
    :return: a numpy array with n columns
    """
    data = []
    with open(infile, "r") as file:
        for lc, line in enumerate(file):
            if lc < skip_n: continue
            line = line.strip()  # Remove leading/trailing whitespace and newline characters
            if not line or line.startswith("#"):
                continue  # Skip empty lines and lines starting with "#"
            fields = line.split()
            if len(fields) >= n:
                data.append(fields[:n])
                if len(fields) > n:
                    print("WARNING: line has more columns than the specified readin # of columns")
            else:
                raise ReadingException('line has fewer columns than the specified reading # of columns')
    return np.array(data)

def read_n_lines(infile, n):
    """
    Reading the first n line of the input file

    :type infile: str
    :param infile: input file to read

    :type n: int
    :param n: number of lines to read

    :rtype: list
    :return: a list of raw lines (inclue white space or whatsoever)
    """
    data = []
    with open(infile, "r") as file:
        for lc, line in enumerate(file):
            if lc < n: data.append(line)
            else: break
    return data
def xyz_parse(infile, read_types=False, q_opt=False):
    """
    Simple wrapper function for grabbing the coordinates and
    elements from a xyz file

    :type infile: str
    :param infile: filename of the xyz

    :type read_types: bool
    :param read_types: whether to read atom types from fifth column

    :type q_opt: bool
    :param q_opt: whether to read total charge of molecule from line 1

    :rtype: (list, array, int)
    :return: elements, geometry arrays, total charge of the molecule
    """
    # Iterate over the remainder of contents and read the
    # geometry and elements into variable. Note that the
    # first two lines are considered a header

    lines = read_n_lines(infile, 2)

    # reading number of atoms and initialize
    fields = lines[0].split()
    if len(fields) < 1:
        raise ReadingException("ERROR in xyz_parse: {} is missing atom number\
              information".format(infile))
    N_atoms = int(fields[0])

    # Get charge
    if q_opt:
        fields = lines[1].split()
        q = 0
        if "q" in fields:
            try:
                q = int(fields[fields.index("q") + 1])
            except:
                print("Charge specification misformatted in {}.\
                      Defaulting to q=0.".format(infile))
    # parse body
    n_column = 4 if not read_types else 5
    data = read_n_column(infile, n_column, skip_n=2)
    # Consistency check
    if data.shape[0] != N_atoms:
        raise ReadingException("ERROR in xyz_parse: {} atom number in header doesn't match coordinates.".format(infile))

    Elements = data[:,0].astype(str)
    Geometry = data[:,range(1,4)].astype(float)
    if read_types:
        Atom_types = data[:, 4]

    if q_opt:
        return Elements, Geometry, q
    else:
        return Elements, Geometry


# Description: Parses taffi.db files and returns a dictionary with the parameters and modes
def parse_FF_params(FF_files, FF_dict={"masses": {}, "charges": {}, "bonds": {}, "angles": {}, "dihedrals": {},
                                       "dihedrals_harmonic": {}, "vdw": {}}):
    modes_from_FF = []
    for i in FF_files:
        with open(i, 'r') as f:
            for lines in f:
                fields = lines.split()
                if len(fields) == 0: continue
                if fields[0].lower() == "atom":   FF_dict["masses"][fields[1]] = float(fields[3])
                if fields[0].lower() == "charge": FF_dict["charges"][fields[1]] = float(fields[2])
                if fields[0].lower() == "bond":
                    modes_from_FF += [(fields[1], fields[2])]
                    modes_from_FF += [(fields[2], fields[1])]
                    FF_dict["bonds"][(fields[1], fields[2])] = [fields[3], float(fields[4]), float(fields[5])]
                if fields[0].lower() == "angle":
                    modes_from_FF += [(fields[1], fields[2], fields[3])]
                    modes_from_FF += [(fields[3], fields[2], fields[1])]
                    FF_dict["angles"][(fields[1], fields[2], fields[3])] = [fields[4], float(fields[5]),
                                                                            float(fields[6])]
                if fields[0].lower() in ["dihedral", "torsion"]:
                    modes_from_FF += [(fields[1], fields[2], fields[3], fields[4])]
                    modes_from_FF += [(fields[4], fields[3], fields[2], fields[1])]
                    if fields[5] == "opls":
                        FF_dict["dihedrals"][(fields[1], fields[2], fields[3], fields[4])] = [fields[5]] + [float(i) for
                                                                                                            i in fields[
                                                                                                                 6:10]]
                    elif fields[5] == "harmonic":
                        FF_dict["dihedrals_harmonic"][(fields[1], fields[2], fields[3], fields[4])] = [fields[5]] + [
                            float(fields[6]), int(float(fields[7])), int(float(fields[8]))]
                    elif fields[5] == "quadratic":
                        FF_dict["dihedrals_harmonic"][(fields[1], fields[2], fields[3], fields[4])] = [fields[5]] + [
                            float(fields[6]), float(fields[7])]
                if fields[0].lower() == "vdw":
                    FF_dict["vdw"][(fields[1], fields[2])] = [fields[3], float(fields[4]), float(fields[5])]
                    FF_dict["vdw"][(fields[2], fields[1])] = [fields[3], float(fields[4]), float(fields[5])]

    return FF_dict, modes_from_FF

def main():
    elements, geometry  = xyz_parse('test.xyz')
    print(elements)
    print(geometry)


if __name__ == '__main__':
    main()