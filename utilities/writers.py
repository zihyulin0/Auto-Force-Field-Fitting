"""
Writer function for multiple file types
"""

def write_xyz(file_prefix, elements, geometry, additional_col=None):
    # Open file for writing and write header
    outfile_name = file_prefix + '.xyz'
    with open(outfile_name, 'w') as fid:
        fid.write('{}\n\n'.format(len(elements)))

        for count_i, i in enumerate(elements):
            fid.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} '.format(i, geometry[count_i, 0], geometry[count_i, 1],
                                                                            geometry[count_i, 2]))
            if additional_col is not None:
                fid.write('{}'.format(additional_col[count_i]))
            fid.write('\n')
    return outfile_name

def write_modelist(name,modes,bond_mat=[]):
    """
    # Wrapper function for write commands for *modes files
    """
    with open(name,'w') as f:
        for count_i,i in enumerate(modes.keys()):
            if len(i) == 1:
                modetype='atom'
            elif len(i) == 2:
                modetype='bond'
            elif len(i) == 3:
                modetype='angle'
            elif len(modes[i]["modes"][0]) == 4 and 2 in [ j[modes[i]["modes"][0][1],modes[i]["modes"][0][2]] for j in bond_mat ]:
                modetype='harmonic_dihedral'
            elif len(modes[i]["modes"][0]) == 4:
                modetype='dihedral'
            f.write("\n{} start\n".format(modetype))
            f.write('{}\n'.format(len(modes[i]["modes"])))
            f.write('{}\n'.format(" ".join([ "{:<60s}".format(str(j)) for j in modes[i]["gens"] ])))
            f.write("{}\n".format(" ".join([ "{:<60s}".format("_".join([ str(k) for k in j])) for j in modes[i]["modes"] ])))
            for count_j,j in enumerate(modes[i]["atomtypes"][0]):
                f.write('{}\n'.format(" ".join([ "{:<60s}".format(modes[i]["atomtypes"][k][count_j]) for k in range(len(modes[i]["atomtypes"])) ])))
            f.write("{} end\n".format(modetype))
    return
