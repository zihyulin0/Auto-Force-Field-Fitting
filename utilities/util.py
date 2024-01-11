import os, fnmatch
from utilities.writers import write_xyz
import subprocess
from subprocess import PIPE

def find_files(files,recursive=False):
    """
    wrapper for grabbing file locations, with wildcard support.


    """
    wc_files  = [ i for i in files if "*" in i ]
    files = [ i for i in files if "*" not in i ]
    if recursive:
        files += [ dp+"/"+f for i in wc_files for dp,dn,fn in os.walk('.') for f in fn if fnmatch.fnmatch(f,i) ]
    else:
        for i in wc_files:
            path = '/'.join(i.split('/')[:-1])
            if len(path) == 0:
                path = "."
            files += [ path+"/"+f for f in os.listdir(path) if fnmatch.fnmatch(f,i) ]

    # Handle "./" condition
    for count_i,i in enumerate(files):
        if i[:2] == "./":
            files[count_i] = files[count_i][2:]

    return list(set(files))

def GetInchi(elements, geometry):
    """
    # convert geo2inchi key

    # Open file for writing and write header
    # create a xyz file for openbabel input
    """
    write_xyz('inchi', elements, geometry)

    # create inchikey
    p = subprocess.Popen("obabel -ixyz inchi.xyz -oinchikey".split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
    inchikey, err = p.communicate()
    inchikey = str(inchikey, 'utf-8').split('\n')[0]  ##might be a problem if there are characters utf-8 can't handle
    # remove xyz file
    os.remove("inchi.xyz")

    return inchikey