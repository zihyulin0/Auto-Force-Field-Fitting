#!/bin/env python                                                                                                                                                              
import sys,os,argparse,subprocess,shutil,time,matplotlib,glob,getpass,json,fnmatch

# For plotting (Agg called needed for cluster image generation)
matplotlib.use('Agg') 
from pylab import *
import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit,minimize,lsq_linear
from copy import deepcopy
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/FF_functions')
from monitor_jobs import *
from parse_all import *
from file_parsers import *
import frag_gen_inter,frag_gen,xyz_to_orca,paramgen,merge_FF,extract_intramolecular_params,gen_md_for_sampling,gen_jobs_for_charges 
import gen_md_for_CHARMM,place_charge,generate_connolly
import restart_scans
import codecs,json
# Add TAFFY Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from transify import *
from adjacency import *
from file_parsers import *
from id_types import *
from kekule import *


def main(argv):
    


    user = getpass.getuser()
    run_dir = os.getcwd()

    parser = argparse.ArgumentParser(description='Driver script for taffi benchmark simulations. Expects a configuation file formatted with space delimited keyword arguments')

    #optional arguments                                                                                                                                                        
    parser.add_argument('-xyz', dest='xyz', default='fitcharge.xyz',
                        help = 'xyz file name, default: fitcharge.xyz')

    args=parser.parse_args()    
    
    element,geo = xyz_parse("fitcharge.xyz")
   
    with open('fitcharge.log','r') as f:
        charge_polar_dict = {}
        flag = 0
        for lc,lines in enumerate(f):
            fields = lines.split()
            if len(fields) < 4 and flag == 0: continue
          
            if flag == 0  and fields[0] == 'After' and fields[1] == 'polarizability' and fields[2] == 'scaling':
               flag = 1
               continue

            if flag == 0  and fields[0] == 'TOTAL' and fields[1] == 'POLARIZABILITY' :
               flag = 2
               continue

            if flag == 1 :
               if len(fields) == 0:
                  flag = 0
                  continue 
               atom_name = fields[0]
               charge_polar_dict[atom_name] = {'charge':None,'polar':'NA'}
               charge_polar_dict[atom_name]['charge'] = float(fields[1])
               if len(fields) == 3:
                  charge_polar_dict[atom_name]['polar'] = float(fields[2])
               continue

            if flag == 2 and len(fields) == 5:
               polar = float(fields[4])
               flag = 0 
               break

    print("Molecular polarizability: {:10.6f} (Angstrom^3)".format(polar))
               
    adj = Table_generator(element,geo)
    hybrid = Hybridization_finder(element,adj,force_opt=True)
    # one bond depth atom type (this makes making neighbor list very easy)
    atomT = id_types(element,adj,2)
    with open('fit_CHARMM.db','w') as f:

        f.write("\n# Charge definitions\n#\n{:10s} {:41s} {:20s}\n".format("#","Atom_type","Charge"))
        for i in set(atomT):
            idx = [ count_k for count_k,k in enumerate(atomT) if k == i ]
            atomnames = list(charge_polar_dict.keys()) 
            f.write("{:10s} {:40s} {:< 20.6f}\n".format("charge",str(i),charge_polar_dict[atomnames[idx[0]]]['charge']))
        f.write("\n# Atomic Polarizability\n#\n{:10s} {:41s} {:20s}\n".format("#","Atom_type","polarizability(A^3)"))
        for i in set(atomT):
            idx = [ count_k for count_k,k in enumerate(atomT) if k == i ]
            atomnames = list(charge_polar_dict.keys()) 
            if charge_polar_dict[atomnames[idx[0]]]['polar'] != 'NA':
               f.write("{:10s} {:40s} {:< 20.6f}\n".format("polar",str(i),charge_polar_dict[atomnames[idx[0]]]['polar']))
            else:
               f.write("{:10s} {:40s} {:< 20.6f}\n".format("polar",str(i),0.0))

    with open('drude.dff','w') as f:
         f.write("# units: kJ/mol, A, deg\n"+\
                 "# kforce is in the form k/2 r_D^2\n"+\
                 "# type  m_D/u   q_D/e    k_D   alpha/A3  thole\n")
         for key in charge_polar_dict:
            if charge_polar_dict[key]['polar'] != 'NA':
               f.write("{:5s} {:7.1f} {:7.1f} {:7.1f} {:10.6f} {:5.1f}\n".format(key,0.4,-1.0,4184.0,charge_polar_dict[key]['polar'],2.6))


            
         
         
    quit()
def delete_gbw():

    fileList = glob.glob('./*.gbw')

    for filePath in fileList:
       try:
          os.remove(filePath)
       except:
          print("Error while deleting file : ", filePath)

def delete_wildcard(pattern):
    fileList = glob.glob('./'+pattern+'.*')
      
    for filePath in fileList:
       try:
          os.remove(filePath)
       except:
          print("Error while deleting file : ", filePath)

def find_string(textfile,finding):
   # finding has to a byte string, ex:b'blabla'
   import mmap

   with open(textfile, 'rb', 0) as file, \
     mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ) as s:
      if s.find(finding) != -1:
         return True
      else: 
         return False 
      
   

def subproc_string(s):

    s = s.split()

    # Replace , arguments with spaces
    for count_i,i in enumerate(s):
        if "^" in i:
            s[count_i] = s[count_i].replace('^',' ')
    
    return s

def read_json(jsonfile):
    if os.path.isfile(jsonfile) is False:
      print("Error: json file: {} not found".format(jsonfile))
      quit()
    obj_text = codecs.open(jsonfile, 'r').read()
    jload = json.loads(obj_text)
    for key in jload:
      jload[key]['geo'] = np.array(jload[key]['geo'])
    
    return jload
def read_alljson(jsonfile):
    if os.path.isfile(jsonfile) is False:
      print("Error: json file: {} not found".format(jsonfile))
      quit()
    obj_text = codecs.open(jsonfile, 'r').read()
    jload = json.loads(obj_text)
    
    return jload
class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder, "a",buffering = 1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass

if __name__ == "__main__":
    main(sys.argv[1:])
