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
import frag_gen_inter,frag_gen,xyz_to_orca,paramgen,merge_FF,extract_intramolecular_params,gen_md_for_sampling,gen_jobs_for_charges 
import restart_scans,dihedral_restart
import codecs,json
from file_parsers import xyz_parse

def main(argv):
    

    # Read all.json (contains model compound and dependency information for the batch 
    status = read_alljson('./all.json')
    g = 1
    mc_keys = [ _ for _ in status["mc"].keys() if int(g) in status["mc"][_]["gen"] ]
    print(" ".join(mc_keys))



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
