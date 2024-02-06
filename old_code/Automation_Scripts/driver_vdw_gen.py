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
import gen_jobs_for_vdw
from monitor_jobs import *
from parse_all import *
import frag_gen_inter,frag_gen,xyz_to_orca,paramgen,merge_FF,extract_intramolecular_params,gen_md_for_sampling,gen_jobs_for_charges 
import restart_scans,dihedral_restart
import codecs,json
from file_parsers import xyz_parse

# This script only did configuration generation, this is for Drude when the parametrized charge parameters already exist in the provided db,
# Drude only need the configurations

def main(argv):
    
    c = parse_configuration('config_negishi.txt')

    # Make config file absolute
    config = os.path.abspath('config_negishi.txt')

    # Get user/run_dir/start_log
    c["user"] = getpass.getuser()
    c["run_dir"] = os.getcwd()
    c["c_path"] = config

    # Initialize Params_for_Batch.db, the local set of
    # fixed parameters for use in the fitting procedure.
    if not os.path.isfile('Params_for_Batch.db'):
        shutil.copyfile(c['ff'],'Params_for_Batch.db')
    c['ff'] = "{}/Params_for_Batch.db".format(c['run_dir'])


    # Read all.json (contains model compound and dependency information for the batch 
    status = read_alljson('{}/all.json'.format(c["run_dir"]))

    """
    # Generate fake db for each generate so that gen_jobs_for_vdw knows what pairs of atomtypes to extract
    gens = sorted(status["gens"].keys())

    no_vdw = 'params_novdw.db'
    
    vdw_FF = get_all_vdw(c['ff'])

    for g in gens:
      mc_keys_previous = [ _ for _ in status["mc"].keys() if int(g)-1 in status["mc"][_]["gen"] ]
      f =  open('gen-{}.db'.format(g),'w') 
      with open(no_vdw,'r') as k:
         for lc,lines in enumerate(k):
            f.write(lines)
      for i in mc_keys_previous:
         for j in set(status['mc'][i]['atom_types']):
            f.write(vdw_FF[j])
      f.close()
    quit()
    """

    fit_charges(c,status)


    return

def get_all_vdw(FF):
   vdw_FF = {} # key: atype, value: that line
   with open(FF,'r') as f:
      for lc,lines in enumerate(f):
         fields = lines.split()
         if len(fields) > 1 and fields[0] == 'vdw':
            vdw_FF[fields[1]] = lines

   return vdw_FF

# Wrapper for fitting the charges for the generation. 
# This function performs (1) MD sampling, (2) QC job submission, (3) charge fitting to the output files. 
def fit_charges(c,status):


    # Find inchi keys for this generation
    #mc_keys = [ _ for _ in status["mc"].keys() if int(g) in status["mc"][_]["gen"] ]

    mc_keys = status["mc"].keys()


    #####################################################
    # Step 1 - Run md sampling within each inchi folder #
    #####################################################
    """
    jobids = []
    for inchi in mc_keys:

        # Get number of atoms in the molecule and number of molecules required for the simulation
        N_atoms = len(status["mc"][inchi]["atom_types"])
        N_mol = int(1000/N_atoms) + 1


        if not os.path.isdir('{}/vdw'.format(inchi)):
            os.mkdir('{}/vdw'.format(inchi))


        # Move into inchi folder and generate MD inputs within charges folder 
        os.chdir(inchi) # ***for gen_md.py has to go into the folder
        if not os.path.isdir('vdw/md-1'):
            #sublist = ('{}.xyz -N {} -T 298  -T_A 400 -t 20E5 -t_A 1E5 -sigma_scale 0.6 -o vdw/md-1  -gens {} --molecule -mixing_rule wh -FF '.format(inchi,N_mol,c['gens'])).split()
            sublist = ('{}.xyz -N {} -T 298  -T_A 400 -t 20E5 -t_A 1E5 -sigma_scale 1.0 -o vdw/md-1  -gens {} --molecule -mixing_rule wh -FF '.format(inchi,N_mol,c['gens'])).split()
            sublist.append('{} '.format(c['ff']))
            gen_md_for_sampling.main(sublist)

        # Move into charges folder and submit the MD job to the compute nodes (overwrites old job if equil.lammpstrj is missing)
        os.chdir('vdw/md-1')
        if not os.path.isfile("equil.lammpstrj"):
            substring = "python3 {}/Automation_Scripts/lammps_submit.py -f md-1.in.init -p {} -t {} -ppn {} -q {} -sched {} -size {} -o md -path_to_exe {} -npp {} --silent --overwrite"
            substring = substring.format(c["taffi_path"],c["charges_md_procs"],c["charges_md_wt"],c["charges_md_ppn"],c["charges_md_q"],c["charges_md_sched"],c["charges_md_size"],c["lammps_exe"],c['charges_md_npp'])
            output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
            output = str(output,'utf-8')
            jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
        os.chdir(c['run_dir'])

    # Wait until jobs complete 
    monitor_jobs(jobids,c["user"])
    print("COMPLETE: MD-sampling for  vdw  calculations")
    """


    ##########################################
    # Step 2 - Run QC Jobs to Obtain Charges #  
    ##########################################
    jobids = []
    gens = sorted(status["gens"].keys())

    for g in gens:
        mc_keys = [ _ for _ in status["mc"].keys() if int(g) in status["mc"][_]["gen"] ]

        for inchi in mc_keys:
           print(inchi)

           # Generate QC jobs for this model compound
           if not os.path.isfile("{}/vdw/md-1/equil.lammpstrj".format(inchi)):
               print("ERROR: MD-sampliing for vdw did not complete for {}".format(inchi))
               quit()
           
           # Generate QC jobs for this model compound
           if not os.path.isdir("{}/vdw/configs".format(inchi)):

               os.chdir('{}'.format(inchi))
               c['nod3'] = "--no_D3"
               sublist = "vdw -r_min_scale 1.5 -every 10 -N 2000 -p {} {} -QC_types dft -f {} --remove_cross -FF".format(c['vdw_qc_procs'],c['nod3'],c['functional']).split()
               sublist.append('{}/gen-{}.db'.format(c['run_dir'],g))
               gen_jobs_for_vdw.main(sublist)
           os.chdir(c['run_dir'])

    mc_keys = status["mc"].keys()
    for inchi in mc_keys:

       os.chdir('{}/vdw'.format(inchi))

       substring = "python3 {}/Automation_Scripts/orca_submit.py  -p {} -t {} -ppn {} -q {} -sched {} -size {} -path_to_exe {}  --silent"
       substring = substring.format(c["taffi_path"],c["vdw_qc_procs"],c["vdw_qc_wt"],c["vdw_qc_ppn"],c["vdw_qc_q"],c["vdw_qc_sched"],c["vdw_qc_size"],c["orca_exe"])
       output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
       output = str(output,'utf-8')
       jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
       os.chdir(c['run_dir'])

       # Wait until jobs complete
       monitor_jobs(jobids,c['user'])


       os.chdir(c['run_dir'])

    return

# Submits the vdw script for each inchi key in the generation
# the vdw MD, sampling, and fitting is particularly involved so it is outsourced to the dedicated vdw.py script
# XXX Lin has indicated there may be a subprocess related issue with these jobs. 7/20/20
def fit_vdw(g,c,status):

    # Find inchi keys for this generation
    mc_keys = [ _ for _ in status["mc"].keys() if int(g) in status["mc"][_]["gen"] ]

    process = {}
    jobids = []
    for inchi in mc_keys:

        # Skip complete jobs. Encountered because model compounds can be in multiple generations. 
        if(os.path.isfile(inchi+'/final_vdw.db')):
            continue

        # prepare subdirectory and intermediate_params.db with params from after charges.
        if not os.path.isdir('{}/vdw'.format(inchi)):
            os.mkdir('{}/vdw'.format(inchi))
        if not os.path.isfile('{}/vdw/intermediate_params.db'.format(inchi)):
            sublist = "{}/vdw/intermediate_params.db gen-{}_after_charges/params-DFT.db -only".format(inchi,g).split()
            sublist.append("atom bond angle torsion charge")
            merge_FF.main(sublist)

        # Launch the subprocess for parameterizing the vdw params
        os.chdir(inchi) ## vdw_loop has to be run in inchikey folder 
        command = ('{}/Automation_Scripts/vdw_loop.py -c {} -FF'.format(c["taffi_path"],c["c_path"])).split()
        command.append(c['ff'])
        process[inchi] = subprocess.Popen(command)
        os.chdir(c['run_dir']) # Back to run dir

    # Wait until all process complete
    running = {}
    for key in process:
        running[key] = process[key].poll() 
    while ( None in  list(running.values())):
        for key in process:
            running[key] = process[key].poll()
        time.sleep(20)

    # Check completion
    for inchi in mc_keys:
        if os.path.isdir(inchi+'/vdw') is False:
            print("An error occured duing the VDW fitting: folder i\"{}/vdw\" does not exist.".format(inchi))
            quit()
        elif os.path.isfile(inchi+'/vdw/status.txt') is False:
            print("An error occured duing the VDW fitting: file \"{}/vdw/status.txt\" does not exist, cannot determine convergence status.".format(inchi)) 
            quit()
        elif (find_string(inchi+'/vdw/status.txt',b'complete')) is False:
            print("An error occured duing the VDW fitting: \"{}/vdw\" is incomplete.".format(inchi))
            quit()

    # remove *.gbw files (take up a LOT of memory)
    delete_gbw()  

    # XXX: need a fix here to recognize the stopped vdw job when there's no atomtypes for that inchi key to parametrize
    # this normally happens when provided db file already has most of the parameters
    # Check that the final_vdw.db database was created
    for inchi in mc_keys:
        if(os.path.isfile(inchi+'/final_vdw.db')) is False:
            print("An error occured after the vdw parameterizations, {}/final_vdw.db was not found.".format(inchi))
            quit()

    # Update all.json
    status["gens"][g]["vdw"] = True
    json.dump(status,codecs.open(c["run_dir"]+'/all.json', 'w', encoding='utf-8'),indent=10)
    print("COMPLETE: VDW fitting for gen: {}".format(g))
    return

# XXX ADD FUNCTION DESCRIPTIONS. 7/14/20
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

# finding has to a byte string, ex:b'blabla'
def find_string(textfile,finding):
   import mmap
   with open(textfile, 'rb', 0) as file, mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ) as s:
       if s.find(finding) != -1:
           return True
       else: 
           return False 
      
# Replace ^ arguments with spaces
def subproc_string(s):
    s = s.split()
    for count_i,i in enumerate(s):
        if "^" in i:
            s[count_i] = s[count_i].replace('^',' ')
    
    return s

def read_json(jsonfile):
    if os.path.isfile(jsonfile) is False:
      print("Error: json file: {} not found".format(jsonfile))
      quit()
    obj_text = codecs.open(jsonfile, 'r').read()
    return json.loads(obj_text)

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
