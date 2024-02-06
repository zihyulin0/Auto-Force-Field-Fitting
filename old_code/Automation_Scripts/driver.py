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
    
    parser = argparse.ArgumentParser(description='Driver script for taffi benchmark simulations. Expects a configuation file formatted with space delimited keyword arguments')

    #optional arguments                                                                                                                                                        
    parser.add_argument('-c', dest='config', default='config.txt',
                        help = 'The program expects a configuration file from which to assign various run conditions. (default: config.txt in the current working directory)')
    parser.add_argument('-o',dest='output',default='taffi.log',
                        help = 'log name for stdout')

    args=parser.parse_args()    
    
    # parse configuration dictionary (c)
    print("parsing configurations...")
    c = parse_configuration(args.config)

    # Make config file absolute
    args.config = os.path.abspath(args.config)

    # Get user/run_dir/start_log
    c["user"] = getpass.getuser()
    c["run_dir"] = os.getcwd()
    c["c_path"] =args.config
    sys.stdout = Logger(c["run_dir"]+'/'+args.output)

    # Initialize Params_for_Batch.db, the local set of
    # fixed parameters for use in the fitting procedure.
    if not os.path.isfile('Params_for_Batch.db'):
        shutil.copyfile(c['ff'],'Params_for_Batch.db')
    c['ff'] = "{}/Params_for_Batch.db".format(c['run_dir'])

    # Create all.json file and model compounds if not present
    if os.path.isfile(c["run_dir"]+'/all.json') is False:
        gen_models(c)

    # Read all.json (contains model compound and dependency information for the batch 
    status = read_alljson('{}/all.json'.format(c["run_dir"]))

    # Run geometry optimizations
    if not status["geoopt"]:
        geometry_optimizations(c,status)

    # Run mode scans
    if not status["modescan"]:
        mode_scans(c,status)

    # Main loop for fitting intermolecular params conforming to dependencies (i.e., by generation).
    gens = sorted(status["gens"].keys())
    print("Entering generation fitting loop...")
    for g in gens:
        
        # Run initial fits
        if not status["gens"][g]["initial_fit"]:
            fit_params(g,c,status)

        # Run charge fitting
        if not status["gens"][g]["charges"]:
            fit_charges(g,c,status)

        # Run after_charge fitting
        if not status["gens"][g]["after_charge"]:
            fit_params(g,c,status)

        # Run vdw fitting
        if not status["gens"][g]["vdw"]:
            fit_vdw(g,c,status)

        # Run final param fitting
        if not status["gens"][g]["final_fit"]:
            fit_params(g,c,status)

    return

# Submits and monitors the frag_gen_all.py job
def gen_models(c):
    jobids = []
    sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
    sublist.append('python3 {}/FF_functions/frag_gen_all.py \'*.xyz\' -FF \'{}\' -gens {}'.format(c["taffi_path"],c["ff"],c['gens']))
    sublist= sublist + ('frag_gen_all -p 1 -t {} -q {} -ppn {}'.format(c["param_fit_wt"],c["param_fit_q"],c["param_fit_ppn"])).split()
    output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    output = str(output,'utf-8')
    jobids += [ output.split("\n")[-2].split()[-1]]
    monitor_jobs(jobids,c["user"])
    return

# Submits and monitors the geometry_optimizations jobs
# (1) Coverts the high-level xyz files generated by frag_gen_all.py to orca geometry optimization jobs, (2) handles resubmission conditions. 
def geometry_optimizations(c,status):        

    # Loop over inchi keys and generate run folders and input files
    print( "\trunning geometry optimizations...")    
    for key in status["mc"].keys():
        other_opt = {}
        other_opt['no_freq'] = 1
        other_opt['proc_num'] = c['param_geoopt_procs']
        other_opt['random_factor'] = 0.0
        other_opt['charge'] = c['charge']
        other_opt['functional'] = c['functional']
        other_opt['basis'] = c['basis']
        if c['functional'] == 'wB97X-D3':
            other_opt['D3_option'] = ""
        e,g = xyz_parse('{}/{}.xyz'.format(key,key))
        xyz_to_orca.fun(c,key,{"elements":e,"geo":g},other_opt)

    # Submit the geometry optimizations
    jobids = []
    substring = "{}/Automation_Scripts/orca_submit.py -p {} -d geoopt -t {} -ppn {} -q {} -sched {} -size {} -o optimizations -path_to_exe {}  --silent"
    substring = substring.format(c["taffi_path"],c["param_geoopt_procs"],c["param_geoopt_wt"],c["param_geoopt_ppn"],c["param_geoopt_q"],c["param_geoopt_sched"],c["param_geoopt_size"],c["orca_exe"])
    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    output = str(output,'utf-8')
    jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
    monitor_jobs(jobids,c["user"])

    # Remove job submission files
    delete_wildcard('optimization')

    # Remove *.gbw files
    delete_gbw()

    # Check if geoopt jobs are incomplete/need resubmission
    count = 0
    resub_flag = 1
    while resub_flag == 1:
        resub_flag = 0
        if count >= 5:
            print("ERROR: five restarts have been attempted and the geoopt jobs still haven't converged. Something bad is happening.")
            quit()

        # Check status of jobs
        for key in status["mc"].keys():
            geofolders = [ key+'/Intra/'+i for i in os.listdir(key+'/Intra') if os.path.isdir(key+'/Intra/'+i) == True and fnmatch.fnmatch(i,"*geoopt*")] 
            geoopt_output = [ os.path.isfile(i+'/geoopt.out') for i in geofolders ]
            if geofolders == []:
                print("ERROR: {}/Intra/*_geoopt* folder wasn't found".format(key))
                quit()
            if False in geoopt_output:
                print("WARNING: {}/Intra/*_geoopt*/geoopt.out wasn't found".format(key))
                resub_flag = 1
            flag = [ find_string(i+'/geoopt.out',b'THE OPTIMIZATION HAS CONVERGED') for i in geofolders]
            if False in flag:
                print("WARNING: {}/Intra/*_geoopt* didn't run to completion...".format(key))
                resub_flag = 1

        # Attempt resubmission if the resub_flag has been triggered
        if(resub_flag):
            count += 1
            print("resub of geoopt has been triggered")
            jobids = []
            substring = "{}/Automation_Scripts/orca_submit.py -d geoopt -p {} -t {} -ppn {} -q {} -sched {} -size {} -o optimizations -path_to_exe {}  --silent --resubmit"
            substring = substring.format(c["taffi_path"],c["param_geoopt_procs"],c["param_geoopt_wt"],c["param_geoopt_ppn"],c["param_geoopt_q"],c["param_geoopt_sched"],c["param_geoopt_size"],c["orca_exe"])
            output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
            output = str(output,'utf-8')
            jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]

            # Wait until jobs complete
            monitor_jobs(jobids,c["user"])
    
    # Toggle flag in status
    status["geoopt"] = True
    json.dump(status,codecs.open(c["run_dir"]+'/all.json', 'w', encoding='utf-8'),indent=10)
    return

# Wrapper for performing the quantum chemistry scans for the intramolecular modes. 
# This function handles generation of the mode scans, QC job submission, completion checks, restart and resubmission conditions associated with the intramolecular modes. 
def mode_scans(c,status):

    # Loop over inchi keys and generate run folders and input files
    print( "\trunning mode scans...")    
    for key in status["mc"].keys():

        # Generate bond/angle mode scan job if needed
        if not os.path.isdir(key+'/Intra/bonds_angles'):            
            command_list = ('-p {} -q {} -gens {} -theory dft -f {} -b {} {} -inchi {} -modes'.format(c['param_ba_procs'],status["mc"][key]["q"],c['gens'],c['functional'],c['basis'],c['nod3'],key)).split()
            command_list.append("bonds angles")
            paramgen.main(command_list)

        # Generate dihedral mode scan job if needed
        if not os.path.isdir(key+'/Intra/dihedrals'):
            command_list = ('-p {} -q {} -gens {} -theory dft -f {} -b {} -d_step 10.0 {} -inchi {} --scan -modes'.format(c['param_d_procs'],status["mc"][key]["q"],c['gens'],c['functional'],c['basis'],c['nod3'],key)).split() 
            command_list.append("dihedrals")
            paramgen.main(command_list)

    # Submit scans
    print("submitting mode scans...")
    substring = "python3 {}/Automation_Scripts/orca_submit.py -f *in -d bonds_angles -p {} -t {}  -q {} -sched {} -ppn {} -size {} -path_to_exe {} --silent -o bonds_angles"
    substring = substring.format(c["taffi_path"],c["param_ba_procs"],c["param_ba_wt"],c["param_ba_q"],c["param_ba_sched"],c["param_ba_ppn"],c["param_ba_size"],c["orca_exe"])
    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    output = str(output,'utf-8')
    jobids = [ m.split()[-1] for m in output.split("\n")[:-1]]

    substring = "python3 {}/Automation_Scripts/orca_submit.py -f *in -d dihedrals -p {} -t {}  -q {} -sched {} -ppn {} -size {} -path_to_exe {} --silent -o dihedrals"
    substring = substring.format(c["taffi_path"],c["param_d_procs"],c["param_d_wt"],c["param_d_q"],c["param_d_sched"],c["param_d_ppn"],c["param_d_size"],c["orca_exe"])
    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    output = str(output,'utf-8')
    jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]

    # Wait until jobs complete
    monitor_jobs(jobids,c["user"])

    # Remove job submission files
    delete_wildcard('bonds_angles')
    delete_wildcard('dihedrals')

    # Resubmission loop
    count = 0
    resub_flag = 1
    while resub_flag == 1:
        resub_flag = 0

        # Break if too many attempts
        if count >= 5:
            print("\tERROR: five restarts have been attempted and the mode scans still haven't converged. Something bad is happening.")
            quit()

        # Remove old test, REDO, and RESTART folders (sometimes the result of an interuption mid-cycle)
        # NOTE: REDO is used for bonds/angle restarts and RESTART is used for incomplete dihedral scans. Maybe there is a better system where all of these have the same name but are contextually distinguished. 
        for i in [ "{}/{}".format(dp,d0) for dp,d,f in os.walk('.') for d0 in d if fnmatch.fnmatch(d0,"test") or fnmatch.fnmatch(d0,"RESTART") or fnmatch.fnmatch(d0,"REDO") ]:
            shutil.rmtree(i)

        # check bond angle submissions with extract_intramolecular_paraams.py
        substring = "python {}/FF_functions/extract_intramolecular_params.py -gens {} --find_restarts -lammps_exe {} -FF {} -o test -modes bonds^angles^harm_dihedrals".format(c["taffi_path"],c['gens'],c['lammps_exe'],c['ff'])
        output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()
        flag = find_string("test/extract_intramolecular.log",b"Generating the input files for a reoptimization based on the lowest energy")   #search string must be byte object, hence the "b"
        shutil.rmtree('test')

        # Handle flexible dihedral resubmission
        substring = "python {}/FF_functions/dihedral_restart.py . --verbose".format(c["taffi_path"])
        output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0]
        if (output.split("\n")[0]) != '':
            D_flag = len(output.split("\n"))
        else:
            D_flag = 0
                     
        # Handles dihedral resubmission for incomplete scans
        if (D_flag != 0):
            print("\tSome flexible dihedral scans did not complete. Resuming (attempt {})...".format(count))

            # Submit restarted dihedral scans
            D_ids_sub = []
            substring = "python3 {}/Automation_Scripts/orca_submit.py -f *in -d RESTART -p {} -t {}  -q {} -sched {} -ppn {} -size {} -path_to_exe {} --silent -o dihedrals_restart"
            substring = substring.format(c["taffi_path"],c["param_d_procs"],c["param_d_wt"],c["param_d_q"],c["param_d_sched"],c["param_d_ppn"],c["param_d_size"],c["orca_exe"])
            print(substring)
            output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf-8').communicate()[0]
            D_ids_sub += [ m.split()[-1] for m in output.split("\n")[:-1]]
               
        # Handles bond/angle resubmission/reoptimization @ better minima
        if flag is True:
            print("Some bond/angle jobs failed to converge. Resubmitting (attempt {})...".format(count))
            jobids = []
            substring = "python3 {}/Automation_Scripts/orca_submit.py -f *in -d REDO -p {} -t {}  -q {} -sched {} -ppn {} -size {} -o ba_geo_resub -path_to_exe {} --silent"
            substring = substring.format(c["taffi_path"],c["param_geoopt_procs"],c["param_geoopt_wt"],c["param_geoopt_q"],c["param_geoopt_sched"],c["param_geoopt_ppn"],c["param_geoopt_size"],c["orca_exe"])
            output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf-8').communicate()[0]
            jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
            monitor_jobs(jobids,c["user"])

            # Generate and run new mode scans
            restart_scans.main(("Intra -f {} -b {} -gens {} {}".format(c['basis'],c['functional'],c['gens'],c['nod3'])).split())
            substring = "python3 {}/Automation_Scripts/orca_submit.py -f *in  -p {} -t {}  -q {} -sched {} -ppn {} -size {} -o ba_scan_resub -path_to_exe {} --silent"
            substring = substring.format(c["taffi_path"],c["param_geoopt_procs"],c["param_geoopt_wt"],c["param_geoopt_q"],c["param_geoopt_sched"],c["param_geoopt_ppn"],c["param_geoopt_size"],c["orca_exe"])
            output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf-8').communicate()[0]
            jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
            monitor_jobs(jobids,c["user"])

        # Run the dihedral post processing if jobs were submitted.
        if (D_flag != 0):  

           # Wait until jobs complete
           monitor_jobs(D_ids_sub,c["user"])

           # Run post processing
           dihedral_restart.main(". --post".split()) 

        if (D_flag == 0 ) and flag is False:
           break


    # Remove *.gbw files 
    delete_gbw() # XXX Update to act recursively. 7/14/20

    # Toggle flag in status
    status["modescan"] = True
    json.dump(status,codecs.open(c["run_dir"]+'/all.json', 'w', encoding='utf-8'),indent=10)
    return

# Wrapper for fitting the intramolecular parameters for the generation.
# This function handles the three cases for which parameterizations occur each generation (initial_fit, after_charge, final_fit)
def fit_params(g,c,status):

    # Figure out whether initial_fit, after_charge, or final_fit stage needs to be performed
    stage = next( _ for _ in ["initial_fit","after_charge","final_fit"] if not status["gens"][g][_] )

    # Find inchi keys for this generation
    mc_keys = [ _ for _ in status["mc"].keys() if int(g) in status["mc"][_]["gen"] ]

    # case: initial fit
    # No charges or VDW parameters for the gen are available, so the fit is performed based on charges from the equilibrium geometry and UFF params
    if stage == "initial_fit":

        sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
        sublist.append('{}python {}/FF_functions/extract_intramolecular_params.py -f \'{}\' -o gen-{}_initial_params -FF \'{}\' -min_cycles 5 -gens {} -max_cycles 10 -N_sc 1 -lammps_exe {} --mixing_rule lb -batch {}'.format(\
        c["module_string"],c["taffi_path"]," ".join(mc_keys),g,c["ff"],c['gens'],c["lammps_exe"],g))
        sublist= sublist + ('initial_params  -p 1 -t {} -q {} -ppn {}'.format(c["param_fit_wt"],c["param_fit_q"],c["param_fit_ppn"])).split()
        output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
        output = str(output,'utf-8')
        jobids = [ output.split("\n")[-2].split()[-1]]

        # Wait until jobs complete
        monitor_jobs(jobids,c['user'])
        # Check completion and update status/all.json
        if os.path.isdir("gen-{}_initial_params".format(g)) is False or os.path.isfile("gen-{}_initial_params/params-DFT.db".format(g)) is False:
            print("An error occured during the gen {} initial parameter fit.".format(g))
            quit()
        else:
            status["gens"][g]["initial_fit"] = True
            json.dump(status,codecs.open(c["run_dir"]+'/all.json', 'w', encoding='utf-8'),indent=10)

    # case: after_charge
    # Charge force field files are read. No VDW parameters for the gen are available, so the fit is performed based on the UFF params
    if stage == "after_charge":

        sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
        sublist.append('{}python {}/FF_functions/extract_intramolecular_params.py -f \'{}\' -o gen-{}_after_charges'.format(c["module_string"],c["taffi_path"]," ".join(mc_keys),g)+\
                       ' -FF \'{} {}\''.format(c["ff"]," ".join(["{}/charges/CHELPG_calcs/charges/fit_charges.db".format(_) for _ in mc_keys ]))+\
                       ' -min_cycles 5 -max_cycles 10 -N_sc 1 -gens {} -lammps_exe {} --mixing_rule lb -batch {}'.format(c['gens'],c["lammps_exe"],g))
        sublist= sublist + ('after_charges  -p 1 -t {} -q {} -ppn {}'.format(c["param_fit_wt"],c["param_fit_q"],c["param_fit_ppn"])).split()
        output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
        output = str(output,'utf-8')
        jobids = [ output.split("\n")[-2].split()[-1]]

        # Wait until jobs complete
        monitor_jobs(jobids,c['user'])

        # Check completion and update status/all.json        
        if os.path.isdir("gen-{}_after_charges".format(g)) is False or os.path.isfile("gen-{}_after_charges/params-DFT.db".format(g)) is False:
            print("An error occured during the gen {} after charges parameter fit.".format(g))
            quit()
        else:
            status["gens"][g]["after_charge"] = True
            json.dump(status,codecs.open(c["run_dir"]+'/all.json', 'w', encoding='utf-8'),indent=10)

    # case: final_fit
    # Charge and VDW force field files are read. Everything is available to perform the fit, mixing rule is changed. 
    if stage == "final_fit":

        sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
        sublist.append('{}python {}/FF_functions/extract_intramolecular_params.py -f \'{}\' -o gen-{}_final_params'.format(c["module_string"],c["taffi_path"]," ".join(mc_keys),g)+\
                       ' -FF \'{} {}\''.format(c["ff"]," ".join(["{}/charges/CHELPG_calcs/charges/fit_charges.db {}/final_vdw.db".format(_,_) for _ in mc_keys ]))+\
                       ' -min_cycles 5 -max_cycles 10 -N_sc 3 -lammps_exe {} -gens {} --mixing_rule wh -batch {}'.format(c["lammps_exe"],c['gens'],g))
        sublist= sublist + ('final_params  -p 1 -t {} -q {} -ppn {}'.format(c["param_fit_wt"],c["param_fit_q"],c["param_fit_ppn"])).split()
        output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
        output = str(output,'utf-8')
        jobids = [ output.split("\n")[-2].split()[-1]]

        # Wait until jobs complete
        monitor_jobs(jobids,c['user'])

        # Check completion and update status/all.json        
        if os.path.isdir("gen-{}_final_params".format(g)) is False or os.path.isfile("gen-{}_final_params/params-DFT.db".format(g)) is False:
            print("An error occured during the gen {} final parameter fit.".format(g))
            quit()
        else:
            status["gens"][g]["final_fit"] = True
            json.dump(status,codecs.open(c["run_dir"]+'/all.json', 'w', encoding='utf-8'),indent=10)

        # Merge final params for this generation with Params_for_Batch.db
        option = {}
        option['master'] = '{}/Params_for_Batch.db'.format(c['run_dir'])
        option['new_params'] = '{}/gen-{}_final_params/params-DFT.db'.format(c['run_dir'],g)
        option['no_ct'] = True
        merge_FF.fun(option)

    return

# Wrapper for fitting the charges for the generation. 
# This function performs (1) MD sampling, (2) QC job submission, (3) charge fitting to the output files. 
def fit_charges(g,c,status):

    # Absolute paths of FF are needed for the subdirectory calls
    print("fitting gen-{} charges".format(g))
    init_FF = os.path.abspath("gen-{}_initial_params/params-DFT.db".format(g))

    # Find inchi keys for this generation
    mc_keys = [ _ for _ in status["mc"].keys() if int(g) in status["mc"][_]["gen"] ]

    #####################################################
    # Step 1 - Run md sampling within each inchi folder #
    #####################################################
    jobids = []
    for inchi in mc_keys:

        # Get number of atoms in the molecule and number of molecules required for the simulation
        N_atoms = len(status["mc"][inchi]["atom_types"])
        N_mol = int(1000/N_atoms) + 1

        # Move into inchi folder and generate MD inputs within charges folder 
        os.chdir(inchi) # ***for gen_md.py has to go into the folder
        if not os.path.isdir('charges'):
            sublist = ('{}.xyz -N {} -T 298  -T_A 400 -t 1E5 -t_A 1E5 -charge_scale 1.0 -o charges -q {} -gens {} --molecule -FF '.format(inchi,N_mol,status["mc"][inchi]["q"],c['gens'])).split()
            sublist.append('{} {}'.format(c['ff'],init_FF))
            gen_md_for_sampling.main(sublist)

        # Move into charges folder and submit the MD job to the compute nodes (overwrites old job if equil.lammpstrj is missing)
        os.chdir('charges')
        if not os.path.isfile("equil.lammpstrj"):
            substring = "python3 {}/Automation_Scripts/lammps_submit.py -f charges.in.init -p {} -t {} -ppn {} -q {} -sched {} -size {} -o charges -path_to_exe {} -npp {} --silent --overwrite"
            substring = substring.format(c["taffi_path"],c["charges_md_procs"],c["charges_md_wt"],c["charges_md_ppn"],c["charges_md_q"],c["charges_md_sched"],c["charges_md_size"],c["lammps_exe"],c['charges_md_npp'])
            output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
            output = str(output,'utf-8')
            jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
        os.chdir(c['run_dir'])

    # Wait until jobs complete 
    monitor_jobs(jobids,c["user"])
    print("COMPLETE: MD-sampling for gen-{} partial-charge calculations".format(g))

    ##########################################
    # Step 2 - Run QC Jobs to Obtain Charges #  
    ##########################################
    jobids = []
    for inchi in mc_keys:

        # Generate QC jobs for this model compound
        if not os.path.isfile("{}/charges/equil.lammpstrj".format(inchi)):
            print("ERROR: MD-sampliing for partial charges did not complete for {}".format(inchi))
            quit()
        
        # Generate QC jobs for this model compound
        if not os.path.isdir("{}/charges/CHELPG_calcs".format(inchi)):
            os.chdir("{}/charges".format(inchi))
            substring = 'equil.lammpstrj charges.map -every 10 -N 200 -mol_list 0 -p 8 -QC_types dft -o CHELPG_calcs -q {} -f {} {}'
            substring = substring.format(status["mc"][inchi]["q"],c['functional'],c['nod3'])
            gen_jobs_for_charges.main(substring.split())

            os.chdir('CHELPG_calcs')
            substring = "python3 {}/Automation_Scripts/orca_submit.py -p {} -t {} -ppn {} -q {} -sched {} -size {} -o charges -path_to_exe {}  --silent"
            substring = substring.format(c["taffi_path"],c["charges_qc_procs"],c["charges_qc_wt"],c["charges_qc_ppn"],c["charges_qc_q"],c["charges_qc_sched"],c["charges_qc_size"],c["orca_exe"])
            output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
            output = str(output,'utf-8')
            jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]

        os.chdir(c['run_dir'])

    # Wait until jobs complete 
    monitor_jobs(jobids,c["user"])
    print("COMPLETE: QC calculations for gen-{} partial-charges".format(g))

    ####################################
    # Step 3  - Perform charge fitting #
    ####################################      
    jobids = []
    for inchi in mc_keys:

        if os.path.isdir("{}/charges/CHELPG_calcs".format(inchi)) is False or os.path.isfile("{}/charges/CHELPG_calcs/configs/0/0_charges.vpot".format(inchi)) is False:
            print("ERROR: QC for partial chages did not complete for {}".format(inchi))
            quit()
            
        # If the fit_charges.db file does not exist, then remove the fit directory and attempt to refit the charges 
        os.chdir(inchi+'/charges') 
        if(os.path.isdir('CHELPG_calcs/charges')) and ((os.path.isfile('CHELPG_calcs/charges/fit_charges.db')) is False) :
            shutil.rmtree('CHELPG_calcs/charges')

        # Fit the charges
        if not os.path.isdir('CHELPG_calcs/charges'):

            # Have to use sublist or the split method would split the python function call to two
            sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
            sublist.append('{}/FF_functions/extract_charges.py CHELPG_calcs -FF \'{}\' -o charges -w_hyper 0.0 -w_dipole 0.1 --two_step -w_qtot 1.0 -q {} -gens {} --keep_min'.format(c["taffi_path"],c["ff"],c["charge"],c["gens"]))
            sublist= sublist + ('charge_parse  -p 1 -t {} -q {} -ppn {}'.format(c["param_fit_wt"],c["param_fit_q"],c["param_fit_ppn"])).split()
            output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
            output = str(output,'utf-8')
            jobids += [ output.split("\n")[-2].split()[-1]]
        os.chdir(c['run_dir'])

    # Wait until jobs complete 
    monitor_jobs(jobids,c["user"])

    # Check completion
    for inchi in mc_keys:
        if os.path.isfile("{}/charges/CHELPG_calcs/charges/fit_charges.db".format(inchi)) is False: 
            print("An error occured during the gen {} after charge fitting of {}.".format(g,inchi))
            quit()

    # Update all.json
    status["gens"][g]["charges"] = True
    json.dump(status,codecs.open(c["run_dir"]+'/all.json', 'w', encoding='utf-8'),indent=10)
    print("COMPLETE: Fitting partial charges for gen-{} partial-charges".format(g))
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
        #command = ('{}/Automation_Scripts/vdw_loop.py -c {} -FF'.format(c["taffi_path"],c["c_path"])).split()
        command = ('{}/Automation_Scripts/vdw_loop_negishi.py -c {} -FF'.format(c["taffi_path"],c["c_path"])).split()
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

# Description: finds the number of disconnected subnetworks in the 
#              adjacency matrix, which corresponds to the number of 
#              separate molecules.
#
# Inputs:      adj_mat: numpy array holding a 1 in the indices of bonded
#                        atom types. 
#
# Returns:     mol_count: scalar, the number of molecules in the adj_mat
def mol_count(adj_mat):
    
    # Initialize list of atoms assigned to molecules and a counter for molecules
    placed_idx = []    
    mol_count = 0

    # Continue the search until all the atoms have been assigned to molecules
    while len(placed_idx)<len(adj_mat):

        # Use sequential elements of the adj_mat as seeds for the spanning network search
        for count_i,i in enumerate(adj_mat):

            # Only proceed with search if the current atom hasn't been placed in a molecule
            if count_i not in placed_idx:

                # Increment mol_count for every new seed and add the seed to the list of placed atoms
                mol_count += 1               
                placed_idx += [count_i]
                
                # Find connections
                idx = [ count_j for count_j,j in enumerate(i) if j==1 and count_j not in placed_idx ]
                
                # Continue until no new atoms are found
                while len(idx) > 0:
                    current = idx.pop(0)
                    if current not in placed_idx:
                        placed_idx += [current]
                        idx += [ count_k for count_k,k in enumerate(adj_mat[current]) if k == 1 and count_k not in placed_idx ]
    return mol_count


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
