#!/bin/env python                                                                                                                                                              
import sys,os,argparse,subprocess,shutil,time,matplotlib,glob,fnmatch,getpass
import numpy as np
from matplotlib import pyplot as plt
# For plotting (Agg called needed for cluster image generation)
matplotlib.use('Agg') 
from pylab import *
from scipy.stats import linregress
from scipy.optimize import curve_fit,minimize,lsq_linear
from copy import deepcopy
import random
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/FF_functions')
import db_to_charmm,xyz_to_pdb,density_converter,add_time_pdb,read_charmm_energy,gen_itp_newdrude,crd_to_xyz
import add_drude_pdb
from monitor_jobs import *

def main(argv):


    parser = argparse.ArgumentParser(description='Driver script for taffi benchmark simulations. Expects a configuation file formatted with space delimited keyword arguments')

    #optional arguments                                                                                                                                                        
    parser.add_argument('-c', dest='config', default='config.txt',
                        help = 'The program expects a configuration file from which to assign various run conditions. (default: config.txt in the current working directory)')
    parser.add_argument('-nonpolar_NMA', dest='nonpolar_NMA', default='',
                        help = 'polar normal mode analysis need nonpolar NMA to identify obolete modes, for polar NMA')
    parser.add_argument('--polar', dest='polar_flag', default=False, action='store_const', const=True,
                        help = 'When present, will generate Drude FF. (default: False)')
    parser.add_argument('--shake', dest='shake_flag', default=False, action='store_const', const=True,
                        help = 'When present, will use shake to contrain h bonds, for MD (default: False)')
    parser.add_argument('--create', dest='create_flag', default=False, action='store_const', const=True,
                        help = 'When present, will only generate run file but not analysis. for MD (default: False)')

    sys.stdout = Logger('FOM_driver')
    print("PROGRAM CALL: python {}\n".format(' '.join([i for i in sys.argv])))

    # parse configuration dictionary (c)
    print("parsing benchmark configurations...")
    global global_c
    c = parse_configuration(parser.parse_args())
    global_c  = c

    # Make parse inputs
    args=parser.parse_args(argv)

    # run geometry optimizations and frequency analysis
    print("\nentering normal mode phase...")
    run_normal_analysis(c,args.nonpolar_NMA,args.polar_flag)

    # run md
    print("\nentering molecular dynamics phase...")
    #run_md(c,args.polar_flag,args.shake_flag,args.create_flag)

    quit()

# Wrapper function for the normal mode analysis
# This use SHAKE to contraint bonds no matter what
def run_normal_analysis(c,nonpolar_NMA,polar_flag):
    
    # Make the Normal_Modes 
    if os.path.isdir("Normal_Modes") is False:
        os.mkdir("Normal_Modes")

    # Loop over the xyz files, copy them, create canonical conformer, and generate orca input files
    print("\tgenerating input files for  geometry optimizations...")
    for i in c["xyz"]:
        name=".".join(i.split('.')[:-1])
        if os.path.isdir("Normal_Modes/{}_freq".format(name)) is False:
            
            # transify the geometry
            subprocess.check_output(subproc_string("python {}/Lib/transify.py {} {}.opt.xyz".format(c["taffi_path"],i,name)), encoding='utf8')


            # obminimize output pdf instead of xyz and transify already done obminimize within. Lin 10/26/22
            """
            # use openbabel to minimize the geometry via the UFF force-field
            output = subprocess.Popen(subproc_string("obminimize -ff uff -sd -c 1e-20 -n 100000 straightened.xyz  > {}.opt.xyz 2> /dev/null".format(name)), stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0]
            with open("{}.opt.xyz".format(name),'w') as f:
                for lines in output:
                    fields = lines.split()
                    if len(fields) > 0 and fields[0] == "WARNING:":
                        continue
                    else:
                        f.write(lines)
            """
                    
            if os.path.isfile("{}.opt.xyz".format(name)) is False:
                print("ERROR: file {} failed at openbabel optimization...".format(i))
                quit()

            output = subprocess.Popen(subproc_string("python {}/FF_functions/xyz_to_orca_FOM.py {}.opt.xyz -f {} -p {} -o {}_freq.in {}".format(c["taffi_path"],name,c["functional"],c["fom_geoopt_procs"],name,c["d3_option"])),stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0]
            os.mkdir("Normal_Modes/{}_freq".format(name))
            shutil.move("{}.opt.xyz".format(name),"Normal_Modes/{}_freq/{}.opt.xyz".format(name,name))
            shutil.move("{}_freq.in".format(name),"Normal_Modes/{}_freq/.".format(name))
            #os.remove("straightened.xyz")

    # Submit the geometry optimizations
    jobids = []
    print("\trunning geometry optimizations...")
    #substring = "python {}/Automation_Scripts/orca_submit.py -f *freq.in -p {} -t {} -ppn {} -q {} -sched {} -o frequencies -path_to_exe {} -size {} --silent"
    substring = "python {}/Automation_Scripts/orca_submit.py -d Normal_Modes/ -p {} -t {} -ppn {} -q {} -sched {} -o frequencies -path_to_exe {} -size {} --silent"
    substring = substring.format(c["taffi_path"],c["fom_geoopt_procs"],c["fom_geoopt_wt"],c["fom_geoopt_ppn"],c["fom_geoopt_q"],c["fom_geoopt_sched"],c["orca_exe"],c["fom_geoopt_size"])
    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0]
    jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
    monitor_jobs(jobids,'lin1209')

    # Remove submission files
    for f in glob.glob("frequencies.*"):
        os.remove(f)

    # Calculate the classical modes
    jobids=[]

    f = open('Normal_Modes/summary.txt','w')
    f.write("{:<20s} {:<20s} {:<20s}\n".format('name','dipole(Debye)','polarizability(A^3)'))
    
    for count_i,i in enumerate(c["xyz"]):
        name=".".join(i.split('.')[:-1])
        name_sub=name.replace('(','-')
        name_sub=name_sub.replace(')','-')
        name_safe=esc_str(name)

        # Check for incomplete
        if os.path.isdir("Normal_Modes/{}_minimized".format(name)) is True and os.path.isfile("Normal_Modes/{}_minimized/minimized.xyz".format(name)) is False:
            shutil.rmtree("Normal_Modes/{}_minimized".format(name))

        # Run the minimization and mode calculation
        if os.path.isdir("Normal_Modes/{}_minimized".format(name)) is False:

            # Check the requisite data is in place. 
            if os.path.isdir("Normal_Modes/{}_freq".format(name)) is False or os.path.isfile("Normal_Modes/{}_freq/{}_freq.out".format(name,name)) is False:
                print("An error occured when trying to process {}_freq. Exiting...".format(name))
                quit()

            current_dir = os.getcwd()
            xyzfile = '{}/Normal_Modes/{}_freq/geo_opt.xyz'.format(current_dir,name)
            min_folder = '{}/Normal_Modes/{}_minimized'.format(current_dir,name)
            qc_folder = '{}/Normal_Modes/{}_freq'.format(current_dir,name)

            # Assemble the submission string for the batch submission file
            if polar_flag:
               substring ="python {}/FF_functions/charmm_minimize_molecule.py {} -o {} --polar -gens {} -FF {} > /dev/null \n"
            else:
               substring ="python {}/FF_functions/charmm_minimize_molecule.py {} -o {} -gens {} -FF {} > /dev/null \n"
            substring=substring.format(c["taffi_path"],xyzfile,min_folder,global_c['gens'],global_c['fom_ff'])                

            substring+="if [ ! -f {}_minimized/minimized.xyz ]; then\n"
            substring+="    echo 'An error occured during the minimization of {}.xyz'\n"
            substring+="    exit\n"
            substring+="fi\n\n"
            substring = substring.format(name_safe,name)
            substring+="# Calculate the FF normal modes, and create an aligned geometry\n"
            if polar_flag:
               if nonpolar_NMA == '':
                  print('ERROR: polar NMA needs nonpolar result to identify obsolte modes to compare with QC')
                  quit()
               substring+="python {}/Parsers/parse_normal_modes_charmm.py -main_out {}/min.log -QC_out {}/{}_freq.out -nonpolar_out {}/{}_minimized/min.log -o {}_FF   > /dev/null \n"
               substring = substring.format(c["taffi_path"],min_folder,qc_folder,name_safe,nonpolar_NMA,name_safe,name_safe)
            else:
               substring+="python {}/Parsers/parse_normal_modes_charmm.py -main_out {}/min.log -QC_out {}/{}_freq.out -nonpolar_out {}/min.log -o {}_FF   > /dev/null \n"
               substring = substring.format(c["taffi_path"],min_folder,qc_folder,name_safe,min_folder,name_safe)
            substring+="python {}/Parsers/align_points.py {}_freq/geo_opt.xyz {}_minimized/minimized.xyz -o {}_aligned > /dev/null  \n"
            substring = substring.format(c["taffi_path"],name_safe,name_safe,name_safe)

            # Write the batch script and submit the job
            shell_submit(substring,"{}.submit".format(count_i),1,1,c["fom_geoopt_q"],c["fom_geoopt_ppn"],batch=c["batch"],run_dir="{}/{}".format(os.getcwd(),"Normal_Modes"))
            output = subprocess.Popen([c["sub_cmd"],"{}.submit".format(count_i)],stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0].strip("\r\n")
            jobids += [ output.split()[-1] ]
         
         

    # Wait until the jobs complete
    monitor_jobs(jobids,'lin1209')
    for i in c["xyz"]:
        name=".".join(i.split('.')[:-1])
        name_sub=name.replace('(','-')
        name_sub=name_sub.replace(')','-')
        name_safe=esc_str(name)

        mu = get_mu("Normal_Modes/{}_minimized/mu.log".format(name)) 
        polar = get_polar("Normal_Modes/{}_minimized/min.log".format(name))
        f.write("{:<20s} {:<20.6f} {:<20.6f}\n".format(name,mu,polar))

    """
    # Remove the submission scripts
    for i in c["xyz"]:
        name=".".join(i.split('.')[:-1])
        name_sub=name.replace('(','-')
        name_sub=name_sub.replace(')','-')
        if os.path.isfile("{}.submit".format(name_sub)) is True:
            for f in glob.glob("{}.submit*".format(name_sub)):
                os.remove(f)
    """
    f.close()

    # Plot the normal mode comparisons
    if os.path.isfile("Normal_Modes/nm.txt") is False:
        current = os.getcwd()
        os.chdir("Normal_Modes")
        output = subprocess.Popen("python {}/Parsers/plot_normal_modes.py > nm.log".format(c["taffi_path"]),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True, encoding='utf8').communicate()[0].strip("\r\n")
        os.chdir(current)

    return

def get_mu(filename):
   with open(filename,'r') as f:
      for lc,lines in enumerate(f):
         fields = lines.split()
         if fields[0] == 'MU': mu = float(fields[1])
   return mu
         

def charmm_minimize_molecule(xyz,polar_flag,ofolder):

   # xyz,ofolder are better in absolute path
   if os.path.isdir(ofolder) is False:
      os.makedirs(ofolder)

   os.chdir(ofolder)

   sublist = ('{} -o single.pdb'.format(xyz)).split()
   xyz_to_pdb.main(sublist)
   
   if polar_flag:
      sublist = ('{}   -o ff  -gens {} -mixing_rule wh --polar -FF '.format(xyz,global_c['gens'])).split()
      sublist.append('{}'.format(global_c['fom_ff']))
      db_to_charmm.main(sublist)
   else:
      sublist = ('{}   -o ff  -gens {} -mixing_rule wh  -FF '.format(xyz,global_c['gens'])).split()
      sublist.append('{}'.format(global_c['fom_ff']))
      db_to_charmm.main(sublist)


   # this inp has been tested with tip3p, see /depot/bsavoie/apps/charmm/test/tip3p_test/polar
   # the result match with /depot/bsavoie/apps/charmm/test/tip3p_test/
   with open('min.inp','w') as f: 
      f.write("* setup structures\n")
      f.write("*\n\n") 

      f.write("! this is necessart to work with non-interger total charge\n")
      f.write("bomlev -1\n")
      f.write("read rtf card name ff/lin.rtf\n")
      f.write("read para card name ff/lin.prm\n")
      f.write("open unit 91 card read name single.pdb\n")
      f.write("read sequ pdb unit 91\n")
      f.write("set residue LIN\n\n")

      f.write("generate @residue first none last none setup warn drude dmass 0.4 !show\n\n")

      f.write("rewind unit 91\n")
      f.write("read coor pdb unit 91\n")
      f.write("coor sdrude\n")
      f.write("coor shake\n")
      f.write("coor print\n\n")


      f.write("set cutoff 999\n")
      f.write("set ctofnb = @cutoff\n")
      f.write("Calc ctofnb = @ctofnb - 2.0\n")
      f.write("set ctonnb = @ctofnb\n")
      f.write("Calc ctonnb = @ctonnb - 1.0\n\n")

      f.write("! shake constraint\n")
      f.write("SHAKE bonh param nofast -\n")
      f.write("select .not. type D* end -\n")
      f.write("select .not. type D* end\n\n")

      f.write("DrudeHardWall L_WALL 0.2 ! set the hard wall at 0.2 Angstrom\n\n")

      f.write("MINI SD nstep 5000 nprint 50 -\n")
      f.write("atom vatom vdistance cdie eps 1.0 vswitch fshift wmin 1.2 -\n")
      f.write("cutnb @cutoff ctofnb @ctofnb ctonnb @ctonnb -\n")
      f.write("e14fac 0.0\n\n") 
   
      f.write('coor orient\n')
      f.write('coor print\n')
      f.write("! First with DSCF to eliminate too high freq from Drude\n")
      f.write("! this should be a more correct setting, but the result seems problematic\n")
      f.write("!VIBRAN\n")
      f.write("!DIAGONALIZE DSCF Finite\n")
      f.write("!END\n\n")
      f.write("! calculate molecular polarizability\n")
      f.write("! This won't work with Diagonalize option on\n")
      f.write("! This is from example in documentation\n")
      f.write("shake off\n")
      f.write("VIBRAN\n")
      f.write("DIAGONALIZE\n")
      f.write("PRINT NORMAL VECTORS DIPOLES SELECT ALL END\n")
      f.write("END\n\n")

      f.write("! calculate dipole\n")
      f.write("COOR DIPOLE\n")
      f.write("SET MU1 = ?RDIP\n\n")

      f.write("OPEN WRITE CARD UNIT 33 NAME mu.log\n")
      f.write("WRITE TITLE UNIT 33\n")
      f.write("* MONOMER:\n")
      f.write("*   MU    @MU1\n")
      f.write("*\n")


      f.write("OPEN WRITE CARD UNIT 10 NAME mm_min.crd\n")
      f.write("WRITE COOR CARD UNIT 10\n")

   shell_command('charmm < min.inp > min.log \n wait','exe')
   sublist = ('mm_min.crd -o minimized.xyz'.split())
   crd_to_xyz.main(sublist)

   return


def shell_command(command,shname):
    with open('{}.sh'.format(shname),'w') as f:
         f.write("#!/bin/bash\n")
         f.write(command)
    command = 'sh {}.sh'.format(shname).split()
    process = subprocess.Popen(command)
    status = process.poll()                                                                                                                                                                        
    while status == None:
       status = process.poll()
       time.sleep(3)

    return

# Wrapper function for the md and figure of merit analysis
def run_md(c,polar_flag=False,shake_flag=False,create_flag=False):

    # Make the md directory if it doesn't already exist
    if os.path.isdir("FOM") is False:
        os.mkdir("FOM")
        
    # Copy the xyz files if necessary
    for i in c["xyz"]:
        if os.path.isfile("FOM/{}".format(i)) is False:
            shutil.copy(i,"FOM/{}".format(i))

    # Generate the single molecule and condensed phase input files
    print("\tgenerating input files for gas-phase and condensed-phase simulations...")
    folder_gas_list=[] 
    folder_list=[] 
    for i in c["xyz"]:

        # Save base name and versions to variables
        name=".".join(i.split('.')[:-1])
        name_sub=name.replace('(','-')
        name_sub=name_sub.replace(')','-')
        name_safe=esc_str(name)
        
        # Grab the number of atoms in this molecule
        with open(i,'r') as f:
            for lc,lines in enumerate(f):
                N_atoms = int(lines.split()[0])
                break

        # Calculate the number of molecules to be run in the condensed phase simulation
        N_mol = 128 

        # Get the simulation temperatures
        T = get_temps(c["fom_temps"],name)
        Vmol = get_rho(c["fom_dens"],name,T)
 
        # If the gas-phase simulation folder doesn't exist
        if os.path.isdir("FOM/{}_gas".format(name)) is False:
            os.mkdir("FOM/{}_gas".format(name))
            
            # Generate FF for CHARMM
            os.chdir("FOM/{}_gas".format(name))
            if os.path.isdir('FF') is False:
               if polar_flag:
                  sublist = ('../../{}   -o FF  -gens {} -mixing_rule wh --polar -FF '.format(i,c['gens'])).split()
                  sublist.append('{}'.format(c['fom_ff']))
                  db_to_charmm.main(sublist)
               else:
                  sublist = ('../../{}   -o FF  -gens {} -mixing_rule wh  -FF '.format(i,c['gens'])).split()
                  sublist.append('{}'.format(c['fom_ff']))
                  db_to_charmm.main(sublist)
            os.chdir("../..")

        # Generate five trajectories at each temperature
        for t in T:
            if os.path.isdir("FOM/{}_gas/{}_K".format(name,t)) is False:
                 os.mkdir("FOM/{}_gas/{}_K".format(name,t))

            for j in range(1,c["fom_n_traj"]+1):
                 if os.path.isdir("FOM/{}_gas/{}_K/run_{}".format(name,t,j)) is False:
                     os.mkdir("FOM/{}_gas/{}_K/run_{}".format(name,t,j)) 
                     os.chdir("FOM/{}_gas/{}_K/run_{}".format(name,t,j))
                     prepare_gas(i,t,j,polar_flag,shake_flag)
                     os.chdir("../../../../")
                 folder_gas_list.append('FOM/{}_gas/{}_K/run_{}'.format(name,t,j))
                        
        # If the condensed-phase simulation folder doesn't exist
        if os.path.isdir("FOM/{}_cond".format(name)) is False:
            os.mkdir("FOM/{}_cond".format(name))

            # Generate FF for CHARMM
            os.chdir("FOM/{}_cond".format(name))
            if os.path.isdir('FF') is False:
               if polar_flag:
                  sublist = ('../../{} -o FF -q {} -gens {} -mixing_rule wh --polar -FF '.format(i,c['charge'],c['gens'])).split()
                  sublist.append('{}'.format(c['fom_ff']))
                  db_to_charmm.main(sublist)
               else:
                  sublist = ('../../{} -o FF -q {} -gens {} -mixing_rule wh -FF '.format(i,c['charge'],c['gens'])).split()
                  sublist.append('{}'.format(c['fom_ff']))
                  db_to_charmm.main(sublist)
            os.chdir("../..")
            
        #Generate five trajectories at each temperature
        for t in T:
            if os.path.isdir("FOM/{}_cond/{}_K".format(name,t)) is False:
                os.mkdir("FOM/{}_cond/{}_K".format(name,t))
            for j in range(1,c["fom_n_traj"]+1):
                """
                ttttt = os.getcwd()
                os.chdir("FOM/{}_cond/{}_K/run_{}".format(name,t,j))
                write_sh(j,t,global_c['fom_md_cond_q'],global_c['fom_md_cond_wt'])
                os.chdir(ttttt)
                """

                
                if os.path.isdir("FOM/{}_cond/{}_K/run_{}".format(name,t,j)) is False:
                    os.mkdir("FOM/{}_cond/{}_K/run_{}".format(name,t,j)) 
                    os.chdir("FOM/{}_cond/{}_K/run_{}".format(name,t,j))
                    prepare_cond(N_mol,i,Vmol[t],t,j,polar_flag,shake_flag)
                    os.chdir("../../../../")
                folder_list.append('FOM/{}_cond/{}_K/run_{}'.format(name,t,j))

    bundle_number = 30
  
    sub_gas_folder = [folder_gas_list[i * bundle_number:(i + 1) * bundle_number] for i in range((len(folder_gas_list) + bundle_number - 1) // bundle_number )]  
    with open('bundle_gas.sh','w') as f:
      f.write("#!/bin/bash\n")
      for fc,sub_folder in enumerate(sub_gas_folder):
         f.write("list{}=(".format(fc))
         for i in sub_folder:
            f.write("'{}' ".format(i))
         f.write(")\n")
         f.write("for i in \"${{list{}[@]}}\"; do\n".format(fc))
         f.write("  cd $i\n")
         f.write("  sbatch gas.sh\n")
         f.write("  cd ../../../../\n")
         f.write("done\n")
         f.write("sleep 5\n")
    subprocess.call("chmod 777 bundle_gas.sh", shell=True)
                        
    sub_cond_folder = [folder_list[i * bundle_number:(i + 1) * bundle_number] for i in range((len(folder_list) + bundle_number - 1) // bundle_number )]  
    with open('bundle.sh','w') as f:
      f.write("#!/bin/bash\n")
      f.write("count=1\n")
      for fc,sub_folder in enumerate(sub_cond_folder):
         f.write("list{}=(".format(fc))
         for i in sub_folder:
            f.write("'{}' ".format(i))
         f.write(")\n")
         f.write("for i in \"${{list{}[@]}}\"; do\n".format(fc))
         f.write("  screen -X -S ${count}_FOM quit \n")
         f.write("  screen -dmS ${count}_FOM\n")
         f.write("  screen -S ${count}_FOM -X stuff \"ana3 ; cd $i ; sh extend.sh &> extend.log \\n\" \n")
         f.write("  count=$((count+1)) \n")
         #f.write("  cd $i\n")
         #f.write("  sh extend.sh &> extend.log &\n")
         #f.write("  cd ../../../../\n")
         f.write("done\n")
         f.write("sleep 5\n")
      f.write("wait\n")
    subprocess.call("chmod 777 bundle.sh", shell=True)
    with open('bundle_total.sh','w') as f:
      f.write("#!/bin/bash\n")
      f.write("count=1\n")
      f.write("list=(")
      for i in folder_list:
         f.write("'{}' ".format(i))
      f.write(")\n")
      f.write("for i in \"${list[@]}\"; do\n")
      f.write("  screen -X -S ${count}_FOM quit \n")
      f.write("  screen -dmS ${count}_FOM\n")
      f.write("  screen -S ${count}_FOM -X stuff \"ana3 ; cd $i ; sh extend.sh &> extend.log \\n\" \n")
      f.write("  count=$((count+1)) \n")
      #f.write("  cd $i\n")
      #f.write("  sh extend.sh &> extend.log &\n")
      #f.write("  cd ../../../../\n")
      f.write("done\n")
      f.write("wait\n")
    subprocess.call("chmod 777 bundle_total.sh", shell=True)
                    
    if create_flag:
      quit()


    # Combine the output trajectories
    print("\tcombining trajectories...")
    #N_start = int(np.ceil(float(c["fom_t_equil"])/float(c["fom_t_ext"])))   # N_start is based upon how long the user wants a pre-equilibration
    #N_end = N - 1
    N_start = 2
    N_mol = 128
    time_per_trj = global_c['fom_t_ext']
    bundle_string = ""
    bundle_count = 0
    b_count = 0
    jobids = []
    with open('track.txt','w') as f:
        f.write("path      N_trj\n")
    for i in c["xyz"]:

        # Save base name and versions to variables
        name=".".join(i.split('.')[:-1])
        name_sub=name.replace('(','-')
        name_sub=name_sub.replace(')','-')
        name_safe=esc_str(name)

        # Get the simulation temperatures
        T = get_temps(c["fom_temps"],name)

        # Loop over temperatures and trajectories
        for t in T:
            for j in range(1,c["fom_n_traj"]+1):
                current_dir = os.getcwd()

                # Find the last complete trajectory
                os.chdir("FOM/{}_cond/{}_K/run_{}".format(name,t,j))
                pattern = '*.sys.crd'
                trjfile = [ files for files in os.listdir('.') if fnmatch.fnmatch(files,pattern) ] 
                os.chdir(current_dir)
                if len(trjfile) < 3:
                  with open(current_dir+'/track.txt','a') as f:
                     f.write("{:<40s} fail\n".format('{}_cond/{}_K/run_{}'.format(name,t,j)))
                else:
                  trjfile = sorted(trjfile,key = crd_num)
                  N_end = int(trjfile[-1].split('.')[0])
                  with open(current_dir+'/track.txt','a') as f:
                     f.write("{:<40s} {:>4d}\n".format('{}_cond/{}_K/run_{}'.format(name,t,j),N_end))
                  #os.chdir("FOM/{}_cond/{}_K/run_{}/analysis".format(name,t,j))
                  #polar_tmp(i,t)
                  #os.chdir(current_dir)

                # Combine the production trajectories 
                if os.path.isdir("FOM/{}_cond/{}_K/run_{}/analysis".format(name,t,j)) is False:
                    if len(trjfile) >= 3:
                       os.makedirs("FOM/{}_cond/{}_K/run_{}/analysis".format(name,t,j))
                       os.chdir("FOM/{}_cond/{}_K/run_{}/analysis".format(name,t,j))
                       prepare_analysis(c,t,i,N_mol,N_end,N_start,time=time_per_trj,polar_flag=polar_flag)
                       os.chdir(current_dir)
                       substring="sh exe.sh &\n" ## 10 ps
                       substring=substring.format(c["taffi_path"],t)
                       bundle_string+="\ncd {}/FOM/{}_cond/{}_K/run_{}/analysis\n".format(os.getcwd(),name_safe,t,j)+substring

                       # Assemble the submission string for the batch submission file
                       b_count += 1

                if os.path.isfile("FOM/{}_gas/{}_K/run_{}/energy.txt".format(name,t,j)) is False:
                    substring="{}/FF_functions/read_charmm_energy.py -log gas.log --silent --gas &\n" ## 10 ps
                    substring=substring.format(c["taffi_path"],t)
                    bundle_string+="\ncd {}/FOM/{}_gas/{}_K/run_{} \n".format(os.getcwd(),name_safe,t,j)+substring
                    # Assemble the submission string for the batch submission file
                    b_count += 1

                if b_count == int(c["fom_geoopt_ppn"]):
                    bundle_string+="\nwait\n"
                    shell_submit(bundle_string,"analysis_bundle_{}.submit".format(bundle_count),c["fom_geoopt_ppn"],c["fom_geoopt_wt"],c["fom_geoopt_q"],c["fom_geoopt_ppn"],batch=c["batch"],run_dir="{}".format(os.getcwd()))
                    output = subprocess.Popen([c["sub_cmd"],"analysis_bundle_{}.submit".format(bundle_count)],stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0].strip("\r\n")
                    output = output.split()[-1]
                    jobids += [output]
                    bundle_string = ""
                    bundle_count += 1
                    b_count = 0
    # submit remainder of jobs
    if b_count != 0:
        bundle_string+="\nwait\n"
        shell_submit(bundle_string,"analysis_bundle_{}.submit".format(bundle_count),c["fom_geoopt_ppn"],c["fom_geoopt_wt"],c["fom_geoopt_q"],c["fom_geoopt_ppn"],batch=c["batch"],run_dir="{}".format(os.getcwd()))
        output = subprocess.Popen([c["sub_cmd"],"analysis_bundle_{}.submit".format(bundle_count)],stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0].strip("\r\n")
        output = output.split()[-1]
        jobids += [output]
        bundle_string = ""
        bundle_count += 1
        b_count = 0

    # Wait until the parses complete
    user = getpass.getuser()
    monitor_jobs(jobids,user)
                    
    # Parse the diffusion coefficients and dielectric constants
    print("\tparsing diffusion and dielectric constant from trajectories...")
    bundle_string = ""
    bundle_count = 0
    b_count = 0
    jobids = []
    for i in c["xyz"]:

        # Save base name and versions to variables
        name=".".join(i.split('.')[:-1])
        name_sub=name.replace('(','_')
        name_sub=name_sub.replace(')','_')
        name_safe=esc_str(name)

        # Get the simulation temperatures
        T = get_temps(c["fom_temps"],name)

        # Loop over temperatures and trajectories
        for t in T:
            for j in range(1,c["fom_n_traj"]+1):
                current_dir = os.getcwd()

                if os.path.isdir("FOM/{}_cond/{}_K/run_{}/analysis".format(name,t,j)) is False: continue

                # Parse diffusivities

                # Check completeness, remove incomplete folders 
                if os.path.isdir("FOM/{}_cond/{}_K/run_{}/analysis/msd".format(name,t,j)) is True and os.path.isfile("FOM/{}_cond/{}_K/run_{}/analysis/msd/avg-0.txt".format(name,t,j)) is False:
                    shutil.rmtree("FOM/{}_cond/{}_K/run_{}/analysis/msd".format(name,t,j))

                if os.path.isfile("FOM/{}_cond/{}_K/run_{}/analysis/msd/avg-0.txt".format(name,t,j)) is False:

                    # Assemble the submission string for the batch submission file
                    substring="python3 {}/Parsers/msd_parse.lin.py -mols 0 -o msd -every 10 &\n"
                    substring=substring.format(c["taffi_path"],t)
                    bundle_string+="\ncd {}/FOM/{}_cond/{}_K/run_{}/analysis\n".format(os.getcwd(),name_safe,t,j)+substring

                    # Assemble the submission string for the batch submission file
                    b_count += 1
                    if b_count == int(c["fom_geoopt_ppn"]):
                        bundle_string+="\nwait\n"
                        shell_submit(bundle_string,"msd_eps_bundle_{}.submit".format(bundle_count),c["fom_geoopt_ppn"],c["fom_geoopt_wt"],c["fom_geoopt_q"],c["fom_geoopt_ppn"],batch=c["batch"],run_dir="{}".format(os.getcwd()))
                        output = subprocess.Popen([c["sub_cmd"],"msd_eps_bundle_{}.submit".format(bundle_count)],stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0].strip("\r\n")
                        output = output.split()[-1]
                        jobids += [output]
                        bundle_string = ""
                        bundle_count += 1
                        b_count = 0

                # Parse dielectric constants
                if os.path.isfile("FOM/{}_cond/{}_K/run_{}/analysis/traj_dipole.txt".format(name,t,j)) is False:

                    # Assemble the submission string for the batch submission file
                    #substring="python {}/Parsers/dielectric_parse_new.py -map run_{}.map -traj combined.lammpstrj -every 1 -T {} -mols 0 -folder eps &\n"
                    substring="{}/Parsers/dielectric_parse.lin.py -Temp {} &\n" ## 10 ps
                    substring=substring.format(c["taffi_path"],t)
                    bundle_string+="\ncd {}/FOM/{}_cond/{}_K/run_{}/analysis\n".format(os.getcwd(),name_safe,t,j)+substring

                    # Assemble the submission string for the batch submission file
                    b_count += 1
                    if b_count == int(c["fom_geoopt_ppn"]):
                        bundle_string+="\nwait\n"
                        shell_submit(bundle_string,"msd_eps_bundle_{}.submit".format(bundle_count),c["fom_geoopt_ppn"],c["fom_geoopt_wt"],c["fom_geoopt_q"],c["fom_geoopt_ppn"],batch=c["batch"],run_dir="{}".format(os.getcwd()))
                        output = subprocess.Popen([c["sub_cmd"],"msd_eps_bundle_{}.submit".format(bundle_count)],stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0].strip("\r\n")
                        output = output.split()[-1]
                        jobids += [output]
                        bundle_string = ""
                        bundle_count += 1
                        b_count = 0
    
    # submit remainder of jobs
    if b_count != 0:
        bundle_string+="\nwait\n"
        shell_submit(bundle_string,"msd_eps_bundle_{}.submit".format(bundle_count),c["fom_geoopt_ppn"],c["fom_geoopt_wt"],c["fom_geoopt_q"],c["fom_geoopt_ppn"],batch=c["batch"],run_dir="{}".format(os.getcwd()))
        output = subprocess.Popen([c["sub_cmd"],"msd_eps_bundle_{}.submit".format(bundle_count)],stdout=subprocess.PIPE,stderr=subprocess.PIPE, encoding='utf8').communicate()[0].strip("\r\n")
        output = output.split()[-1]
        jobids += [output]
        bundle_string = ""
        bundle_count += 1
        b_count = 0

    # Wait until the parses complete
    user = getpass.getuser()
    monitor_jobs(jobids,user)

    """
    # Check parse for completion
    i_msd = []
    i_eps = []
    for i in c["xyz"]:

        # Save base name and versions to variables
        name=".".join(i.split('.')[:-1])

        # Get the simulation temperatures
        T = get_temps(c["fom_temps"],name)

        # Loop over temperatures and trajectories
        for t in T:
            for j in range(1,c["fom_n_traj"]+1):
                current_dir = os.getcwd()

                # Check completeness, remove incomplete folders 
                if os.path.isdir("FOM/{}_cond/{}_K/run_{}/msd".format(name,t,j)) is True and os.path.isfile("FOM/{}_cond/{}_K/run_{}/msd/avg-0.txt".format(name,t,j)) is False:
                    i_msd += ["FOM/{}_cond/{}_K/run_{}/msd".format(name,t,j)]

                if os.path.isdir("FOM/{}_cond/{}_K/run_{}/eps".format(name,t,j)) is True and os.path.isfile("FOM/{}_cond/{}_K/run_{}/eps/0_ac.txt".format(name,t,j)) is False:
                    i_eps += ["FOM/{}_cond/{}_K/run_{}/eps".format(name,t,j)]

    # Exit if there are incompletes
    if len(i_msd) > 0:
        print("\nThe following msd parses failed:\n\n\t{}".format("\n\t".join([ str(_) for _ in i_msd ])))
    if len(i_eps) > 0:
        print("\nThe following eps parses failed:\n\n\t{}".format("\n\t".join([ str(_) for _ in i_eps ])))        
    if len(i_msd) > 0 or len(i_eps) > 0:
        quit()
    """

    
    # Clean up submission files
    #for f in glob.glob("*.submit*"):
    #    os.remove(f)
    for f in glob.glob("*.err"):
        os.remove(f)
    for f in glob.glob("*.out"):
        os.remove(f)
    

    # Parse thermodynamic figures of merit
    print("\tparsing thermodynamic averages...")
    kb=0.001987204118                           # Boltzman constant in units of kcal/mol K
    kcal_mol_to_J = 4.184 * 1000.0 / 6.0221409E23 # employed for converting the isothermal compressibility to useful units.  
    kcal_to_kJ = 4.184
    pressure = 1.0                    # for thermal expansion coeff (atm)
    atm_to_pa = 101325.0              # atmospheres to pascals
    PV_const = pressure * 1E-30 * atm_to_pa / kcal_mol_to_J # A^-3 to m^-3, conversion to pa, conversion to kcal/mol
    Data = {}
    Data_avg = {}
    if os.path.isdir("FOM/Results") is False:
        os.mkdir("FOM/Results")

    # Create header of thermo file
    """
    if os.path.isfile("FOM/Results/thermo.txt") is False:
        with open("FOM/Results/thermo.txt",'w') as f:
            f.write("{:<40s} {:<10s} {:<16s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s}\n"\
                    .format("compound","T(k)","d(g/cm^3)","std(g/cm^3)","H_val(kJ/mol)","std(kJ/mol)","eps","std","iso(GP^-1)","std(GP^-1)","t_exp(10^-3/K)","std(10^-3/K)"))
            write_flag = True
    else:
        write_flag = False
    """
    with open("FOM/Results/thermo.txt",'w') as f:
         f.write("{:<40s} {:<10s} {:<16s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s}\n"\
                 .format("compound","T(k)","d(g/cm^3)","std(g/cm^3)","H_val(kJ/mol)","std(kJ/mol)","eps","std","iso(GP^-1)","std(GP^-1)","t_exp(10^-3/K)","std(10^-3/K)","eps_inf","std","mobility(10^-5 cm^2/s)","std"))
         write_flag = True

    with open("FOM/Results/thermo_detail.txt",'w') as f:
         f.write("{:<40s} {:<10s} {:<10s} {:<16s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s} {:<15s}\n"\
                 .format("compound","T(k)","ntrj","d(g/cm^3)","std(g/cm^3)","H_val(kJ/mol)","std(kJ/mol)","eps","std","iso(GP^-1)","std(GP^-1)","t_exp(10^-3/K)","std(10^-3/K)","eps_inf","std","mobility(10^-5 cm^2/s)","std"))
         write_flag_detail = True

    # For each molecule that was simulated, parse the density, dielectric, diffusivity, thermal exp coeff, and bulk modulus
    for i in c["xyz"]:

        molw = get_molw(i)

        # Save base name and versions to variables
        name=".".join(i.split('.')[:-1])
        name_sub=name.replace('(','-')
        name_sub=name_sub.replace(')','-')
        name_safe=esc_str(name)

        Data[name] = {}
        Data_avg[name] = {}

        # Grab the number of atoms in this molecule
        with open(i,'r') as f:
            for lc,lines in enumerate(f):
                N_atoms = int(lines.split()[0])
                break

        # Calculate the number of molecules in the simulation
        N_mol = 128
        N_atoms = N_mol * N_atoms

        # Get the simulation temperatures
        T = get_temps(c["fom_temps"],name)

        # Loop over temperatures and trajectories
        for t in T:
            Data[name][t] = {}
            for j in range(1,c["fom_n_traj"]+1):

                if os.path.isfile("FOM/{}_cond/{}_K/run_{}/analysis/energy.txt".format(name,t,j)) is False: continue
                # Read data into dictionary
                Data[name][t][j] = {}
                Data[name][t][j]["U"],Data[name][t][j]["V"] = get_cols("FOM/{}_cond/{}_K/run_{}/analysis/energy.txt".format(name,t,j),[3,6],1)
                Data[name][t][j]["U_gas"] = get_cols("FOM/{}_gas/{}_K/run_{}/energy.txt".format(name,t,j),[3],1)
                Data[name][t][j]["t_msd"],Data[name][t][j]["msd"] = get_cols("FOM/{}_cond/{}_K/run_{}/analysis/msd/avg-0.txt".format(name,t,j),[0,1],1)
                Data[name][t][j]["t_eps"],Data[name][t][j]["eps_ravg"] = get_cols("FOM/{}_cond/{}_K/run_{}/analysis/traj_dipole.txt".format(name,t,j),[0,3],1) 
                Data[name][t][j]["eps"] = Data[name][t][j]["eps_ravg"][-1]
                Data[name][t][j]["eps_inf"] = 1 


                # Convert to arrays
                #Data[name][t][j]["d"] = np.array(Data[name][t][j]["d"])
                Data[name][t][j]["U"] = np.array(Data[name][t][j]["U"])
                Data[name][t][j]["V"] = np.array(Data[name][t][j]["V"])
                Data[name][t][j]["U_gas"] = np.array(Data[name][t][j]["U_gas"])
                Data[name][t][j]["t_msd"] = np.array(Data[name][t][j]["t_msd"])
                Data[name][t][j]["msd"] = np.array(Data[name][t][j]["msd"])
                Data[name][t][j]["t_eps"] = np.array(Data[name][t][j]["t_eps"])
                Data[name][t][j]["eps_ravg"] = np.array(Data[name][t][j]["eps_ravg"])


                # Calculate molecular volume
                Data[name][t][j]["Vmol"] = Data[name][t][j]["V"]/N_mol
                Data[name][t][j]["d"] = vmol_to_dens(Data[name][t][j]["Vmol"],molw)

                # Calculate corrected eps 
                if polar_flag:
                     Data[name][t][j]['polarizability'] = get_polar("FOM/{}_cond/{}_K/run_{}/analysis/polar.log".format(name,t,j))
                     # eps_inf based on Clausius– Mossotti eq, for nonpolar this is 1
                     correct = 4*np.pi*Data[name][t][j]['polarizability']/(3*np.mean(Data[name][t][j]["Vmol"]))
                     eps_inf = (-1-2*correct)/(correct-1)
                     Data[name][t][j]["eps"] = Data[name][t][j]["eps"]-1+eps_inf
                     Data[name][t][j]["eps_inf"] = eps_inf

                # Calculate isothermal compressibility
                Data[name][t][j]["iso_comp"] = np.std(Data[name][t][j]["V"])**(2.0)/(np.mean(Data[name][t][j]["V"])*kb*float(t))*1E-30*1E9/kcal_mol_to_J 

                # # Calculate thermal expansion coefficient (OLD)
                # tmp_dev_V = Data[name][t][j]["V"] - mean(Data[name][t][j]["V"])
                # tmp_dev_H = Data[name][t][j]["H"] - mean(Data[name][t][j]["H"])
                # Data[name][t][j]["therm_exp"] = mean(tmp_dev_V*tmp_dev_H)/(kb*float(t)**(2.0)*mean(Data[name][t][j]["V"]))*1E3

                # Calculate thermal expansion coefficient
                Data[name][t][j]["therm_exp"] = ( np.mean(Data[name][t][j]["V"]*Data[name][t][j]["U"]) - np.mean(Data[name][t][j]["V"])*np.mean(Data[name][t][j]["U"]) + PV_const * np.std(Data[name][t][j]["V"])**(2.0) ) \
                                         / (kb*float(t)**(2.0)*np.mean(Data[name][t][j]["V"])) * 1E3


                print("{} {} {}: {} {} {} {}".format(name,t,j,Data[name][t][j]["iso_comp"],Data[name][t][j]["therm_exp"],Data[name][t][j]["eps"],Data[name][t][j]["eps_inf"]))

                # Save thermodynamic mean and std for each run to a dedicated mean_thermo.txt file
                if os.path.isfile("FOM/{}_cond/{}_K/run_{}/mean_thermo.txt".format(name,t,j)) is True:
                    os.remove("FOM/{}_cond/{}_K/run_{}/mean_thermo.txt".format(name,t,j))
                if os.path.isfile("FOM/{}_cond/{}_K/run_{}/mean_thermo.txt".format(name,t,j)) is False:
                    with open("FOM/{}_cond/{}_K/run_{}/mean_thermo.txt".format(name,t,j),'w') as f:
                        f.write("{:<20s}  {:<20s} {:<20s}\n".format("property","mean","stdev"))
                        #f.write("{:<20s} {:< 20.8f} {:< 20.8f}\n".format("Vmol",np.mean(Data[name][t][j]["Vmol"]),np.std(Data[name][t][j]["Vmol"])))
                        f.write("{:<20s} {:< 20.8f} {:< 20.8f}\n".format("density",np.mean(Data[name][t][j]["d"]),np.std(Data[name][t][j]["d"])))
                        f.write("{:<20s} {:< 20.8f} {:< 20.8f}\n".format("U",np.mean(Data[name][t][j]["U"]),np.std(Data[name][t][j]["U"])))

            # Make plot of average internal energy vs time across the trajectories
            if os.path.isfile("FOM/Results/{}_cond_{}_K_U_composite.png".format(name,t)) is False:
                make_plot("FOM/Results/{}_cond_{}_K_U_composite.png".format(name,t),\
                          [ np.arange(len(Data[name][t][j]["U"]))*0.001 for j in sorted(Data[name][t].keys()) ],\
                          [ Data[name][t][j]["U"] for j in sorted(Data[name][t].keys()) ],\
                          x_label="Time (ns)",y_label="U (kcal/mol)")                      

            # Make plot of density vs time across the trajectories
            """
            if os.path.isfile("FOM/Results/{}_cond_{}_K_vmol_composite.png".format(name,t)) is False:
                make_plot("FOM/Results/{}_cond_{}_K_vmol_composite.png".format(name,t),\
                          [ np.arange(len(Data[name][t][j]["Vmol"]))*0.001 for j in sorted(Data[name][t].keys()) ],\
                          [ Data[name][t][j]["Vmol"] for j in sorted(Data[name][t].keys()) ],\
                          x_label="Time (ns)",y_label="d (g/cm^3)")                      
            """
            if os.path.isfile("FOM/Results/{}_cond_{}_K_d_composite.png".format(name,t)) is False:
                make_plot("FOM/Results/{}_cond_{}_K_d_composite.png".format(name,t),\
                          [ np.arange(len(Data[name][t][j]["d"]))*0.001 for j in sorted(Data[name][t].keys()) ],\
                          [ Data[name][t][j]["d"] for j in sorted(Data[name][t].keys()) ],\
                          x_label="Time (ns)",y_label="d (g/cm^3)")                      

            # Make plot of msd vs time across the trajectories
            if os.path.isfile("FOM/Results/{}_cond_{}_K_msd_composite.png".format(name,t)) is False:
                make_plot("FOM/Results/{}_cond_{}_K_msd_composite.png".format(name,t),\
                          [ Data[name][t][j]["t_msd"] for j in sorted(Data[name][t].keys()) ],\
                          [ Data[name][t][j]["msd"] for j in sorted(Data[name][t].keys()) ],\
                          x_label="Time (ns)",y_label="MSD (A^2)")                      
            
            # Make plot of running average of eps vs time across the trajectories
            if os.path.isfile("FOM/Results/{}_cond_{}_K_eps_composite.png".format(name,t)) is False:
                make_plot("FOM/Results/{}_cond_{}_K_eps_composite.png".format(name,t),\
                          [ Data[name][t][j]["t_eps"] for j in sorted(Data[name][t].keys()) ],\
                          [ Data[name][t][j]["eps_ravg"] for j in sorted(Data[name][t].keys()) ],\
                          x_label="Time (ns)",y_label="eps")                      

        # Save averages of important quatities
        for t in T:
            for j in range(1,c["fom_n_traj"]+1):
                if os.path.isfile("FOM/{}_cond/{}_K/run_{}/analysis/energy.txt".format(name,t,j)) is False: continue
                Data[name][t][j]["Vmol"] = np.mean(Data[name][t][j]["Vmol"])
                Data[name][t][j]["d"] = np.mean(Data[name][t][j]["d"])
                Data[name][t][j]["H_vap"] = ( np.mean(Data[name][t][j]["U_gas"]) - np.mean(Data[name][t][j]["U"])/float(N_mol) + kb*float(t) ) * kcal_to_kJ


        # Calculate averages over the trajectories
        #keys = ["Vmol","H_vap","iso_comp","therm_exp","eps","t_msd","msd"]
        #keys = ["Vmol","H_vap","iso_comp","therm_exp","eps","eps_inf"]
        keys = ["d","H_vap","iso_comp","therm_exp","eps","eps_inf","t_msd","msd"]
        for t in T:
            real_ntrj = 0
            Data_avg[name][t] = {}
            for j in range(1,c["fom_n_traj"]+1):
                if os.path.isfile("FOM/{}_cond/{}_K/run_{}/analysis/energy.txt".format(name,t,j)) is False: continue
                real_ntrj += 1
                
                Data[name][t][j]["msd_slope_lin"],Data[name][t][j]["msd_std_lin"] = calc_derivative(Data[name][t][j]["t_msd"],Data[name][t][j]["msd"])
                Data[name][t][j]["msd_slope_log"],Data[name][t][j]["msd_std_log"] = calc_derivative(Data[name][t][j]["t_msd"],Data[name][t][j]["msd"],logx_opt=True,logy_opt=True)
                Data[name][t][j]["mobility"],Data[name][t][j]["mobility_err"],Data[name][t][j]["mobility_inds"] = calc_mobility(np.array(Data[name][t][j]["msd_slope_log"]),np.array(Data[name][t][j]["msd_slope_lin"]),start=100) 
                with open("FOM/Results/thermo_detail.txt",'a') as f:                
                    f.write("{:<40s} {:<10s} {:<10s} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f}\n"\
                            .format(name,t,str(real_ntrj),Data[name][t][j]["d"],0.0,Data[name][t][j]["H_vap"],0.0,\
                                    Data[name][t][j]["eps"],0.0,Data[name][t][j]["iso_comp"],0.0,Data[name][t][j]["therm_exp"],0.0,Data[name][t][j]["eps_inf"],0.0,Data[name][t][j]["mobility"],0.0))
                for k in keys:

                    if real_ntrj == 1:
                        Data_avg[name][t][k] = Data[name][t][j][k]
                        Data_avg[name][t][k+"_std"] = Data[name][t][j][k]**(2.0)
                    else:
                        if k == 't_msd' or k == 'msd':
                           o_datalen = Data_avg[name][t][k].shape[0]
                           n_datalen = Data[name][t][j][k].shape[0]
                           if o_datalen > n_datalen:
                              Data_avg[name][t][k] = Data_avg[name][t][k][:n_datalen] + Data[name][t][j][k]
                              Data_avg[name][t][k+"_std"] = Data_avg[name][t][k+"_std"][:n_datalen] + Data[name][t][j][k]**(2.0)
                           else:
                              Data_avg[name][t][k] = Data_avg[name][t][k] + Data[name][t][j][k][:o_datalen]
                              Data_avg[name][t][k+"_std"] = Data_avg[name][t][k+"_std"] + Data[name][t][j][k][:o_datalen]**(2.0)
                        else:
                           Data_avg[name][t][k] = Data_avg[name][t][k] + Data[name][t][j][k]
                           Data_avg[name][t][k+"_std"] = Data_avg[name][t][k+"_std"] + Data[name][t][j][k]**(2.0)
            
            # Calculate average and std_dev
            for k in keys:
                if real_ntrj == 0: continue
                #Data_avg[name][t][k] = Data_avg[name][t][k] / float(c["fom_n_traj"])
                Data_avg[name][t][k] = Data_avg[name][t][k] / float(real_ntrj)
                #Data_avg[name][t][k+"_std"] = ( Data_avg[name][t][k+"_std"] / float(c["fom_n_traj"]) - Data_avg[name][t][k]**(2.0) )**(0.5)
                Data_avg[name][t][k+"_std"] = ( np.abs(Data_avg[name][t][k+"_std"] / float(real_ntrj) - Data_avg[name][t][k]**(2.0)) )**(0.5) ## NEED TO BE VALIDATED

            # Check convergence with respect to energy
            """
            for j in range(1,c["fom_n_traj"]+1):
                if j == 1:
                    Data_avg[name][t]["U"] = Data[name][t][j]["U"]
                else:
                    Data_avg[name][t]["U"] = Data_avg[name][t]["U"] + Data[name][t][j]["U"]
            Data_avg[name][t]["U"] = Data_avg[name][t]["U"] / float(c["fom_n_traj"])
            if abs((np.mean(Data_avg[name][t]["U"][:1000]) - np.mean(Data_avg[name][t]["U"][-1000:]))/float(N_atoms)) > 0.01:
                print("{} {} didn't converge ({})".format(name,t,abs((np.mean(Data_avg[name][t]["U"][:1000]) - np.mean(Data_avg[name][t]["U"][-1000:]))/float(N_atoms))))
            """

            if real_ntrj == 0: continue
            

            # Calculate the derivatives
            Data_avg[name][t]["msd_slope_lin"],Data_avg[name][t]["msd_std_lin"] = calc_derivative(Data_avg[name][t]["t_msd"],Data_avg[name][t]["msd"])
            Data_avg[name][t]["msd_slope_log"],Data_avg[name][t]["msd_std_log"] = calc_derivative(Data_avg[name][t]["t_msd"],Data_avg[name][t]["msd"],logx_opt=True,logy_opt=True)
            Data_avg[name][t]["mobility"],Data_avg[name][t]["mobility_err"],Data_avg[name][t]["mobility_inds"] = calc_mobility(np.array(Data_avg[name][t]["msd_slope_log"]),np.array(Data_avg[name][t]["msd_slope_lin"]),start=100) 
#            print "{} {}: {} {} {} {} {}".format(name,t,Data_avg[name][t]["mobility"],Data_avg[name][t]["mobility_err"],Data_avg[name][t]["mobility_inds"],Data_avg[name][t]["therm_exp"],Data_avg[name][t]["iso_comp"])

            if write_flag is True:
                with open("FOM/Results/thermo.txt",'a') as f:                
                    f.write("{:<40s} {:<10s} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f} {:< 15.6f}\n"\
                            #.format(name,t,Data_avg[name][t]["Vmol"],Data_avg[name][t]["Vmol_std"],Data_avg[name][t]["H_vap"],Data_avg[name][t]["H_vap_std"],\
                            #        Data_avg[name][t]["eps"],Data_avg[name][t]["eps_std"],Data_avg[name][t]["iso_comp"],Data_avg[name][t]["iso_comp_std"],Data_avg[name][t]["therm_exp"],Data_avg[name][t]["therm_exp_std"],Data_avg[name][t]["eps_inf"],Data_avg[name][t]["eps_inf_std"],Data_avg[name][t]["mobility"],Data_avg[name][t]["mobility_err"]))
                            .format(name,t,Data_avg[name][t]["d"],Data_avg[name][t]["d_std"],Data_avg[name][t]["H_vap"],Data_avg[name][t]["H_vap_std"],\
                                    Data_avg[name][t]["eps"],Data_avg[name][t]["eps_std"],Data_avg[name][t]["iso_comp"],Data_avg[name][t]["iso_comp_std"],Data_avg[name][t]["therm_exp"],Data_avg[name][t]["therm_exp_std"],Data_avg[name][t]["eps_inf"],Data_avg[name][t]["eps_inf_std"],Data_avg[name][t]["mobility"],Data_avg[name][t]["mobility_err"]))

        """
        # Calculate the derivatives
        
        for t in T:
            Data_avg[name][t]["msd_slope_lin"],Data_avg[name][t]["msd_std_lin"] = calc_derivative(Data_avg[name][t]["t_msd"],Data_avg[name][t]["msd"])
            Data_avg[name][t]["msd_slope_log"],Data_avg[name][t]["msd_std_log"] = calc_derivative(Data_avg[name][t]["t_msd"],Data_avg[name][t]["msd"],logx_opt=True,logy_opt=True)
            Data_avg[name][t]["mobility"],Data_avg[name][t]["mobility_err"],Data_avg[name][t]["mobility_inds"] = calc_mobility(np.array(Data_avg[name][t]["msd_slope_log"]),np.array(Data_avg[name][t]["msd_slope_lin"]),start=100) 
            print("{} {}: {} {} {} {} {}".format(name,t,Data_avg[name][t]["mobility"],Data_avg[name][t]["mobility_err"],Data_avg[name][t]["mobility_inds"],Data_avg[name][t]["therm_exp"],Data_avg[name][t]["iso_comp"]))
        """
    quit()

def vmol_to_dens(vmol,molw):
    
    #       kg->g   g-->mol    mol -> #      m^3->A^3
    Const = 1000  * 1/molw   *  (6.02*1e23) *  1e-30    # unit: #/A^3
    dens = 1/(vmol*Const*1000) # unit: g/cm^3

    return dens

def get_molw(filename):

    mass_dict = {'H':1.00794,'He':4.002602,'Li':6.941,'Be':9.012182,'B':10.811,'C':12.011,'N':14.00674,'O':15.9994,'F':18.9984032,'Ne':20.1797,\
                 'Na':22.989768,'Mg':24.3050,'Al':26.981539,'Si':28.0855,'P':30.973762,'S':32.066,'Cl':35.4527,'Ar':39.948,\
                 'K':39.0983,'Ca':40.078,'Sc':44.955910,'Ti':47.867,'V':50.9415,'Cr':51.9961,'Mn':54.938049,'Fe':55.845,'Co':58.933200,'Ni':58.6934,'Cu':63.546,'Zn':65.39,\
                 'Ga':69.723,'Ge':72.61,'As':74.92159,'Se':78.96,'Br':79.904,'Kr':83.80,\
                 'Rb':85.4678,'Sr':87.62,'Y':88.90585,'Zr':91.224,'Nb':92.90638,'Mo':95.94,'Tc':98.0,'Ru':101.07,'Rh':102.90550,'Pd':106.42,'Ag':107.8682,'Cd':112.411,\
                 'In':114.818,'Sn':118.710,'Sb':121.760,'Te':127.60,'I':126.90447,'Xe':131.29,\
                 'Cs':132.90545,'Ba':137.327,'La':138.9055,'Hf':178.49,'Ta':180.9479,'W':183.84,'Re':186.207,'Os':190.23,'Ir':192.217,'Pt':195.078,'Au':196.96655,'Hg':200.59,\
                 'Tl':204.3833,'Pb':207.2,'Bi':208.98038,'Po':209.0,'At':210.0,'Rn':222.0}
    
    with open(filename,'r') as f:
         flag = 0
         atom_count = 0
         molw = 0
         for lc,lines in enumerate(f):
            fields = lines.split()  
            if lc == 0: 
               natoms = int(fields[0])
               continue
            if lc == 1: 
               flag =1   
               continue
            if flag == 1:
               if atom_count == natoms: break
               if fields[0] not in list(mass_dict.keys()):
                  print("Unknowm atom: {} in the {} file, Exiting...".format(fields[0],filename))
                  quit()
               molw += mass_dict[fields[0]]
               atom_count += 1
    print("{:} has molecular weight of {:3.5f} g/mol".format(filename,molw))

    return molw

def crd_num(x):
    return(int(x.split('.')[0]))

def get_polar(filename):
    with open(filename,'r') as f:
        flag = 0
        for lc,lines in enumerate(f):
            fields = lines.split()
            if len(fields) < 4 and flag == 0: continue
          

            if flag == 0  and fields[0] == 'TOTAL' and fields[1] == 'POLARIZABILITY' :
               flag = 1
               continue

            if flag == 1 and len(fields) == 5:
               polar = float(fields[4])
               flag = 0 
               break
    return polar

   
def polar_tmp(xyz,temp):
   sublist = ('../../../../{} -o single.pdb'.format(xyz)).split()
   xyz_to_pdb.main(sublist)
   with open('polar.inp','w') as f: 
      f.write("* setup structures\n")
      f.write("*\n\n") 

      f.write("! this is necessart to work with non-interger total charge\n")
      f.write("bomlev -1\n")
      f.write("set temp = {}\n".format(temp))
      f.write("read rtf card name ../lin.rtf\n")
      f.write("read para card name ../lin.prm\n")
      f.write("open unit 91 card read name single.pdb\n")
      f.write("read sequ pdb unit 91\n")
      f.write("set residue LIN\n\n")

      f.write("generate @residue first none last none setup warn drude dmass 0.4 !show\n\n")

      f.write("rewind unit 91\n")
      f.write("read coor pdb unit 91\n")
      f.write("coor sdrude\n")
      f.write("coor shake\n")
      f.write("coor print\n\n")

      f.write("! Dimension of the solvent volume\n")
      f.write("set 7 = 9999\n")
      f.write("set 8 @7\n")
      f.write("set 9 @8\n\n")

      f.write("set cutoff 999\n")
      f.write("set ctofnb = @cutoff\n")
      f.write("Calc ctofnb = @ctofnb - 2.0\n")
      f.write("set ctonnb = @ctofnb\n")
      f.write("Calc ctonnb = @ctonnb - 1.0\n\n")

      f.write("! shake constraint\n")
      f.write("SHAKE bonh param nofast -\n")
      f.write("select .not. type D* end -\n")
      f.write("select .not. type D* end\n\n")

      f.write("DrudeHardWall L_WALL 0.2 ! set the hard wall at 0.2 Angstrom\n\n")

      f.write("MINI SD nstep 5000 nprint 50 -\n")
      f.write("atom vatom vdistance cdie eps 1.0 vswitch fshift wmin 1.2 -\n")
      f.write("cutnb @cutoff ctofnb @ctofnb ctonnb @ctonnb -\n")
      f.write("e14fac 0.0\n\n") 
   
      f.write("! calculate molecular polarizability\n")
      f.write("VIBRAN\n")
      f.write("DIAGONALIZE\n")
      f.write("PRINT NORMAL  DIPOLES SELECT all END\n")
      f.write("END\n")
   with open('polar.sh','w') as f:
      f.write("#!/bin/bash\n")
      f.write("charmm < polar.inp > polar.log\n")
   sublist = ('sh polar.sh').split()
   output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
   os.remove('polar.sh')
   
   return
   


def prepare_analysis(c,temp,xyz,nmol,N_end,N_start=1,time=100,polar_flag=False):

   with open('write_pdb.inp','w') as f:
      f.write("* setup write out combined pdb\n") 
      f.write("*\n\n") 

      f.write("! this is necessart to work with non-interger total charge\n")
      f.write("bomlev -1\n\n")

      f.write("read rtf card name ../lin.rtf\n")
      f.write("read para card name ../lin.prm\n\n")

      f.write("set residue = LIN\n")
      f.write("READ SEQUENCE LIN {}\n".format(nmol))
      if polar_flag:
         f.write("generate @residue first none last none setup warn drude dmass 0.4 !show\n\n")
      else:
         f.write("generate @residue first none last none setup warn  !show\n\n")

      f.write("scalar charge show\n\n")

      f.write("open write unit 10 card name traj.pdb\n\n")

      f.write("set starttrj = {}\n".format(N_start))
      f.write("set endtrj = {}\n".format(N_end))
      f.write("calc ntrj = @endtrj - @starttrj +1\n\n")

      f.write("! open trj files\n")
      f.write("set i = 0\n")
      f.write("label openloop\n")
      f.write("   calc openu = @i + 101\n")
      f.write("   calc trj = @starttrj + @i\n")
      f.write("   open read unit @openu file name ../@trj.trj \n")
      f.write("   incr i by 1\n")
      f.write("if @i LT @ntrj goto openloop\n\n")

      f.write("traj first 101  nunit @ntrj\n")
      f.write("set frame 0\n")
      f.write("calc ntot = {}*@ntrj\n".format(time))
      f.write("label loop\n")
      f.write("   traj read\n")
      f.write("   write coor pdb unit 10\n")
      f.write("   incr frame by 1\n")
      f.write("if @frame LT @ntot goto loop\n")


   with open('write_vol.inp','w') as f:
      f.write("* write out boxlength vs time\n") 
      f.write("*\n\n") 

      f.write("! this is necessart to work with non-interger total charge\n")
      f.write("bomlev -1\n\n")

      f.write("read rtf card name ../lin.rtf\n")
      f.write("read para card name ../lin.prm\n\n")

      f.write("set residue = LIN\n")
      f.write("READ SEQUENCE LIN {}\n".format(nmol))
      if polar_flag:
         f.write("generate @residue first none last none setup warn drude dmass 0.4 !show\n\n")
      else:
         f.write("generate @residue first none last none setup warn  !show\n\n")

      f.write("scalar charge show\n\n")

      f.write("set starttrj = {}\n".format(N_start))
      f.write("set endtrj = {}\n".format(N_end))
      f.write("calc ntrj = @endtrj - @starttrj +1\n\n")

      f.write("! open trj files\n")
      f.write("set i = 0\n")
      f.write("label openloop\n")
      f.write("   calc openu = @i + 101\n")
      f.write("   calc trj = @starttrj + @i\n")
      f.write("   open read unit @openu file name ../@trj.trj \n")
      f.write("   incr i by 1\n")
      f.write("if @i LT @ntrj goto openloop\n\n")

      f.write("open write unit 20 card name volume.txt\n")
      f.write("CORREL maxt 20000 \n")
      f.write("ENTER LXCD CELL A\n")
      f.write("traj first 101  nunit @ntrj\n")
      f.write("WRITE LXCD UNIT 20 DUMB TIME\n")
      f.write("* title\n")

   sublist = ('../../../../{} -o single.pdb'.format(xyz)).split()
   xyz_to_pdb.main(sublist)
   with open('polar.inp','w') as f: 
      f.write("* setup structures\n")
      f.write("*\n\n") 

      f.write("! this is necessart to work with non-interger total charge\n")
      f.write("bomlev -1\n")
      f.write("set temp = {}\n".format(temp))
      f.write("read rtf card name ../lin.rtf\n")
      f.write("read para card name ../lin.prm\n")
      f.write("open unit 91 card read name single.pdb\n")
      f.write("read sequ pdb unit 91\n")
      f.write("set residue LIN\n\n")

      f.write("generate @residue first none last none setup warn drude dmass 0.4 !show\n\n")

      f.write("rewind unit 91\n")
      f.write("read coor pdb unit 91\n")
      f.write("coor sdrude\n")
      f.write("coor shake\n")
      f.write("coor print\n\n")

      f.write("! Dimension of the solvent volume\n")
      f.write("set 7 = 9999\n")
      f.write("set 8 @7\n")
      f.write("set 9 @8\n\n")

      f.write("set cutoff 999\n")
      f.write("set ctofnb = @cutoff\n")
      f.write("Calc ctofnb = @ctofnb - 2.0\n")
      f.write("set ctonnb = @ctofnb\n")
      f.write("Calc ctonnb = @ctonnb - 1.0\n\n")

      f.write("! shake constraint\n")
      f.write("SHAKE bonh param nofast -\n")
      f.write("select .not. type D* end -\n")
      f.write("select .not. type D* end\n\n")

      f.write("DrudeHardWall L_WALL 0.2 ! set the hard wall at 0.2 Angstrom\n\n")

      f.write("MINI SD nstep 5000 nprint 50 -\n")
      f.write("atom vatom vdistance cdie eps 1.0 vswitch fshift wmin 1.2 -\n")
      f.write("cutnb @cutoff ctofnb @ctofnb ctonnb @ctonnb -\n")
      f.write("e14fac 0.0\n\n") 
   
      f.write("! calculate molecular polarizability\n")
      f.write("VIBRAN\n")
      f.write("DIAGONALIZE\n")
      f.write("PRINT NORMAL  DIPOLES SELECT excluded END\n")
      f.write("END\n")



   with open('exe.sh','w') as f:
      f.write("#!/bin/bash\n")
      f.write("charmm < write_pdb.inp > write_pdb.log\n")
      f.write("charmm < write_vol.inp > write_vol.log\n")
      f.write("charmm < polar.inp > polar.log\n")
      f.write("{}/FF_functions/add_time_pdb.py\n".format(c["taffi_path"]))
      f.write("{}/FF_functions/read_charmm_energy.py -log \'../md.*.log\' -N_start 2 --silent\n".format(c["taffi_path"]))
      f.write("wait\n")

   """

   if os.path.isdir('eps') is False:
      os.makedirs('eps')
   os.chdir('eps')
   if polar_flag:
      sublist = ('../../../../../{}  -N 1 -o FF -q {} -gens {} -mixing_rule wh --polar -FF '.format(xyz,c['charge'],c['gens'])).split()
   else:
      sublist = ('../../../../../{}  -N 1 -o FF -q {} -gens {} -mixing_rule wh -FF '.format(xyz,c['charge'],c['gens'])).split()
   sublist.append('{}'.format(c['fom_ff']))
   gen_itp_newdrude.main(sublist)
   with open('topol.top','w') as f:
      f.write("#include \"FF/FF.itp\"\n\n") 

      f.write("[ system ]\n")
      f.write("Liquid\n\n")

      f.write("[ molecules ]\n")
      f.write("LIN {}\n".format(nmol))
   shutil.copy('../../initial.pdb','./initial.pdb')
   if polar_flag:
      add_drude_pdb.main(['initial.pdb'])
   with open('em.mdp','w') as f:
      f.write("; RUN CONTROL PARAMETERS\n")
      f.write("integrator = steep\n")
      f.write("; Start time and timestep in ps\n")
      f.write("nsteps = -1\n")
      f.write("emtol = 1000\n")
      f.write("; OUTPUT CONTROL OPTIONS\n")
      f.write("nstxout = 1\n")
      f.write("; Output frequency for energies to log file and energy file\n")
      f.write("nstlog = 100\n")
      f.write("nstcalcenergy = 1\n")
      f.write("nstenergy = 1\n")
      f.write("; NEIGHBORSEARCHING PARAMETERS\n")
      f.write("nstlist = 1\n")
      f.write("ns-type = Grid\n")
      f.write("pbc = xyz\n\n")
      if polar_flag:
         f.write("; DRUDE PARAMETERS\n") 
         f.write("drude = yes\n")
         f.write("drude-mode = SCF\n")
   if polar_flag:
      command = "/home/lin1209/install/gromacs_2/bin/gmx grompp  -f em.mdp  -c initial_drude.pdb -o eps.tpr"
   else:
      command = "/home/lin1209/install/gromacs_2/bin/gmx grompp  -f em.mdp  -c initial.pdb -o eps.tpr"
   output = subprocess.Popen(command.split(),stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[-1]
   os.chdir('..')
   """
         
   return


def prepare_cond(nmol,xyz,Vmol,Temp,run,polar_flag=False,shake_flag=False):

   shutil.copy('../../FF/lin.prm','./lin.prm')
   shutil.copy('../../FF/lin.rtf','./lin.rtf')
   if os.path.isfile('initial.pdb') is False:
      sublist = ('{} -xyz ../../../{} -Vmol {} --use_Vmol'.format(nmol,xyz,Vmol)).split() 
      #sublist = ('{} -xyz ../../../{} -dens {} '.format(nmol,xyz,Vmol)).split() 
      box = density_converter.main(sublist)
   write_setup(polar_flag)
   with open('setup.sh','w') as f:
      f.write("#!/bin/bash\n")
      f.write("charmm < setup.inp > setup.log\n")
   sublist = ('sh setup.sh').split()
   output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
   os.remove('setup.sh')
   
   random.seed()
   rseed = random.randint(100000,999999)
   #box = (Vmol*nmol)**(1/3) + 5 #XXX problematic when input is density instead of vmol
   write_nvt(Temp,nmol,box,rseed,polar_flag,shake_flag)
   write_equil(Temp,nmol,box,rseed,polar_flag,shake_flag)
   write_extend(Temp,nmol,box,rseed,polar_flag,shake_flag)
   write_sh(Temp,run,global_c['fom_md_cond_q'],global_c['fom_md_cond_wt'])

   return

def write_sh(Temp,run,queue,walltime):
   with open('run.sh','w') as f:
      f.write("#!/bin/bash\n\n") 

      f.write("#SBATCH --job-name ext_{}K_{}.1\n".format(Temp,run))
      f.write("#SBATCH -N 1\n")
      f.write("#SBATCH -n 1\n")
      f.write("#SBATCH -t {}:00:00\n".format(walltime))
      f.write("#SBATCH -A {}\n".format(queue))
      f.write("#SBATCH -o charmm.out\n")
      f.write("#SBATCH -e charmm.err\n\n")
                                                                                                                                                                                                            
      # XXX hard coded, need to modified
      f.write("export PATH=/apps/cent7/xalt/bin:/apps/cent7/vmd/vmd-1.9.3/bin:/usr/lib64/qt-3.3/bin:/apps/cent7/intel/impi/2017.1.132/bin64:/apps/cent7/intel/itac/2017.1.024/bin:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/bin/intel64:/apps/cent7/intel/inspector_2017.1.1.484836/bin64:/apps/cent7/intel/advisor_2017.1.1.486553/bin64:/apps/cent7/intel/vtune_amplifier_xe_2017/bin64:/apps/cent7/intel/opencl-1.2-sdk-6.3.0.1904/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/hpss/bin:/opt/hsi/bin:/opt/ibutils/bin:/opt/thinlinc/bin:/bin:/opt/hpss/bin:/opt/hsi/bin:/opt/intel/mic/bin:/sbin:/usr/bin:/usr/lib64/alliance/bin:/usr/lib64/ccache:/usr/lib64/qt-3.3/bin:/usr/libexec/lightdm:/usr/libexec/sdcc:/usr/local/bin:/usr/local/sbin:/usr/lpp/mmfs/bin:/usr/sbin:/usr/site/rcac/scripts:/usr/site/rcac/sbin:/usr/site/rcac/bin:/usr/site/rcac/scripts:/usr/site/rcac/support_scripts:/home/lin1209/bin:/depot/bsavoie/apps/openbabel/bin:/depot/bsavoie/apps/charmm/exec/gnu:/home/lin1209/install/gnome/bin:/home/lin1209/install/charmm/exec/gnu\n\n")

      f.write("export LD_LIBRARY_PATH=/apps/cent7/vmd/vmd-1.9.3/lib:/apps/cent7/intel/impi/2017.1.132/lib64:/apps/cent7/intel/itac/2017.1.024/intel64/slib:/apps/cent7/intel/itac/2017.1.024/lib:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/ipp/lib/intel64:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/mkl/lib/intel64:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/tbb/lib/intel64/gcc4.7:/apps/cent7/intel/opencl-1.2-sdk-6.3.0.1904/lib64:/apps/cent7/intel/debugger_2017/libipt/intel64/lib:/apps/cent7/intel/debugger_2017/iga/lib\n\n")

      f.write("simnum=1\n")
      f.write("rm ${simnum}.*\n")
      f.write("charmm < md-run.inp >> md.${simnum}.log\n")

   subprocess.call("chmod 777 run.sh", shell=True)

   with open('equil.sh','w') as f:
      f.write("#!/bin/bash\n\n") 

      f.write("#SBATCH --job-name ext_{}K_{}.equil\n".format(Temp,run))
      f.write("#SBATCH -N 1\n")
      f.write("#SBATCH -n 1\n")
      f.write("#SBATCH -t {}:00:00\n".format(walltime))
      f.write("#SBATCH -A {}\n".format(queue))
      f.write("#SBATCH -o equil.out\n")
      f.write("#SBATCH -e equil.err\n\n")
                                                                                                                                                                                                            
      f.write("export PATH=/apps/cent7/xalt/bin:/apps/cent7/vmd/vmd-1.9.3/bin:/usr/lib64/qt-3.3/bin:/apps/cent7/intel/impi/2017.1.132/bin64:/apps/cent7/intel/itac/2017.1.024/bin:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/bin/intel64:/apps/cent7/intel/inspector_2017.1.1.484836/bin64:/apps/cent7/intel/advisor_2017.1.1.486553/bin64:/apps/cent7/intel/vtune_amplifier_xe_2017/bin64:/apps/cent7/intel/opencl-1.2-sdk-6.3.0.1904/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/hpss/bin:/opt/hsi/bin:/opt/ibutils/bin:/opt/thinlinc/bin:/bin:/opt/hpss/bin:/opt/hsi/bin:/opt/intel/mic/bin:/sbin:/usr/bin:/usr/lib64/alliance/bin:/usr/lib64/ccache:/usr/lib64/qt-3.3/bin:/usr/libexec/lightdm:/usr/libexec/sdcc:/usr/local/bin:/usr/local/sbin:/usr/lpp/mmfs/bin:/usr/sbin:/usr/site/rcac/scripts:/usr/site/rcac/sbin:/usr/site/rcac/bin:/usr/site/rcac/scripts:/usr/site/rcac/support_scripts:/home/lin1209/bin:/depot/bsavoie/apps/openbabel/bin:/depot/bsavoie/apps/charmm/exec/gnu:/home/lin1209/install/gnome/bin:/home/lin1209/install/charmm/exec/gnu\n\n")

      f.write("export LD_LIBRARY_PATH=/apps/cent7/vmd/vmd-1.9.3/lib:/apps/cent7/intel/impi/2017.1.132/lib64:/apps/cent7/intel/itac/2017.1.024/intel64/slib:/apps/cent7/intel/itac/2017.1.024/lib:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/ipp/lib/intel64:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/mkl/lib/intel64:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/tbb/lib/intel64/gcc4.7:/apps/cent7/intel/opencl-1.2-sdk-6.3.0.1904/lib64:/apps/cent7/intel/debugger_2017/libipt/intel64/lib:/apps/cent7/intel/debugger_2017/iga/lib\n\n")
   
      f.write("charmm < md-nvt.inp >> nvt.log\n")
      f.write("date\n")
      f.write("charmm < md-equil.inp >> md.0.log\n")
      f.write("date\n")

   subprocess.call("chmod 777 equil.sh", shell=True)

   with open('extend.sh','w') as f:
      f.write("#!/bin/bash\n\n")

      f.write("presim=`ls *.sys.crd -lv | tail -n 1 | awk '{print $9}' | awk -F '.' '{printf(\"%d\",$1)}'`\n")
      f.write("if [ -z $presim ]\n") 
      f.write("then\n")
      f.write("   tmp=( $( sbatch equil.sh ) )\n")
      f.write("   ids_sub+=(${tmp[${#tmp[@]}-1]})\n")
      f.write("   /depot/bsavoie/data/Lin/taffi_beta/FF_functions/monitor_jobs.py ${ids_sub[@]}\n")
      f.write("fi\n\n")

      f.write("title=\"ext_{}K_{}\"\n".format(Temp,run))
      f.write("totsim={}\n".format(int(int(global_c['fom_t_tot'])/int(global_c['fom_t_ext']))))
      f.write("for i in $( seq 14 ); do\n")
      f.write("   ids_sub=()   \n") 
      f.write("   sim=`ls *.sys.crd -lv | tail -n 1 | awk '{print $9}' | awk -F '.' '{printf(\"%d\",$1+1)}'`\n")
      f.write("   presim=`ls *.sys.crd -lv | tail -n 1 | awk '{print $9}' | awk -F '.' '{printf(\"%d\",$1)}'`\n")
      f.write("   if [ -z $sim  ]\n")
      f.write("   then\n")
      f.write("      here=`pwd`\n")
      f.write("      echo \"$here failed\"\n")
      f.write("      exit\n")
      f.write("   fi\n")
      f.write("   box=`grep 'DYNA A ' md.${presim}.log | tail -n 1 | awk '{printf(\"%3.5f\",$4)}'` \n")
      f.write("   line=`grep 'set boxval' md-run.inp -n | awk -F ':' '{printf(\"%d\",$1)}'`\n")
      f.write("   sed -i \"${line}c set boxval = ${box}\" md-run.inp\n")
      f.write("   line=`grep 'set simnum' md-run.inp -n | awk -F ':' '{printf(\"%d\",$1)}'`\n")
      f.write("   sed -i \"${line}c set simnum = $sim\" md-run.inp\n")
      f.write("   line=`grep 'simnum=' run.sh -n | awk -F ':' '{printf(\"%d\",$1)}'`\n")
      f.write("   sed -i \"${line}c simnum=$sim\" run.sh\n")
      f.write("   line=`grep $title run.sh -n | awk -F ':' '{printf(\"%d\",$1)}'`\n")
      f.write("   sed -i \"${line}c #SBATCH --job-name ${title}.${sim}\" run.sh\n")
      f.write("   if [ $presim -ge $totsim ] \n")
      f.write("   then \n")
      f.write("      echo \"total simnum: $presim\" \n")
      f.write("      exit \n") 
      f.write("   fi\n")
      f.write("   tmp=( $( sbatch run.sh ) )\n")
      f.write("   ids_sub+=(${tmp[${#tmp[@]}-1]})\n")
      f.write("   /depot/bsavoie/data/Lin/taffi_beta/FF_functions/monitor_jobs.py ${ids_sub[@]}\n")
      f.write("done     \n")
   subprocess.call("chmod 777 extend.sh", shell=True)

   return



def write_nvt(Temp,nmol,box,rseed,polar_flag=False,shake_flag=False):
   with open('md-nvt.inp','w') as f:
      f.write("! this is necessart to work with non-interger total charge\n")
      f.write("bomlev -1\n\n")

      f.write("set simnum = 0\n")
      f.write("set temp = {}\n\n".format(Temp))

      f.write("! read the psf and coordinate file\n")
      f.write("read rtf card name lin.rtf\n")
      f.write("read para card name lin.prm\n\n")

      f.write("set residue = LIN\n")
      f.write("READ SEQUENCE LIN {}\n".format(nmol))
      #f.write("generate @residue first none last none setup warn drude dmass 0.4 !show\n")
      f.write("! noangle nodihedral is to prevent charmm from auto generate angle and dihedral base on bonds,\n")
      f.write("! those explicitly defined in rtf file will still be generated\n")
      f.write("generate @residue first none last none setup noangle nodihedral warn drude dmass 0.4 !show\n")
      f.write("OPEN READ CARD UNIT 10 NAME lin.crd\n")
      f.write("read coor card name lin.crd\n\n")


      if polar_flag:
         f.write("coor sdrude\n")

      if shake_flag:
         f.write("coor shake\n\n")
         f.write("! shake constraint\n")
         f.write("SHAKE bonh param nofast -\n")
         f.write("select .not. type D* end -\n")
         f.write("select .not. type D* end\n\n")
      
      if polar_flag:
         f.write("DrudeHardWall L_WALL 0.2 ! set the hard wall at 0.2 Angstrom\n\n")

      f.write("! SETUP CRYSTAL (DEFINE, BUILD),\n")
      f.write("set boxval = {:6.3f}\n".format(box))
      f.write("crystal define cubic @boxval @boxval @boxval 90. 90. 90\n")
      f.write("crystal build noper 0 cutoff 16\n\n")

      f.write("! Set up images -- center the protein by segment and the solvent by residue\n")
      f.write("image byres xcen 0.0 ycen 0.0 zcen 0.0 sele all end\n\n")

      f.write("! set up nonbond parameters\n")
      f.write("nbond nbxmod 5 VDW vatom vswitch e14fac 0.0 cutnb 14.0 cutim 14 ctofnb 12.0 ctonnb 10.0 -\n")
      f.write("inbfrq -1 imfrq -1 wmin 1.0 -\n")
      f.write("ELEC atom cdiel EWALD PMEWALD KAPPA 0.34 FFTX 32 FFTY 32 FFTZ 32 ORDER 6\n")
      f.write("energy\n\n")

      if polar_flag:
         f.write("tpcontrol nther 2 nstep 20 -\n")
         f.write("ther 1 tau 0.1 tref @temp select .not. type d* end -\n")
         f.write("ther 2 tau 0.005 tref 1 select type d* end -\n\n")
      else:
         f.write("tpcontrol nther 1 nstep 1 -\n")
         f.write("ther 1 tau 0.1 tref @temp select all  end \n")

      f.write("set NSTEPS = 20000\n")
      f.write("set printfreq = 1000\n")
      f.write("OPEN WRITE CARD UNIT 62 NAME nvt.nose ! info for temperaature\n")
      f.write("open unit 35 write form name nvt.rst ! restart file\n")
      f.write("open unit 32 write unfo name nvt.trj ! trajectory\n")
      f.write("open unit 33 write form name nvt.kunit ! Temperature\n")
      f.write("open unit 34 write unfo name nvt.vel ! velocity\n\n")

      f.write("DYNA vv2 start -\n")
      f.write("nstep @NSTEPS time 0.001 -\n")
      f.write("iunread -1 -\n")
      f.write("iunwrite 35 iuncrd 32 iunvel 34 kunit 33 -\n")
      f.write("nsavc @printfreq nsavv @printfreq nprint @printfreq isvfrq @printfreq -\n")
      f.write("inbgrq 1 imgfrq 1 - \n")
      f.write("iprfrq @NSTEPS ntrfrq @NSTEPS - ! what's used for swm4 water\n")
      f.write("firstt @temp -\n")
      f.write("finalt @temp -\n")
      f.write("iasvel 1 iseed {:7d} -\n".format(rseed))
      f.write("iuno 62 nsnos @printfreq\n\n")

      f.write("open write unit 10 card name nvt.crd\n")
      f.write("write coor card unit 10\n")
      f.write("close unit 10\n\n")

      f.write("stop\n")
   return

def write_equil(Temp,nmol,box,rseed,polar_flag=False,shake_flag=False):
   with open('md-equil.inp','w') as f:
      f.write("! this is necessart to work with non-interger total charge\n")
      f.write("bomlev -1\n\n")

      f.write("set simnum = 0\n")
      f.write("set temp = {}\n\n".format(Temp))

      f.write("! read the psf and coordinate file\n")
      f.write("read rtf card name lin.rtf\n")
      f.write("read para card name lin.prm\n\n")

      f.write("set residue = LIN\n")
      f.write("READ SEQUENCE LIN {}\n".format(nmol))
      #f.write("generate @residue first none last none setup warn drude dmass 0.4 !show\n")
      f.write("! noangle nodihedral is to prevent charmm from auto generate angle and dihedral base on bonds,\n")
      f.write("! those explicitly defined in rtf file will still be generated\n")
      f.write("generate @residue first none last none setup noangle nodihedral warn drude dmass 0.4 !show\n")
      f.write("OPEN READ CARD UNIT 10 NAME nvt.crd\n")
      f.write("read coor card name nvt.crd\n\n")

      if polar_flag:
         f.write("coor sdrude\n")
   
      if shake_flag:
         f.write("coor shake\n\n")

         f.write("! shake constraint\n")
         f.write("SHAKE bonh param nofast -\n")
         f.write("select .not. type D* end -\n")
         f.write("select .not. type D* end\n\n")

      if polar_flag:
         f.write("DrudeHardWall L_WALL 0.2 ! set the hard wall at 0.2 Angstrom\n\n")

      f.write("! SETUP CRYSTAL (DEFINE, BUILD),\n")
      f.write("set boxval = {:6.3f}\n".format(box))
      f.write("crystal define cubic @boxval @boxval @boxval 90. 90. 90\n")
      f.write("crystal build noper 0 cutoff 16\n\n")

      f.write("! Set up images -- center the protein by segment and the solvent by residue\n")
      f.write("image byres xcen 0.0 ycen 0.0 zcen 0.0 sele all end\n\n")

      f.write("! set up nonbond parameters\n")
      f.write("nbond nbxmod 5 VDW vatom vswitch e14fac 0.0 cutnb 14.0 cutim 14 ctofnb 12.0 ctonnb 10.0 -\n")
      f.write("inbfrq 1 imfrq 1 wmin 1.0 -\n")
      f.write("ELEC atom cdiel EWALD PMEWALD KAPPA 0.34 FFTX 32 FFTY 32 FFTZ 32 ORDER 6\n")
      f.write("energy\n\n")

      if polar_flag:
         f.write("tpcontrol nther 2 nstep 20 -\n")
         f.write("ther 1 tau 0.1 tref @temp select .not. type d* end -\n")
         f.write("ther 2 tau 0.005 tref 1 select type d* end -\n")
         f.write("baro btau 0.1 pref 1.00\n\n")
      else:
         f.write("tpcontrol nther 1 nstep 1 -\n")
         f.write("ther 1 tau 0.1 tref @temp select all end -\n")
         f.write("baro btau 0.1 pref 1.00\n\n")

      f.write("set NSTEPS = 20000\n")
      f.write("set printfreq = 1000\n")
      f.write("open read unit 91 card name nvt.rst\n")
      f.write("OPEN WRITE CARD UNIT 62 NAME @simnum.nose ! info for temperaature\n")
      f.write("open unit 35 write form name @simnum.rst ! restart file\n")
      f.write("open unit 32 write unfo name @simnum.trj ! trajectory\n")
      f.write("open unit 33 write form name @simnum.kunit ! Temperature\n")
      f.write("open unit 34 write unfo name @simnum.vel ! velocity\n\n")

      f.write("DYNA vv2 start -\n")
      f.write("nstep @NSTEPS time 0.001 -\n")
      f.write("iunread 91 -\n")
      f.write("iunwrite 35 iuncrd 32 iunvel 34 kunit 33 -\n")
      f.write("nsavc @printfreq nsavv @printfreq nprint @printfreq isvfrq @printfreq -\n")
      f.write("inbfrq -1 imgfrq -1 - \n")
      f.write("iprfrq @NSTEPS ntrfrq @NSTEPS - ! what's used for swm4 water\n")
      f.write("firstt @temp -\n")
      f.write("finalt @temp -\n")
      f.write("iasvel 1 iseed {:7d} -\n".format(rseed))
      f.write("iuno 62 nsnos @printfreq\n\n")

      f.write("open write unit 10 card name @simnum.sys.crd\n")
      f.write("write coor card unit 10\n")
      f.write("close unit 10\n\n")

      f.write("stop\n")
   return

def write_extend(Temp,nmol,box,rseed,polar_flag=False,shake_flag=False):
   with open('md-run.inp','w') as f:
      f.write("! this is necessart to work with non-interger total charge\n")
      f.write("bomlev -1\n\n")

      f.write("set simnum = 1\n")
      f.write("calc prevsim = @simnum -1\n")
      f.write("set temp = {}\n\n".format(Temp))

      f.write("! read the psf and coordinate file\n")
      f.write("read rtf card name lin.rtf\n")
      f.write("read para card name lin.prm\n\n")

      f.write("set residue = LIN\n")
      f.write("READ SEQUENCE LIN {}\n".format(nmol))
      #f.write("generate @residue first none last none setup warn drude dmass 0.4 !show\n")
      f.write("! noangle nodihedral is to prevent charmm from auto generate angle and dihedral base on bonds,\n")
      f.write("! those explicitly defined in rtf file will still be generated\n")
      f.write("generate @residue first none last none setup noangle nodihedral warn drude dmass 0.4 !show\n")
      f.write("OPEN READ CARD UNIT 10 NAME @prevsim.sys.crd\n")
      f.write("read coor card name @prevsim.sys.crd\n\n")


      if polar_flag:
         f.write("coor sdrude\n")

      if shake_flag:
         f.write("coor shake\n\n")

         f.write("! shake constraint\n")
         f.write("SHAKE bonh param nofast -\n")
         f.write("select .not. type D* end -\n")
         f.write("select .not. type D* end\n\n")

      if polar_flag:
         f.write("DrudeHardWall L_WALL 0.2 ! set the hard wall at 0.2 Angstrom\n\n")

      f.write("! SETUP CRYSTAL (DEFINE, BUILD),\n")
      f.write("set boxval = {:6.3f}\n".format(box))
      f.write("crystal define cubic @boxval @boxval @boxval 90. 90. 90\n")
      f.write("crystal build noper 0 cutoff 16\n\n")

      f.write("! Set up images -- center the protein by segment and the solvent by residue\n")
      f.write("image byres xcen 0.0 ycen 0.0 zcen 0.0 sele all end\n\n")

      f.write("! set up nonbond parameters\n")
      f.write("nbond nbxmod 5 VDW vatom vswitch e14fac 0.0 cutnb 14.0 cutim 14 ctofnb 12.0 ctonnb 10.0 -\n")
      f.write("inbfrq -1 imfrq -1 wmin 1.0 -\n")
      f.write("ELEC atom cdiel EWALD PMEWALD KAPPA 0.34 FFTX 32 FFTY 32 FFTZ 32 ORDER 6\n")
      f.write("energy\n\n")

      if polar_flag:
         f.write("tpcontrol nther 2 nstep 20 -\n")
         f.write("ther 1 tau 0.1 tref @temp select .not. type d* end -\n")
         f.write("ther 2 tau 0.005 tref 1 select type d* end -\n")
         f.write("baro btau 0.1 pref 1.00\n\n")
      else:
         f.write("tpcontrol nther 1 nstep 1 -\n")
         f.write("ther 1 tau 0.1 tref @temp select all  end -\n")
         f.write("baro btau 0.1 pref 1.00\n\n")

      f.write("set NSTEPS = {}\n".format(int(global_c['fom_t_ext'])*1000))
      f.write("set printfreq = 1000\n")
      f.write("open read unit 91 card name @prevsim.rst\n")
      f.write("OPEN WRITE CARD UNIT 62 NAME @simnum.nose ! info for temperaature\n")
      f.write("open unit 35 write form name @simnum.rst ! restart file\n")
      f.write("open unit 32 write unfo name @simnum.trj ! trajectory\n")
      f.write("open unit 33 write form name @simnum.kunit ! Temperature\n")
      f.write("open unit 34 write unfo name @simnum.vel ! velocity\n\n")

      f.write("DYNA vv2 rest -\n")
      f.write("nstep @NSTEPS time 0.001 -\n")
      f.write("iunread 91 -\n")
      f.write("iunwrite 35 iuncrd 32 iunvel 34 kunit 33 -\n")
      f.write("nsavc @printfreq nsavv @printfreq nprint @printfreq isvfrq @printfreq -\n")
      f.write("inbfrq -1 imgfrq -1 -\n")
      f.write("iprfrq @NSTEPS ntrfrq @NSTEPS - ! what's used for swm4 water\n")
      f.write("firstt @temp -\n")
      f.write("finalt @temp -\n")
      f.write("iasvel 1 iseed {:7d} -\n".format(rseed))
      f.write("iuno 62 nsnos @printfreq\n\n")

      f.write("open write unit 10 card name @simnum.sys.crd\n")
      f.write("write coor card unit 10\n")
      f.write("close unit 10\n\n")

      f.write("stop\n")
   return


def write_setup(polar_flag=False):
   with open('setup.inp','w') as f:
      f.write("* setup structures\n")
      f.write("*\n\n")

      f.write("! this is necessart to work with non-interger total charge\n")
      f.write("bomlev -1\n\n")

      f.write("read rtf card name lin.rtf\n")
      f.write("read para card name lin.prm\n")
      f.write("open unit 91 card read name initial.pdb\n")
      f.write("read sequ pdb unit 91\n")
      f.write("set residue LIN\n\n")

      #f.write("generate @residue first none last none setup warn drude dmass 0.4 !show\n")
      f.write("! noangle nodihedral is to prevent charmm from auto generate angle and dihedral base on bonds,\n")
      f.write("! those explicitly defined in rtf file will still be generated\n")
      f.write("generate @residue first none last none setup noangle nodihedral warn drude dmass 0.4 !show\n")

      f.write("rewind unit 91\n")
      f.write("read coor pdb unit 91\n")

      if polar_flag:
         f.write("coor sdrude\n") 
      f.write("coor shake\n")
      f.write("coor print\n\n") 

      f.write("open write unit 10 card name lin.crd\n")
      f.write("write coor card unit 10\n\n")

      f.write("stop\n")
   return

def prepare_gas(xyz,Temp,run,polar_flag=False,shake_flag=False):

   shutil.copy('../../FF/lin.prm','./lin.prm')
   shutil.copy('../../FF/lin.rtf','./lin.rtf')
   sublist = ('../../../{} -o single.pdb'.format(xyz)).split()
   xyz_to_pdb.main(sublist)
   random.seed()
   rseed = random.randint(100000,999999)
   write_gas_inp(Temp,rseed,polar_flag,shake_flag)
   write_gas_sh(Temp,run,global_c['fom_md_cond_q'],global_c['fom_md_cond_wt'])

   return

def write_gas_sh(Temp,run,queue,walltime):
   with open('gas.sh','w') as f:
      f.write("#!/bin/bash\n\n") 

      f.write("#SBATCH --job-name ext_{}K_{}.gas\n".format(Temp,run))
      f.write("#SBATCH -N 1\n")
      f.write("#SBATCH -n 1\n")
      f.write("#SBATCH -t {}:00:00\n".format(walltime))
      f.write("#SBATCH -A {}\n".format(queue))
      f.write("#SBATCH -o gas.out\n")
      f.write("#SBATCH -e gas.err\n\n")
                                                                            

      # XXX: hard coded, need to modified                                                                                                                                
      f.write("export PATH=/apps/cent7/xalt/bin:/apps/cent7/vmd/vmd-1.9.3/bin:/usr/lib64/qt-3.3/bin:/apps/cent7/intel/impi/2017.1.132/bin64:/apps/cent7/intel/itac/2017.1.024/bin:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/bin/intel64:/apps/cent7/intel/inspector_2017.1.1.484836/bin64:/apps/cent7/intel/advisor_2017.1.1.486553/bin64:/apps/cent7/intel/vtune_amplifier_xe_2017/bin64:/apps/cent7/intel/opencl-1.2-sdk-6.3.0.1904/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/hpss/bin:/opt/hsi/bin:/opt/ibutils/bin:/opt/thinlinc/bin:/bin:/opt/hpss/bin:/opt/hsi/bin:/opt/intel/mic/bin:/sbin:/usr/bin:/usr/lib64/alliance/bin:/usr/lib64/ccache:/usr/lib64/qt-3.3/bin:/usr/libexec/lightdm:/usr/libexec/sdcc:/usr/local/bin:/usr/local/sbin:/usr/lpp/mmfs/bin:/usr/sbin:/usr/site/rcac/scripts:/usr/site/rcac/sbin:/usr/site/rcac/bin:/usr/site/rcac/scripts:/usr/site/rcac/support_scripts:/home/lin1209/bin:/depot/bsavoie/apps/openbabel/bin:/depot/bsavoie/apps/charmm/exec/gnu:/home/lin1209/install/gnome/bin:/home/lin1209/install/charmm/exec/gnu\n\n")

      f.write("export LD_LIBRARY_PATH=/apps/cent7/vmd/vmd-1.9.3/lib:/apps/cent7/intel/impi/2017.1.132/lib64:/apps/cent7/intel/itac/2017.1.024/intel64/slib:/apps/cent7/intel/itac/2017.1.024/lib:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/ipp/lib/intel64:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/mkl/lib/intel64:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/tbb/lib/intel64/gcc4.7:/apps/cent7/intel/opencl-1.2-sdk-6.3.0.1904/lib64:/apps/cent7/intel/debugger_2017/libipt/intel64/lib:/apps/cent7/intel/debugger_2017/iga/lib\n\n")
   
      f.write("charmm < md-gas.inp > gas.log\n")

   subprocess.call("chmod 777 gas.sh", shell=True)


def write_gas_inp(Temp,rseed,polar_flag=False,shake_flag=False):
   with open('md-gas.inp','w') as f: 
      f.write("* setup structures\n")
      f.write("*\n\n") 

      f.write("! this is necessart to work with non-interger total charge\n")
      f.write("bomlev -1\n")
      f.write("set temp = {}\n".format(Temp))
      f.write("read rtf card name lin.rtf\n")
      f.write("read para card name lin.prm\n")
      f.write("open unit 91 card read name single.pdb\n")
      f.write("read sequ pdb unit 91\n")
      f.write("set residue LIN\n\n")

      #f.write("generate @residue first none last none setup warn drude dmass 0.4 !show\n")
      f.write("! noangle nodihedral is to prevent charmm from auto generate angle and dihedral base on bonds,\n")
      f.write("! those explicitly defined in rtf file will still be generated\n")
      f.write("generate @residue first none last none setup noangle nodihedral warn drude dmass 0.4 !show\n")

      f.write("rewind unit 91\n")
      f.write("read coor pdb unit 91\n")
      f.write("coor sdrude\n")
      f.write("coor shake\n")
      f.write("coor print\n\n")

      f.write("! Dimension of the solvent volume\n")
      f.write("set 7 = 999\n")
      f.write("set 8 @7\n")
      f.write("set 9 @8\n\n")

      f.write("set cutoff 99\n")
      f.write("set ctofnb = @cutoff\n")
      f.write("Calc ctofnb = @ctofnb - 2.0\n")
      f.write("set ctonnb = @ctofnb\n")
      f.write("Calc ctonnb = @ctonnb - 1.0\n\n")

      if shake_flag:
         f.write("! shake constraint\n")
         f.write("SHAKE bonh param nofast -\n")
         f.write("select .not. type D* end -\n")
         f.write("select .not. type D* end\n\n")

      if polar_flag:
         f.write("DrudeHardWall L_WALL 0.2 ! set the hard wall at 0.2 Angstrom\n\n")

      f.write("MINI SD nstep 5000 nprint 50 -\n")
      f.write("atom vatom vdistance cdie eps 1.0 vswitch fshift wmin 1.2 -\n")
      f.write("cutnb @cutoff ctofnb @ctofnb ctonnb @ctonnb -\n")
      f.write("e14fac 0.0\n\n") 

      f.write("! set up nonbond parameters\n")
      f.write("nbond atom vatom vdistance cdie eps 1.0 vswitch fshift wmin 1.2 -\n")
      f.write("cutnb @cutoff ctofnb @ctofnb ctonnb @ctonnb -\n")
      f.write("e14fac 0.0\n") 
      f.write("energy\n\n\n")

      if polar_flag:
         f.write("TPCONTROL NTHER 2 NHGAM 5.0 NHGAMD 0.0 -\n") 
         f.write("THER 1  TREF @temp  LANG SELECT .NOT. TYPE D* END  -\n")
         f.write("THER 2  TREF   0.00  LANG SELECT TYPE D* END \n") 
         #f.write("tpcontrol nther 2 nstep 20 -\n")
         #f.write("ther 1 tau 0.1 tref @temp select .not. type d* end -\n")
         #f.write("ther 2 tau 0.005 tref 1 select type d* end -\n\n")
      else:
         f.write("tpcontrol nther 1 nstep 1 -\n")
         f.write("ther 1 tau 0.1 tref @temp select all end \n")

      f.write("set NSTEPS = 20000000\n")
      f.write("set printfreq = 1000\n")
      f.write("OPEN WRITE CARD UNIT 62 NAME gas.nose ! info for temperaature\n")
      f.write("open unit 35 write form name gas.rst ! restart file\n")
      f.write("open unit 32 write unfo name gas.trj ! trajectory\n")
      f.write("open unit 33 write form name gas.kunit ! Temperature\n")
      f.write("open unit 34 write unfo name gas.vel ! velocity\n\n")

      f.write("DYNA vv2 start -\n")
      f.write("nstep @NSTEPS time 0.001 -\n")
      f.write("iunread -1 -\n")
      f.write("iunwrite 35 iuncrd 32 iunvel 34 kunit 33 -\n")
      f.write("nsavc @printfreq nsavv @printfreq nprint @printfreq isvfrq @printfreq -\n")
      f.write("iprfrq @NSTEPS ntrfrq @NSTEPS - ! what's used for swm4 water\n")
      f.write("cutnb @cutoff ctofnb @ctofnb ctonnb @ctonnb -\n")
      f.write("firstt @temp -\n")
      f.write("finalt @temp -\n")
      f.write("iasvel 1 iseed {:7d} -\n".format(rseed))
      f.write("iuno 62 nsnos @printfreq\n\n")

      f.write("open write unit 10 card name gas.crd\n")
      f.write("write coor card unit 10\n")
      f.write("close unit 10\n\n")

      f.write("stop\n")

   return

def gas_submit(name,num,Temp,c):

    sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
    sublist.append('module load openmpi/3.1.4 \n/depot/bsavoie/apps/charmm/exec/gnu/charmm < ../md-gas.inp > gas.log')
    sublist= sublist + ('{}_gas_{}_{}K  -p 1 -t {} -q {} -ppn {}'.format(name,num,Temp,c["fom_md_gas_wt"],c["fom_md_gas_q"],c["fom_md_gas_ppn"])).split()
    output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    output = str(output,'utf-8')
    jobid = output.split("\n")[-2].split()[-1]

    return jobid


# Write lammps rerun.in.init file for rerunning the combined trajectory without averaging
def write_rerun(input_file,output_file):
    
    with open(input_file,'r') as f:
        with open(output_file,'w') as o:
            for lines in f:
                fields = lines.split()
                # reset run index, necessary for backwards compatibility, since lammps doesn't like to mix rerun with read_restart
                if len(fields) > 3 and fields[0] == "variable" and fields[1] == "run" and fields[2] == "index":
                    o.write("variable    run             index   0\n")
                elif len(fields) >= 2 and fields[0] == "read_restart":
                    o.write("read_data ${data_name}\n")
                elif len(fields) > 3 and fields[0] == "fix" and fields[1] == "averages" and fields[2] == "all" and fields[3] == "ave/time":
                    o.write("fix raw all ave/time 1 1 ${thermo_freq} v_my_temp v_my_rho v_my_vol v_my_pe v_my_edih v_my_evdw v_my_eel v_my_ent v_my_P file combined.thermo.rerun format %20.10g\n\n"+\
                            "# Rerun trajectory\n"+\
                            "rerun combined.lammpstrj dump x y z\n")
                    break
                else:
                    o.write(lines)
    return

# Description: Calculate the mobility based on the first continuous window of N datapoints with a loglog slope greater than a_thresh 
# 
# inputs: y_log (array): slope from the loglog derivative data
#         y_lin (array): slope from the linlin derivative data
#         a_thresh (float): slope used for evaluating the diffusive regime
#         N_thresh (int): number of continuous frames with y_log slope greater than a_thresh used for calculating the derivative.
#         start (int): index to start the search at
def calc_mobility(y_log,y_lin,a_thresh=0.95,N_thresh=10,start=0):

    count = 0
    inds = []
    for i in range(len(y_log)):
        if y_log[i] > a_thresh and i > start:
            count += 1
            inds += [i]
            if count >= N_thresh:
                return np.mean(y_lin[inds])/6.0*1.0E-16*1.0E9,np.std(y_lin[inds]/6.0*1.0E-16*1.0E9),inds
        else:
            count = 0
            inds = []
    return 0.0,0.0,0.0

# Description: Calculates the derivative of a dataset using linear least-squares regression on a running window of the data.
#
# inputs: x (array): x data for calculating the derivative
#         y (array): y data for calculating the derivative
#         abs_flag (boolean): option for taking the absolute value of the data before calculating the derivative
#         logx_opt (boolean): option for taking the log of the x data before calculating the derivative
#         logy_opt (boolean): option for taking the log of the y data before calculating the derivative
#         x_scale (float): scale factor applied to the x data before calculating the derivative
#         y_scale (float): scale factor applied to the y data before calculating the derivative
#         mode (string): forward or central windowing used for calculating the derivative
#         window (int): number of datapoints to use in the running window for calculating the derivative. In central mode, the 
#                       regression is applied to +/- this number of datapoints about the point that the derivative is being evaluated.
#                       In forward mode, + this number of datapoints about the point that the derivative is being evaluated.
def calc_derivative(x,y,abs_flag=False,logx_opt=False,logy_opt=False,x_scale=1.0,y_scale=1.0,mode="forward",window=10):

    valid_modes = ["forward","central"]
    if mode not in valid_modes:
        print("ERROR in calc_derivative: mode must be set to one of the following {}".format(valid_modes))
        quit()

    # Apply user-defined scaling options
    x_0 = deepcopy(x)
    y_0 = deepcopy(y)
    if abs_flag:
        x_0 = np.abs(x_0)
        y_0 = np.abs(y_0)
    if logx_opt:
        for ind, x_tmp in enumerate(x_0): # ADDED
            if x_tmp == 0:                # ADDED
                x_0[ind] = 1E-16          # ADDED
        x_0 = np.log10(x_0)
    if logy_opt:
        for ind, y_tmp in enumerate(y_0): # ADDED
            if y_tmp == 0:                # ADDED
                y_0[ind] = 1E-16          # ADDED
        y_0 = np.log10(y_0)
 
    # Apply scaling
    y_0 = y_0*y_scale
    x_0 = x_0*x_scale
    N = len(x_0)
    
    slopes = []
    stds = []

    # Calculate the derivative
    if mode == "central":
        for i in range(len(x_0)):

            # starting indices can only do forward windows
            if i < window:
                inds = list(range(0,i)) + list(range(i,i+window+1))

            # ending indices can only do backwards windows
            elif i+window >= N:
                inds = list(range(i-window,i)) + list(range(i,N))

            # the rest can sample the full window
            else:
                inds = list(range(i-window,i+window+1))

            slope, intercept, r_value, p_value, std_err = linregress(x_0[inds],y_0[inds])
            slopes += [slope]
            stds += [std_err]

    # Calculate the derivative
    if mode == "forward":
        for i in range(len(x_0)):

            # ending indices can only do fractional windows
            if i+window >= N:
                inds = list(range(i,N))

            # the rest can sample the full window (N-2 because derivative must involve at least 2 points
            elif i < N-2:
                inds = list(range(i,i+window+1))

            # calculate derivative
            if i < N-2:
                slope, intercept, r_value, p_value, std_err = linregress(x_0[inds],y_0[inds])

                slopes += [slope]
                stds += [std_err]

    return slopes,stds

# # Print FOM summary the densities and print summary of values
# kb=0.001987204118
# T=298
# kcal2kJ=4.184
# printf "\n%40s   %20s   %40s %20s" "Molecule" "density (g/cm^3)" "diff (1E5 * cm^2/s)" "H_vap (kcal/mol)" > summary.txt
# for i in *cond/; do
#     name=$( echo ${i} | awk -F '_cond' '{ print $1 }' )
#     density=$(${TAFFI_PATH}Parsers/avg_thermo.sh ${name}_cond/thermo.avg 1000 60000 3)
#     diff=$( awk '{ if ( $1 == "The" && $2 == "Self-Diffusion" && $3 == "Coefficient" ) { printf("%12.6f",$9*1.0E5) }}' ${name}_cond/diff/msd_parse.log)
#     U_vap=$( ${TAFFI_PATH}Parsers/avg_thermo.sh ${name}_gas/thermo.equil.avg 1000 0 5 )
#     U_bulk=$( ${TAFFI_PATH}Parsers/avg_thermo.sh ${name}_cond/thermo.avg 1000 0 5 )
#     N_mol=$( awk '{ if ( NR == 1 ) { print $1; exit } }' ${name}.xyz )
#     N_mol=$( echo "1000/${N_mol}" | bc )
#     U_bulk=$( echo "${U_bulk} / ${N_mol}" | bc -l )
#     H_vap=$( echo "(${U_vap} - ${U_bulk} + ${kb} * ${T} ) * ${kcal2kJ}" | bc -l )
#     printf "\n%40s %12.4f %40s %23.4f" ${name} ${density} ${diff} ${H_vap} >> summary.txt
# done
# echo "" >> summary.txt

    quit()

# Wrapper function for parsing the run value from a lammps file. Used to determine job completion.
def get_run_N(input_file):
    with open(input_file,'r') as f:
        for lines in f:
            fields = lines.split()
            if len(fields) >= 4 and fields[0] == "variable" and fields[1] == "run" and fields[2] == "index":
                return int(fields[3])
    return None

# Wrapper function for parsing the list of simulation temperatures from the temperature file
# input: t_file: text file holding list of molecule names and temperatures
#        name:  molecule name to return the temperatures for
def get_temps(t_file,name):
    if t_file is None:
        T = ["298.15"]
    else:
        T = None
        with open(t_file,'r') as f:
            for lines in f:
                fields = lines.split()
                if len(fields) > 0:
                    if fields[0] == name:
                        T = [ str(_) for _ in fields[1:] ]
                        break
        if T is None:
            T = ["298.15"]
    return T

# Wrapper function for parsing the list of simulation temperatures from the temperature file
# input: t_file: text file holding list of molecule names and temperatures
#        name:  molecule name to return the temperatures for
def get_rho(rho_file,name,T):
    rho = {}
    if rho_file is None:
        print('ERROR: density file not supplied, Exiting...')
        quit()
    else:
        for t in T:
           tmp_rho = None
           with open(rho_file,'r') as f:
               for lines in f:
                   fields = lines.split()
                   if len(fields) > 0:
                       if fields[0] == name and fields[1] == t:
                           tmp_rho = float(fields[2])
                           break
           if tmp_rho is None:
               tmp_rho = 200
           tmp_val = [ rho[key] for key in rho]
           if tmp_rho is None and len(tmp_val) != 0:
               tmp_rho = tmp_val[0]
           rho[t] = tmp_rho
    return rho

# Function for writing shell based input file
#def shell_submit(cmd,output,p,t,q,ppn=0,a=None,run_dir=None,batch="pbs"):
def shell_submit(cmd, output, p, t, a=None, ppn=0, batch="slurm", run_dir=None):
    N_nodes = int(np.ceil(float(p)/float(ppn)))
    if "min" in str(t):
        t = t.strip("min")
        min_flag = True
    else:
        min_flag = False
    if run_dir is None:
        run_dir = os.getcwd()
    with open(output,'w') as f:
        f.write("#!/bin/bash\n")
        f.write("#\n")

        if batch == "pbs":
            f.write("#PBS -N {} \n".format(output))
            if a is not None:
                f.write("#PBS -A {}\n".format(a))
            if ppn == 0:
                f.write("#PBS -l nodes={}\n".format(N_nodes))
            else:
                f.write("#PBS -l nodes={}:ppn={}\n".format(N_nodes,ppn))
            if min_flag is True:
                f.write("#PBS -l walltime=00:{}:00\n".format(t))
            else:
                f.write("#PBS -l walltime={}:00:00\n".format(t))
            f.write("#PBS -q {}\n".format(q))
            f.write("#PBS -S /bin/sh\n")

        if batch == "slurm":
            f.write("#SBATCH --job-name={}\n".format(output))
            f.write("#SBATCH -o {}out \n".format(output[:-6]))
            f.write("#SBATCH -e {}err \n".format(output[:-6]))
            if a is not None:
                f.write("#SBATCH -A {}\n".format(a))
            if ppn == 0:
                f.write("#SBATCH -N {}\n".format(N_nodes))
            else:
                f.write("#SBATCH -N {}\n".format(N_nodes))
                f.write("#SBATCH -n {}\n".format(ppn*N_nodes))
            if min_flag is True:
                f.write("#SBATCH -t  00:{}:00\n".format(t))
            else:
                f.write("#SBATCH -t {}:00:00\n".format(t))

        f.write("\n# cd into the submission directory\n")
        f.write("cd -- {}\n".format(run_dir))
        f.write("echo Working directory is {}\n".format(run_dir))
        f.write("echo Running on host `hostname`\n")
        f.write("echo Time is `date`\n\n")

        f.write("# Copy the local path\n")
        f.write("PATH={}\n".format(os.getenv("PATH")))
        f.write("sleep 2\n\n# Run Commands\n{}\nsleep 2\n".format(cmd))
    subprocess.call("chmod 777 {}".format(output), shell=True)

    return
                

# Inserts escape characters before parentheses so that the script can                                                                                                                                                                         
# handle filenames/folders that involve these                                                                                                                                                                                                 
def esc_str(a):
    a = a.replace('(','\\(')
    a = a.replace(')','\\)')
    return a

def num_frames(traj):
    N_frames = 0
    with open(traj,'r') as f:
        for lines in f:
            if "ITEM: TIMESTEP" in lines:
                N_frames += 1
    return N_frames

def run_dilatometry(c):

    # Make dilatometry directory
    if os.path.isdir("dilatometry") is False:
        os.mkdir("dilatometry")

    # Loop over pressures and create inputs if they don't exist
    print("\tgenerating dilatometry inputs...")
    for i in c["pressures"]:

        # Make the pressure folder if it doesn't exist
        if os.path.isdir("dilatometry/{}".format(i)) is False:
            os.mkdir("dilatometry/{}".format(i))

        # Make the runs if they don't exist
        for j in range(1,int(c["n_dil_runs"])+1):

            if os.path.isdir("dilatometry/{}/run_{}".format(i,j)) is False:
                substring="python {}/FF_functions/gen_md_for_sampling.py {} -N {} -gens {} -FF {} -pair_styles {} "+\
                          "-onefourscale_coul {} -onefourscale_lj {} -t {} -t_A {} -t_ext {} -T {} -T_A {} -T_D {} "+\
                          "-N_dil {} --tail --impropers --velocities -P {} -o {} --force_read -f {}"
                substring = substring.format(c["taffi_dir"],c["xyz"],c["num_mol"],c["gens"],c["ff"],','.join(c["pair_styles"].split()),c["14coul"],c["14lj"],c["t_dil_equil"],c["t_dil_anneal"],c["t_dil_ext"],\
                                             c["temp_dil_start"],c["temp_dil_start"],c["temp_dil_end"],c["dil_intervals"],i,"dilatometry/{}/run_{}".format(i,j),c["print_freq"])
                subprocess.check_output(subproc_string(substring), encoding='utf8')

    # Run initializations
    jobids = []
    print("\trunning dilatometry initializations...")
    for i in c["pressures"]:
        for j in range(1,int(c["n_dil_runs"])+1):
            if os.path.isfile("dilatometry/{}/run_{}/extend.end.restart".format(i,j)) is False:
                substring="python {}/Automation_Scripts/lammps_submit.py -f {}.in.init -d {},{},{} -p {} -ppn {} -size {} -q {} -sched {} -t {} -a {} --overwrite -path_to_exe {} -o {} -shell {} --silent"
                substring = substring.format(c["taffi_dir"],"run_{}".format(j),"dilatometry",i,"run_{}".format(j),c["procs"],c["ppn"],c["size"],c["queue"],c["sched"],c["wt"],\
                                             c["acct"],c["lammps_exe"],"{}_run_{}".format(i,j),','.join(c["shell_cmd"].split()))
                output = subprocess.check_output(subproc_string(substring), encoding='utf8')
                jobids += [output.split()[3]]

    monitor_jobs(jobids)

    # Run extensions
    for t in range(int(c["dil_intervals"])):
        
        jobids = []
        print("\trunning dilatometry interval {}...".format(t))
        for i in c["pressures"]:
            for j in range(1,int(c["n_dil_runs"])+1):

                # If the output trajectory for this interval is not present then attempt to run the job
                if os.path.isfile("dilatometry/{}/run_{}/{}.sys.lammpstrj".format(i,j,t)) is False:

                    # If the prerequisite restart file is missing, print diagnostic to user and exit
                    if os.path.isfile("dilatometry/{}/run_{}/extend.end.restart".format(i,j)) is False:
                        print("ERROR in dilatometry phase: the initialization run ({}) did not complete.".format("dilatometry/{}/run_{}".format(i,j)))
                        quit()

                    # If the trajectory from the previous interval is missing, print diagnostic to user and exit.
                    elif t > 0 and os.path.isfile("dilatometry/{}/run_{}/{}.sys.lammpstrj".format(i,j,t-1)) is False:
                        print("ERROR in dilatometry phase: the extension run ({}) did not complete.".format("dilatometry/{}/run_{}/{}.sys.lammpstrj".format(i,j,t)))
                        quit()

                    # Else, everything is in place and the job is run.
                    else:
                        substring="python {}/Automation_Scripts/lammps_submit.py -f {}.in.init -d {},{},{} -p {} -ppn {} -size {} -q {} -sched {} -t {} -a {} --overwrite -path_to_exe {} -o {} -shell {} --silent"
                        substring = substring.format(c["taffi_dir"],"extend","dilatometry",i,"run_{}".format(j),c["procs"],c["ppn"],c["size"],c["queue"],c["sched"],c["wt"],\
                                                     c["acct"],c["lammps_exe"],"{}_run_{}".format(i,j),','.join(c["shell_cmd"].split()))
                        output = subprocess.check_output(subproc_string(substring), encoding='utf8')
                        jobids += [output.split()[3]]

        monitor_jobs(jobids)
        
        # Check that all of the runs completed
        for i in c["pressures"]:
            for j in range(1,int(c["n_dil_runs"])):

                # If the trajectory from the previous interval is missing, print diagnostic to user and exit.
                if os.path.isfile("dilatometry/{}/run_{}/{}.sys.lammpstrj".format(i,j,t)) is False:
                    print("ERROR in dilatometry phase: the extension run ({}) did not complete.".format("dilatometry/{}/run_{}/{}.sys.lammpstrj".format(i,j,t)))
                    quit()

    # Combine trajectories
    print("\tcombining dilatometry trajectories...")
    current_dir = os.getcwd()
    for i in c["pressures"]:
        for j in range(1,int(c["n_dil_runs"])+1):

            # Remove incomplete combinations if present (e.g., if the use aborted the script execution in this step)
            if os.path.isfile("dilatometry/{}/run_{}/combined.lammpstrj".format(i,j)):
                os.remove("dilatometry/{}/run_{}/combined.lammpstrj".format(i,j))
            if os.path.isfile("dilatometry/{}/run_{}/combined.thermo".format(i,j)):
                os.remove("dilatometry/{}/run_{}/combined.thermo".format(i,j))

            # Combine the trajectory and thermo files if they haven't been combined
            if os.path.isfile("dilatometry/{}/run_{}/{}_run_{}_combined.lammpstrj".format(i,j,i,j)) is False:

                # Run combine script
                os.chdir("dilatometry/{}/run_{}".format(i,j))
                substring="{}/Automation_Scripts/combine_ext.sh *.sys.lammpstrj 0 {}".format(c["taffi_dir"],int(c["dil_intervals"])-1)
                output = subprocess.check_output(subproc_string(substring), encoding='utf8')
                substring="{}/Automation_Scripts/combine_ext.sh *.thermo.avg 0 {} --thermo".format(c["taffi_dir"],int(c["dil_intervals"])-1)
                output = subprocess.check_output(subproc_string(substring), encoding='utf8')
                shutil.move("combined.lammpstrj","{}_run_{}_combined.lammpstrj".format(i,j))
                shutil.move("combined.thermo","{}_run_{}_combined.thermo".format(i,j))
                os.chdir(current_dir)

    # Make run folders if they don't exist 
    if os.path.isdir("Results") is False:
        os.mkdir("Results")
    if os.path.isdir("Results/dilatometry") is False:
        os.mkdir("Results/dilatometry")    

    # Process dilatometry data and calculate approximate Tg values
    data={}
    Tg={}
    print("\tprocessing dilatometry data...")
    for i in c["pressures"]:
        data[i] = {}
        for j in range(1,int(c["n_dil_runs"])+1):            
            data[i][j] = { "T":[],"d":[],"U":[],"H":[] }
            if os.path.isfile("Results/dilatometry/{}_run_{}_smooth.thermo".format(i,j)) is False:

                # Smooth the data
                substring="python {}/smooth_data.py -f {} -o Results/dilatometry/{}_run_{}_smooth.thermo -start 0 -end {} -head 2 -x_col 1 -y_col 2,4,8 -x_min {} -x_max {} -x_step 10 -w 100 "+\
                          "-labels T(K),density(g/cm^3),PE(kcal/mol),H(kcal/mol) -y_scale 1.0,{},{}"
                substring = substring.format(c["colops_dir"],"dilatometry/{}/run_{}/{}_run_{}_combined.thermo".format(i,j,i,j),i,j,int(float(c["t_dil_ext"])*float(c["dil_intervals"])/float(c["print_freq"])),\
                                             c["temp_dil_end"],c["temp_dil_start"],1.0/float(c["num_mol"]),1.0/float(c["num_mol"]))
                output = subprocess.check_output(subproc_string(substring), encoding='utf8')

            # Collect the smoothed data
            data[i][j]["T"],data[i][j]["d"],data[i][j]["U"],data[i][j]["H"] = get_cols("Results/dilatometry/{}_run_{}_smooth.thermo".format(i,j),[0,1,2,3],1)

            # Make the Tvd TvU and TvH plots (if necessary)
            T_raw,d_raw,U_raw,H_raw = get_cols("dilatometry/{}/run_{}/{}_run_{}_combined.thermo".format(i,j,i,j),[1,2,4,8],2)
            U_raw = np.array(U_raw)/float(c["num_mol"])
            H_raw = np.array(H_raw)/float(c["num_mol"])

            if os.path.isfile("Results/dilatometry/{}_run_{}_Tvd_comparison.pdf".format(i,j)) is False:
                make_plot("Results/dilatometry/{}_run_{}_Tvd_comparison.pdf".format(i,j),[T_raw,data[i][j]["T"]],[d_raw,data[i][j]["d"]],x_label="T (K)",y_label="d (g/cm^3)")            
            if os.path.isfile("Results/dilatometry/{}_run_{}_TvU_comparison.pdf".format(i,j)) is False:
                make_plot("Results/dilatometry/{}_run_{}_TvU_comparison.pdf".format(i,j),[T_raw,data[i][j]["T"]],[U_raw,data[i][j]["U"]],x_label="T (K)",y_label="U (kcal/mol)")
            if os.path.isfile("Results/dilatometry/{}_run_{}_TvH_comparison.pdf".format(i,j)) is False:
                make_plot("Results/dilatometry/{}_run_{}_TvH_comparison.pdf".format(i,j),[T_raw,data[i][j]["T"]],[H_raw,data[i][j]["H"]],x_label="T (K)",y_label="H (kcal/mol)")

        # Make plots with all runs in the same figure
        if os.path.isfile("Results/dilatometry/{}_Tvd_runs.pdf".format(i)) is False:        
            make_plot("Results/dilatometry/{}_Tvd_runs.pdf".format(i),[ data[i][_]["T"] for _ in sorted(data[i].keys()) ],[ data[i][_]["d"] for _ in sorted(data[i].keys()) ],x_label="T (K)",y_label="d (g/cm^3)")                    
        if os.path.isfile("Results/dilatometry/{}_TvU_runs.pdf".format(i)) is False:
            make_plot("Results/dilatometry/{}_TvU_runs.pdf".format(i),[ data[i][_]["T"] for _ in sorted(data[i].keys()) ],[ data[i][_]["U"] for _ in sorted(data[i].keys()) ],x_label="T (K)",y_label="U (kcal/mol)")
        if os.path.isfile("Results/dilatometry/{}_TvH_runs.pdf".format(i)) is False:
            make_plot("Results/dilatometry/{}_TvH_runs.pdf".format(i),[ data[i][_]["T"] for _ in sorted(data[i].keys()) ],[ data[i][_]["H"] for _ in sorted(data[i].keys()) ],x_label="T (K)",y_label="H (kcal/mol)")

        # Calculate averages and stdevs across the runs            
        for k in ["d","U","H"]:
            for j in range(1,int(c["n_dil_runs"])+1):

                # Initialize arrays if this is the first being read
                if j == 1:
                    x_mean = np.array(data[i][j]["T"])
                    y_mean = np.array(data[i][j][k])
                    y_std  = np.array(data[i][j][k])**(2.0)

                # Else, add to arrays
                else:
                    x_mean += np.array(data[i][j]["T"])
                    y_mean += np.array(data[i][j][k])
                    y_std  += np.array(data[i][j][k])**(2.0)

            # Calculate quantities
            x_mean = x_mean / float(int(c["n_dil_runs"]))
            y_mean = y_mean / float(int(c["n_dil_runs"]))
            y_std  = ( y_std  / float(int(c["n_dil_runs"])) - y_mean**(2.0) )**(0.5)

            # Make plot and calculate derivatives
            if k == "d":
                x1,y1,x2,y2,x_disc = calc_discontinuity(x_mean,y_mean)
                Tg[i] = x_disc
                if os.path.isfile("Results/dilatometry/{}_Tv{}_mean.pdf".format(i,k)) is False:
                    make_plot("Results/dilatometry/{}_Tv{}_mean.pdf".format(i,k),[x_mean,x1,x2],[y_mean,y1,y2],y_errs=[y_std,None,None],x_label="T (K)",y_label="d (g/cm^3)",color_inds=[0,3,3],\
                              linewidths=[5,2,2],v_line=[x_disc],labels=["mean","lin_left","lin_right"])
            if k == "H" and os.path.isfile("Results/dilatometry/{}_Tv{}_mean.pdf".format(i,k)) is False:
                make_plot("Results/dilatometry/{}_Tv{}_mean.pdf".format(i,k),[x_mean],[y_mean],y_errs=[y_std],x_label="T (K)",y_label="H (kcal/mol)",labels=["mean"])
            if k == "U" and os.path.isfile("Results/dilatometry/{}_Tv{}_mean.pdf".format(i,k)) is False:
                make_plot("Results/dilatometry/{}_Tv{}_mean.pdf".format(i,k),[x_mean],[y_mean],y_errs=[y_std],x_label="T (K)",y_label="U (kcal/mol)",labels=["mean"])                    
    return
                
# Calculates the Tg based on an algorithm that splits the data in two halves and fits the 
# low temperature region to a linear function and high temperature region to a quadratic function.
# Tg is taken as the split point that minimizes the error of the fit.
def calc_discontinuity(x,y):

    # Collected sorted values as arrays
    x,y = list(zip(*sorted([ (i,y[count_i]) for count_i,i in enumerate(x) ])))
    x = np.array(x)
    y = np.array(y)
    N = len(x)

    # Calculate linear fits on each half of the data
    best_err = None
    for i in range(2,N-2):

        # Left fit (linear only)
        A = np.array([ [ 1.0,_] for _ in x[:i] ])
        intercept1,slope1 = lstsq(A,y[:i],rcond=None)[0]
        err_low = np.mean((y[:i] - (slope1*x[:i] + intercept1))**(2.0))

        # Right fit (linear and quandratic terms)
        A = np.array([ [ 1.0,_,_**(2.0)] for _ in x[i:] ])
        intercept2,slope2,quad2 = lstsq(A,y[i:],rcond=None)[0]
        err_high = np.mean((y[i:] - (slope2*x[i:] + intercept2 + quad2*x[i:]**2.0 ))**(2.0))

        # Calculate total error with weighting based on the number of points in each set
        tot_err = err_low * float(i)/float(N) + err_high * float(N-i)/float(N)
        if best_err is None or tot_err < best_err:
            best_err = tot_err
            crossing_point = (intercept2-intercept1)/(slope1-slope2)
            crossing_ind = i
            b_s1 = slope1
            b_i1 = intercept1
            b_s2 = slope2
            b_i2 = intercept2

            b_q1 = 0
            b_q2 = quad2

    # Calculate the interpolated sets
    y1 = np.array([ (b_i1 + b_s1*x[j] + b_q1*x[j]**(2.0)) for j in range(crossing_ind) ])
    y2 = np.array([ (b_i2 + b_s2*x[j] + b_q2*x[j]**(2.0)) for j in range(crossing_ind,N) ])
    x1 = np.array([ x[j] for j in range(crossing_ind) ])
    x2 = np.array([ x[j] for j in range(crossing_ind,N) ])
    
    return x1,y1,x2,y2,x[crossing_ind]

# Fits linear data with a weight function
def disc_err(params,x,y,w):
    return np.mean(((params[0]*x+params[1]) - y)**2.0 * w)

def disc_quad(params,x,y,w):
    return np.mean((params[0] + params[1]*x + params[2]*x**(2) - y)**2.0 * x )

# Reads in column data from file
def get_cols(name,cols,head):
    data = { i:[] for i in cols }
    with open(name,'r') as f:
        for lc,lines in enumerate(f):
            if lc >= head:
                fields = lines.split()
                for i in cols:
                    data[i] += [float(fields[i])]
    return [ data[i] for i in cols ]

# Formats a subprocess string for execution by subprocess.check_output
def subproc_string(s):

    s = s.split()

    # Replace , arguments with spaces
    for count_i,i in enumerate(s):
        if "^" in i:
            s[count_i] = s[count_i].replace('^',' ')
    
    return s

# Function that sleeps the script until jobids are no longer in a running or pending state in the queue
#def monitor_jobs(jobids):
#    
#    current_jobs = check_queue()
#    while True in [ i in current_jobs for i in jobids ]:
#        time.sleep(60)
#        current_jobs = check_queue()  
#    return

# Returns the pending and running jobids for the user as a list
def check_queue():

    # The first time this function is executed, find the user name and scheduler being used. 
    if not hasattr(check_queue, "user"):

        # Get user name
        check_queue.user = subprocess.check_output("echo ${USER}", shell=True, encoding='utf8').strip("\r\n")

        # Get batch system being used

        squeue_tmp = subprocess.Popen(['which', 'squeue'], stdout=subprocess.PIPE,stderr=subprocess.STDOUT, encoding='utf8').communicate()[0].strip("\r\n")
        qstat_tmp = subprocess.Popen(['which', 'qstat'], stdout=subprocess.PIPE,stderr=subprocess.STDOUT, encoding='utf8').communicate()[0].strip("\r\n")
        check_queue.sched =  None
        if "no squeue in" not in squeue_tmp:
            check_queue.sched = "slurm"
        elif "no qstat in" not in qstat_tmp:
            check_queue.sched = "pbs"
        else:
            print("ERROR in check_queue: neither slurm or pbs schedulers are being used.")
            quit()

    # Get running and pending jobs using the slurm scheduler
    if check_queue.sched == "slurm":

        # redirect a squeue call into output
        output = subprocess.check_output("squeue -l", shell=True, encoding='utf8')

        # Initialize job information dictionary
        jobs = []
        id_ind = None
        for count_i,i in enumerate(output.split('\n')):            
            fields = i.split()
            if len(fields) == 0: continue                
            if id_ind is None and "JOBID" in fields:
                id_ind = fields.index("JOBID")
                if "STATE" not in fields:
                    print("ERROR in check_queue: Could not identify STATE column in squeue -l output.")
                    quit()
                else:
                    state_ind = fields.index("STATE")
                if "USER" not in fields:
                    print("ERROR in check_queue: Could not identify USER column in squeue -l output.")
                    quit()
                else:
                    user_ind = fields.index("USER")
                continue

            # If this job belongs to the user and it is pending or running, then add it to the list of active jobs
            if id_ind is not None and fields[user_ind] == check_queue.user and fields[state_ind] in ["PENDING","RUNNING"]:
                jobs += [fields[id_ind]]

    # Get running and pending jobs using the pbs scheduler
    elif check_queue.sched == "pbs":

        # redirect a qstat call into output
        output = subprocess.check_output("qstat -f", shell=True, encoding='utf8')

        # Initialize job information dictionary
        jobs = []
        job_dict = {}
        current_key = None
        for count_i,i in enumerate(output.split('\n')):
            fields = i.split()
            if len(fields) == 0: continue
            if "Job Id" in i:

                # Check if the previous job belongs to the user and needs to be added to the pending or running list. 
                if current_key is not None:
                    if job_dict[current_key]["State"] in ["R","Q"] and job_dict[current_key]["User"] == check_queue.user:
                        jobs += [current_key]
                current_key = i.split()[2]
                job_dict[current_key] = { "State":"NA" , "Name":"NA", "Walltime":"NA", "Queue":"NA", "User":"NA"}
                continue
            if "Job_Name" == fields[0]:
                job_dict[current_key]["Name"] = fields[2]
            if "job_state" == fields[0]:
                job_dict[current_key]["State"] = fields[2]
            if "queue" == fields[0]:
                job_dict[current_key]["Queue"] = fields[2]
            if "Resource_List.walltime" == fields[0]:
                job_dict[current_key]["Walltime"] = fields[2]        
            if "Job_Owner" == fields[0]:
                job_dict[current_key]["User"] = fields[2].split("@")[0]

        # Check if the last job belongs to the user and needs to be added to the pending or running list. 
        if current_key is not None:
            if job_dict[current_key]["State"] in ["R","Q"] and job_dict[current_key]["User"] == check_queue.user:
                jobs += [current_key]

    return jobs

# Plotting script
def make_plot(name,x_sets,y_sets,y_errs=[],dims=(4.5,4),x_label="x_label",y_label="y_label",labels=None,log_x=False,log_y=False,x_min=None,x_max=None,y_min=None,y_max=None,\
              tick_scale=1.0,x_tick_inc=None,y_tick_inc=None,diagonal=False,no_legend=False,linewidths=[3.0],color_inds=[],v_line=[]):
    
    color_list = [(0.05,0.35,0.75),(0.05,0.8,0.6),(0.9,0.3,0.05),(0.05,0.05,0.05),(0.9,0.6,0.05),(0.9,0.5,0.7),(0.35,0.7,0.9),(0.95,0.9,0.25)]   # Supposedly more CB-friendly
    if len(x_sets) != len(y_sets):
        print("ERROR in make_plot: len(x_sets) must equal len(y_sets)")
        quit()
    if False in [ len(i) == len(y_sets[count_i]) for count_i,i in enumerate(x_sets) ]:
        print("ERROR in make_plot: dimension mismatch in x and y sets.")
        quit()
    if labels is None:
        labels=[ i for i in range(len(x_sets)) ]
    elif len(labels) != len(x_sets):
        print("ERROR in make_plot: len(labels) is not equal to len(x_sets) and len(y_sets).")
        quit()
    elif len(y_errs) != 0 and len(y_errs) != len(y_sets):
        print("ERROR in make_plot: len(y_errs) is not equal to len(x_sets) and len(y_sets).")
        quit()
    while len(color_inds) < len(x_sets):
        color_inds += [len(color_inds)]
    while len(linewidths) < len(x_sets):
        linewidths += [linewidths[-1]]

    # Generate a plot with all individual datasets
    fig = plt.figure(figsize=dims)
    ax = plt.subplot(111)
    plot_handles = []

    # Plot error region
    for count_i,i in enumerate(y_errs):
        if i is not None:
            ax.fill_between(x_sets[count_i],\
                            (np.array(y_sets[count_i])-np.array(i)),\
                            (np.array(y_sets[count_i])+np.array(i)),\
                            facecolor=color_list[color_inds[count_i]],edgecolor='none',alpha=0.25)

    for count_i,i in enumerate(x_sets):
        plot_handle, = ax.plot(np.array(i),np.array(y_sets[count_i]),\
                               linestyle="-",linewidth=linewidths[count_i],\
                               color=color_list[color_inds[count_i]],\
                               label=labels[count_i])
        
        plot_handle.set_dash_capstyle('round')
        plot_handle.set_dash_joinstyle('round')
    #ax.set_ylabel(y_label.decode('utf-8'),fontsize=32,labelpad=10,fontweight='bold')
    #ax.set_xlabel(x_label.decode('utf-8'),fontsize=32,labelpad=10,fontweight='bold')
    ax.set_ylabel(y_label,fontsize=32,labelpad=10,fontweight='bold')
    ax.set_xlabel(x_label,fontsize=32,labelpad=10,fontweight='bold')

    # Set Logarithmic Bounds
    if log_y:
        ax.set_yscale('log')
    if log_x:
        ax.set_xscale('log')
                               
    # Set limits
    c_y_min,c_y_max = ax.get_ylim()
    c_x_min,c_x_max = ax.get_xlim()
    if x_min is not None: c_x_min = x_min
    if x_max is not None: c_x_max = x_max
    if y_min is not None: c_y_min = y_min
    if y_max is not None: c_y_max = y_max
    ax.set_xlim([c_x_min,c_x_max])
    ax.set_ylim([c_y_min,c_y_max])

    # Set tick incremement
    if x_tick_inc is not None:
        plt.xticks(np.arange(c_x_min, c_x_max+x_tick_inc/1E8, x_tick_inc))
    if y_tick_inc is not None:
        plt.yticks(np.arange(c_y_min, c_y_max+y_tick_inc/1E8, y_tick_inc))
                               
    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=24*tick_scale,pad=10,direction='out',width=3,length=6)
    ax.tick_params(axis='both', which='minor',labelsize=24*tick_scale,pad=10,direction='out',width=2,length=4)
    [j.set_linewidth(3) for j in ax.spines.values()]
                               
    # Add diagonal
    if diagonal is True:
        ax.plot(ax.get_xlim(), ax.get_ylim(), ls="-", c=".3")

    # Add vertical lines
    for i in v_line:
        ax.plot([i,i],ax.get_ylim(), ls="--", c=".3")

    # Put a legend to the right of the current axis and use modified save call to account for the legend artist
    if no_legend is False:

        # Place legend outside to the right by default
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1,0.5))
        plt.savefig(name, dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')            

    # If no legend is requested then place to use ordinary tight_layout and save calls
    else:
        plt.tight_layout()
        plt.savefig(name, bbox_inches=0,dpi=300)
    plt.close(fig)
    return

# Function for keeping tabs on the validity of the user supplied inputs
def parse_configuration(args):

    # Convert inputs to the proper data type
    if os.path.isfile(args.config) is False:
        print("ERROR in python_driver: the configuration file {} does not exist.".format(args.config))
        quit()

    # Process configuration file for keywords
    keywords = [ "FOM_GEOOPT_PROCS", "FOM_GEOOPT_WT", "FOM_GEOOPT_Q", "FOM_GEOOPT_SCHED", "FOM_GEOOPT_PPN", "FOM_GEOOPT_SIZE",\
                 "FOM_MD_GAS_PROCS", "FOM_MD_GAS_WT", "FOM_MD_GAS_Q", "FOM_MD_GAS_NPP", "FOM_MD_GAS_SCHED", "FOM_MD_GAS_PPN",\
                 "FOM_MD_GAS_SIZE", "FOM_MD_COND_PROCS", "FOM_MD_COND_WT", "FOM_MD_COND_Q", "FOM_MD_COND_NPP", "FOM_MD_COND_SCHED",\
                 "FOM_MD_COND_PPN", "FOM_MD_COND_SIZE", "FOM_TEMPS", "LAMMPS_EXE", "ORCA_EXE", "TAFFI_PATH", "fom_ff", 'batch', "FOM_N_TRAJ", "acct", \
                 "fom_t_ext", "fom_t_tot", "fom_t_a","fom_restarts", "fom_t_equil","CHARGE","GENS","FOM_DENS","FUNCTIONAL","D3_option" ]
    keywords = [ _.lower() for _ in keywords ]

    list_delimiters = [ "," ]  # values containing any delimiters in this list will be split into lists based on the delimiter
    space_delimiters = [ "&" ] # values containing any delimiters in this list will be turned into strings with spaces replacing delimiters
    configs = { i:None for i in keywords }    
    with open(args.config,'r') as f:
        for lines in f:
            fields = lines.split()
            
            # Delete comments
            if "#" in fields:
                del fields[fields.index("#"):]

            # Parse keywords
            l_fields = [ _.lower() for _ in fields ] 
            for i in keywords:
                if i in l_fields:

                    # Parse keyword value pair
                    ind = l_fields.index(i) + 1
                    if len(fields) >= ind + 1:
                        configs[i] = fields[ind]

                        # Handle delimiter parsing of lists
                        for j in space_delimiters:
                            if j in configs[i]:
                                configs[i] = " ".join([ _ for _ in configs[i].split(j) ])
                        for j in list_delimiters:
                            if j in configs[i]:
                                configs[i] = configs[i].split(j)
                                break
                                
                    # Break if keyword is encountered in a non-comment token without an argument
                    else:
                        print("ERROR in python_driver: enountered a keyword ({}) without an argument.".format(i))
                        quit()
    if configs["functional"] == 'wB97X-D3':
       configs["d3_option"] = "--no_D3"
    else:
       configs["d3_option"] = ''

    # Set defaults if None
    if configs["taffi_path"] is None:
        configs["taffi_path"] = '/'.join(os.path.abspath(__file__).split('/')[:-2])
    if configs["fom_n_traj"] is None:
        configs["fom_n_traj"] = 5
    else:
        configs["fom_n_traj"] = int(configs["fom_n_traj"])
    if configs["fom_restarts"] is None:
        configs["fom_restarts"] = 0
    else:
        configs["fom_restarts"] = int(configs["fom_restarts"])
        if configs["fom_restarts"] < 0:
            print("ERROR: fom_restarts value must be >= 0.")
            quit()

    # Collect the xyz files
    configs["xyz"] = [ f for f in os.listdir('.') if os.path.isfile(f) and ".xyz" in f ]

    # Check that batch is an acceptable system
    if configs["batch"] not in ["pbs","slurm"]:
        print("ERROR in FOM_driver: only pbs and slurm are acceptable arguments for the batch variable.")
        quit()
    elif configs["batch"] == "pbs":
        configs["sub_cmd"] = "qsub"
    elif configs["batch"] == "slurm":
        configs["sub_cmd"] = "sbatch"

    return configs

# Create logger to save stdout to logfile
class Logger(object):

    def __init__(self,logname):
        self.terminal = sys.stdout
        self.log = open(logname+'.log', "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass

if __name__ == "__main__":
    main(sys.argv[1:])
