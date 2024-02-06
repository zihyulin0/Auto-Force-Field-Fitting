#!/bin/env python                                                                                                                                                              
import sys,os,argparse,subprocess,shutil,time,matplotlib,glob,getpass,json,fnmatch
import random

# For plotting (Agg called needed for cluster image generation)
matplotlib.use('Agg') 
from pylab import *
import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit,minimize,lsq_linear
from copy import deepcopy
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/FF_functions')
import log_to_db,merge_FF_drude,db_to_charmm,relax_drude_loop
from monitor_jobs import *
from parse_all import *
from file_parsers import *
import place_charge,generate_connolly,generate_connolly_new,gen_md_for_CHARMM_new
import codecs,json
# Add TAFFY Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from transify import *
from adjacency import *
from file_parsers import *
from id_types import *
from kekule import *


def main(argv):
    

    global user,run_dir,c
    user = getpass.getuser()
    current_dir = os.getcwd()

    

    parser = argparse.ArgumentParser(description='Driver script for taffi benchmark simulations. Expects a configuation file formatted with space delimited keyword arguments')

    #optional arguments                                                                                                                                                        
    parser.add_argument('-c', dest='config', default='config.txt',
                        help = 'The program expects a configuration file from which to assign various run conditions. (default: config.txt in the current working directory)')
    parser.add_argument('-alpha_scale', dest='alpha_scale', default=1,
                        help = 'options for charmm to rescale atomic polarizability when fitting, default: 1')
    parser.add_argument('-outFF', dest='outFF', default='fit_CHARMM.db',
                        help = 'output db filename, default: fit_CHARMM.db')
    parser.add_argument('-folder', dest='folder', default='',
                        help = 'inchikey folder name')
    parser.add_argument('-o',dest='output',default='polar.log',
                        help = 'log name for stdout')
    parser.add_argument('-tag',dest='tag',default='',
                        help = '_tag name default: None')
    parser.add_argument('-nonpolar_dir',dest='nonpolar_dir',default='',
                        help = 'nonpolar directory')
    parser.add_argument('--write_bash', dest='write_bash', default=False, action='store_const', const=True,
                        help = 'When this flag is supplied, will only write out bash with generation info') 

    args=parser.parse_args()    

    if args.nonpolar_dir == '':
      print('ERROE: nonpolar directory not specified, Exiting....')
      quit()

    if args.write_bash:
       # Read all.json (contains model compound and dependency information for the batch 
       status = read_alljson('./all.json')
       # Find inchi keys for this generation
       mc_keys = [ _ for _ in status["mc"].keys() ]
       # Find inchi keys for this generation
       #mc_keys = [ _ for _ in status["mc"].keys() if int(g) in status["mc"][_]["gen"] ]
       gens = sorted(status["gens"].keys())
       with open('drude_gens.sh','w') as f:
         f.write("#!/bin/bash\n")
         for g in gens:
            keys = [ _ for _ in status["mc"].keys() if int(g) == status["mc"][_]["gen"][0] ] 
            f.write("gen{}=(".format(g))
            for i in keys:
               f.write('\'{}\' '.format(i))
            f.write(')\n')
         f.write("\nfor i in \"${gen1[@]}\"; do \n")
         f.write("   screen -X -S $i quit \n")
         f.write("   screen -dmS $i \n")
         f.write("   screen -S $i -X stuff \"ana3 ; ~/taffi_beta/Automation_Scripts/drude_configs.py -nonpolar_dir {} -c {} -folder $i \\n\" \n".format(args.nonpolar_dir,args.config))
         f.write("done \n")
       subprocess.call("chmod u+x drude_gens.sh",shell=True)
       quit()
    
    sys.stdout = Logger(current_dir+'/'+args.folder+'/'+args.output)
    
    # parse configuration dictionary (c)
    print("parsing configurations...")
    c = parse_configuration(args.config)
    # Make config file absolute
    args.config = os.path.abspath(args.config)
    if (c['module_string'] == None):
      c['module_string'] = ''


    #fit_Drude(run_dir,args.folder,c,args,user,args.tag)
    fit_Drude_custom(current_dir,args.folder,args,args.tag)

    return


# Normal Drude fitting, fit all heavy atoms
def fit_Drude_custom(driver_dir,folder,args,tag):

    global run_dir

    run_dir = '{}/{}'.format(driver_dir,folder)
    os.chdir(run_dir) 
 
    nonpolar = args.nonpolar_dir 
    c_config ='{}/{}/charges/CHELPG_calcs/configs'.format(nonpolar,folder)
    nconf = sorted([int(f) for f in os.listdir(c_config) if os.path.isdir(os.path.join(c_config, f))])[-1]+1
    if nconf > 100: nconf = 100
   
    # placing perturbation charges 
    jobids = []
    print("Start placing perturbation charges...")

    # try five times incase something went wrong
    for retry in range(0,5):
       for count_i,i in enumerate(list(range(0,nconf))):
         fitfolder = 'fitcharge_{}'.format(i)
         if os.path.isdir(run_dir+'/'+fitfolder) is False:
            os.makedirs(run_dir+'/'+fitfolder)
         xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
         element,geo = xyz_parse(xyz)
         jobids += place_perturbation_submit(fitfolder,xyz,dplus=retry*100)
       if jobids == []: break
       monitor_jobs(jobids,user)
    clean_pcsubmit(nconf)

    # do QM calculation of polarizability and dipole for the origianl xyz for later FOM to compare to
    print("Starting QM polar calculation...")
    for count_i,i in enumerate(list(range(0,nconf))):
      fitfolder = 'fitcharge_{}'.format(i)
      xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
      element,geo = xyz_parse(xyz)
      polar_orca(fitfolder,geo,element)
  
    os.chdir(run_dir) 
    substring = "{}/Automation_Scripts/orca_submit.py -f polar.in -p {}  -t {} -ppn {} -q {} -sched {} -size {} -o optimizations -path_to_exe {}  --silent"
    substring = substring.format(c["taffi_path"],c['charges_qc_procs'],c["charges_qc_wt"],c["charges_qc_ppn"],c["charges_qc_q"],c["charges_qc_sched"],c["charges_qc_size"],c["orca_exe"])
    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    output = str(output,'utf-8')
    orcaids = [ m.split()[-1] for m in output.split("\n")[:-1]]
    
    print("Starting CHELPG potential calculation...")
    for retry in range(0,5):
       jobids = []
       for q in [-0.5,0.5]:
          for count_i,i in enumerate(list(range(0,nconf))):
            fitfolder = 'fitcharge_{}'.format(i)
            xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
            element,geo = xyz_parse(xyz)
            total_pc,pc_layer1,pc_layer2,pc_layer3 = get_pc(fitfolder)
            jobids += calc_CHELPG(fitfolder,xyz,total_pc,pc_layer1,pc_layer2,pc_layer3,geo,element,q,'_{}'.format(q))
       monitor_jobs(jobids,user)
       clean_chelpg(nconf)

    print("Starting reevaluation of potential at Connolly grid...")
    for retry in range(0,5):
       jobids = []
       for q in [-0.5,0.5]:
          for count_i,i in enumerate(list(range(0,nconf))):
            fitfolder = 'fitcharge_{}'.format(i)
            xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
            total_pc,pc_layer1,pc_layer2,pc_layer3 = get_pc(fitfolder)
            connolly_grid = get_connolly_grid(fitfolder,xyz)
            jobids += calc_connolly_pot(fitfolder,total_pc,connolly_grid,c,'_{}'.format(q))
       monitor_jobs(jobids,user)
       clean_connolly(nconf)
   
    print("Starting rewrite Connolly potential format...")
    for q in [-0.5,0.5]:
       for count_i,i in enumerate(list(range(0,nconf))):
         fitfolder = 'fitcharge_{}'.format(i)
         xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
         total_pc,pc_layer1,pc_layer2,pc_layer3 = get_pc(fitfolder)
         connolly_grid = get_connolly_grid(fitfolder,xyz)
         rewrite_pot_wrap(fitfolder,total_pc,connolly_grid,'_{}'.format(q))
   
    print("Starting writing charmm data...") 
    # Combine data to data_ dir
    for q in [-0.5,0.5]:
       for count_i,i in enumerate(list(range(0,nconf))):
         fitfolder = 'fitcharge_{}'.format(i)
         xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
         element,geo = xyz_parse(xyz)
         total_pc,pc_layer1,pc_layer2,pc_layer3 = get_pc(fitfolder)
         write_charmm_data(geo,fitfolder,fitfolder,list(range(0,total_pc)),'_{}'.format(q))
       write_datadir('_{}'.format(q))

    # tar configs files to save space
    for count_i,i in enumerate(list(range(0,nconf))):
      fitfolder = 'fitcharge_{}'.format(i)
      tar_file(fitfolder)

    monitor_jobs(orcaids,user)
    clean_polar(nconf)
    
    """ 
      bond_lp = get_samples('{}/{}'.format(run_dir,fitfolder))
      pc_layer = [pc_layer1,pc_layer2,pc_layer3]
      sample = []
      for count_j,j in enumerate(pc_layer):
         start_pc = int(sum([k for count_k,k in enumerate(pc_layer) if count_k < count_j]))
         end_pc = int(sum([k for count_k,k in enumerate(pc_layer) if count_k <= count_j]))
         sample += list(range(start_pc+bond_lp[count_j],end_pc)) 
    
      #change_npert(run_dir,total_pc,fitfolder)
      #change_npert(run_dir,len(sample),fitfolder)
      #write_charmm_data(run_dir,0,0,geo,fitfolder,fitfolder,list(range(0,total_pc)),'_neg')
      #write_charmm_data(run_dir,0,0,geo,fitfolder,fitfolder,sample,'')
      os.chdir(run_dir+'/'+fitfolder)
      #polar_orca(geo,element)
      #fit_exe('{}_{}'.format(title,i))
      #continue
      command = ('-logname {}_{}.log -o {}_{}.db -xyz {}'.format(title,i,title,i,xyz)).split()
      tmpchi,tmp_molpolar = log_to_db.main(command)
      QM_dipole,QM_molpolar = read_polarout("polar.out")  
      chi2 += [tmpchi/total_pc]
      tmp = read_db('{}_{}.db'.format(title,i))
      f.write("{:<20s} {:<20.5f} {:<20.5f} {:<20.5f} ".format('{}_{}'.format(title,i),tmpchi/total_pc,tmp_molpolar,QM_molpolar))
      for j in atom_types:
         if j.split('[')[1] != '1':
            f.write("{:<20.5f} {:<20.5f} ".format(tmp[0][j],tmp[1][j])) 
         else:
            f.write("{:<20.5f} ".format(tmp[0][j])) 
      f.write("\n")
    #os.chdir(run_dir)
    #substring = "python3 {}/Automation_Scripts/bundle_submit.py  -f fit*.in -p {} -t {} -ppn {} -q {} -sched {} -size {} -o fit -path_to_exe {}  --silent"
    #substring = substring.format(c["taffi_path"],1,c["charges_qc_wt"],c["charges_qc_ppn"],c["charges_qc_q"],c["charges_qc_sched"],c["charges_qc_ppn"],c["orca_exe"])
    #output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    """ 
    delete_wildcard('','optimization')
    print("Drude QM calculation succeed!")
    quit()
   
    return

def tar_file(folder):
    os.chdir(run_dir+'/'+folder)
    ori_dir = os.getcwd()

    folder_names = ['configs_-0.5','configs_0.5','configs_connolly_-0.5','configs_connolly_0.5']
    for i in folder_names:
      subprocess.call('tar zcvf {}.tar.gz {}'.format(i,i),shell=True)
      subprocess.call('rm -r  {}'.format(i),shell=True)
    subprocess.call('tar zcvf pointcharge.tar.gz pointcharge_layer*'.format(i,i),shell=True)
    subprocess.call('rm -r  pointcharge_layer*',shell=True)
     
    os.chdir(ori_dir)

    return

def clean_polar(nconf):
   for count_i,i in enumerate(list(range(0,nconf))):
      fitfolder = 'fitcharge_{}'.format(i)
      delete_wildcard('{}'.format(fitfolder),'pc_chelpg')
      delete_wildcard('{}'.format(fitfolder),'prop')
      delete_wildcard('{}'.format(fitfolder),'property.txt')
      delete_wildcard('{}'.format(fitfolder),'gbw')
      delete_wildcard('{}'.format(fitfolder),'vpot')
   return


def clean_connolly(nconf):
   folder = ['fitcharge_{}'.format(i) for i in range(0,nconf)]
   folder.append('unperturbed')
   for count_i,i in enumerate(folder):
      for q in [-0.5,0.5]:
         configfolder = 'configs_connolly_{}'.format(q)
         delete_wildcard('{}/{}'.format(i,configfolder),'charge.in',addlayer=True)
         delete_wildcard('{}/{}'.format(i,configfolder),'charge.out',addlayer=True)
         delete_wildcard('{}/{}'.format(i,configfolder),'K.tmp',addlayer=True)

   return

def clean_chelpg(nconf):
   folder = ['fitcharge_{}'.format(i) for i in range(0,nconf)]
   folder.append('unperturbed')
   for count_i,i in enumerate(folder):
      for q in [-0.5,0.5]:
         configfolder = 'configs_{}'.format(q)
         delete_wildcard('{}/{}'.format(i,configfolder),'pc_chelpg',addlayer=True)
         delete_wildcard('{}/{}'.format(i,configfolder),'prop',addlayer=True)
         delete_wildcard('{}/{}'.format(i,configfolder),'property.txt',addlayer=True)
         delete_wildcard('{}/{}'.format(i,configfolder),'charge.in',addlayer=True)
         delete_wildcard('{}/{}'.format(i,configfolder),'charge.out',addlayer=True)

   return

def clean_pcsubmit(nconf):
   for count_i,i in enumerate(list(range(0,nconf))):
      fitfolder = 'fitcharge_{}'.format(i)
      delete_wildcard(fitfolder,'pc_layer')

   return

def delete_wildcard(folder,pattern,addlayer=False):
    if addlayer:
      fileList = glob.glob('{}/{}/*/*{}*'.format(run_dir,folder,pattern))
    else:
      fileList = glob.glob('{}/{}/*{}*'.format(run_dir,folder,pattern))
    for filePath in fileList:
       try:
          os.remove(filePath)
       except:
          print("Error while deleting file : ", filePath)

def read_polarout(filename):
    with open(filename,'r') as f:
        flag = 0
        for lc,lines in enumerate(f):
            fields = lines.split()
            if len(fields) < 2 and flag == 0: continue
          
            if flag == 0  and fields[0] == "Isotropic" and fields[1] == "polarizability" :
               polar = float(fields[3])*0.529177210**3 # A.U to A^3
               continue
            
            if flag == 0  and fields[0] == 'DIPOLE' and fields[1] == 'MOMENT' :
               flag = 1
               continue


            if flag == 1 and len(fields) == 4 and fields[0] == 'Magnitude' and fields[1] == '(Debye)':
               dipole = float(fields[3])
               flag = 0 
               continue
    return dipole,polar

def polar_orca(calc_dir,Geometry,Elements):
   os.chdir('{}/{}'.format(run_dir,calc_dir))
   random_factor = 0.0
   with open("polar.in",'w') as f:
      f.write("# calculation of polarizability\n")
      f.write("! wB97X-D3 def2-TZVP TIGHTSCF CHELPG PMODEL PAL4\n")

      f.write("\n%base \"polar\"\n\n")

      f.write("%elprop Polar 1\n")
      f.write("        dipole true\n")
      f.write("        end\n")
      f.write("* xyz {} {}\n".format(0,1))
      for count_i,i in enumerate(Geometry):
         f.write("  {:3s}".format(Elements[count_i]))
         for j in i:
             f.write(" {:< 16.8f}".format((random.random()*2.0-1.0)*random_factor+j))
         f.write("\n")
      f.write("*\n")
   return

def plot_histo(x,title):

    with open(dataname,'r') as f:
       data = {}
       for lc,lines in enumerate(f):
         fields  = lines.split()
         if lc == 0:
            for count_i,i in enumerate(fields):
               if i == 'polar':
                  data["_".join([fields[count_i-1],fields[count_i]])] = []
               else:
                  data[i] = []
         else:
            for count_key,key in enumerate(data):
               if count_key != 0:
                  data[key].append(float(fields[count_key]))
               else:
                  data[key].append(fields[count_key])

    # plot hisotogram
    fig = plt.figure(figsize=(6,4))
    ax = plt.subplot(111)
    n, bins, patches = ax.hist(x, 100, density=True)
    ax.set_xlabel('seperation(A)')
    ax.set_ylabel('Probability')
    # Format ticks and plot box
    ax.tick_params(axis='both', which='major',labelsize=12,pad=10,direction='out',width=2,length=4)
    ax.tick_params(axis='both', which='minor',labelsize=12,pad=10,direction='out',width=2,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]

    # Save the figure
    plt.savefig(title+'_histo', dpi=300, bbox_inches='tight')
    plt.close(fig)

    return

def plot_custom(x,y,name='toy_2',folder='.'):

    # Generate fit plot
    color_list = [(0.05,0.35,0.75),(0.05,0.8,0.6),(0.9,0.3,0.05),(0.35,0.7,0.9),(0.9,0.5,0.7),(0.9,0.6,0.05),(0.95,0.9,0.25),(0.05,0.05,0.05)]*10   # Supposedly more CB-friendly
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plot_handles = []
    #
    #y_ans = 1.5*x+3
    #plot_handel, = ax.plot(x,y,marker='.',markeredgewidth=0.0,linestyle='None',markersize=20,color=color_list[0])
    #plot_handle, = ax.plot(x,y_ans,ls="-",linewidth=4.0,c=color_list[1]) 
    for count_i,i in enumerate(x.keys()):
       plot_handel, = ax.plot(x[i],y[i],marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,alpha=0.3,color=color_list[count_i],label=i)
       #plot_handel, = ax.plot(x,y,marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,alpha=0.3,color=color_list[0])

    y_min,y_max = ax.get_ylim()
    x_min,x_max = ax.get_xlim()
    
    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Set Labels and Save Name
    ax.set_ylabel("xhi2",fontsize=32,labelpad=10,fontweight='bold')
    ax.set_xlabel("iteration",fontsize=32,labelpad=10,fontweight='bold')
    plot_name = folder+'/'+name+'.png'

    # Generate Legend
    ax.legend(loc='best',frameon=False)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1,0.5))

    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=32,pad=10,direction='out',width=4,length=4)
    ax.tick_params(axis='both', which='minor',labelsize=32,pad=10,direction='out',width=4,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]

    # Save the figure
    plt.savefig(plot_name, dpi=300,bbox_extra_artists=(lgd,),bbox_inches='tight')
    plt.close(fig)


    # Save data to file
    dataname = folder+'/'+name+'.txt'
    niter = len(x[list(x.keys())[0]])
    #with open(dataname,'w') as  f:
    #     f.write('{:<15s} {}\n'.format('iter',' '.join([ "{:<15s}".format(i) for i in x ])))
    #     for i in range(niter):
    #        f.write('{:<15d} {}\n'.format(i,' '.join([ "{:< 15.8f}".format(y[j][i]) for j in x ])))
    return



# Wrtie CHARMM fitcharge script
def write_fitinp(c,FinalFF,args,fitfolder,charge_IG):

    ori_dir = os.getcwd()
    os.chdir(run_dir)

    if os.path.isfile(fitfolder+'/fitcharge.inp') is False:
      #sublist = ('geo_opt.xyz  -o {} -q {} -gens {} -npert 0 -ascale {} -mixing_rule wh --fit -inpname fitcharge -FF '.format(fitfolder,c['charge'],c['gens'],args.alpha_scale)).split()
      sublist = ('geo_opt.xyz  -o {}  -gens {} -npert 0 -ascale {} -mixing_rule wh --fit -inpname fitcharge -FF '.format(fitfolder,c['gens'],args.alpha_scale)).split()
      sublist.append('{} {}'.format(charge_IG,FinalFF))
      gen_md_for_CHARMM_new.main(sublist)
    os.chdir(ori_dir)

    return

def place_perturbation(run_dir,calc_folder,xyz):

    ori_dir = os.getcwd()
    os.chdir(run_dir+'/'+calc_folder)

    # Place perturantion charges
    if os.path.isdir('pointcharge_layer1'):
       pc_layer1 = sorted([int(f.split('.')[1]) for f in os.listdir('pointcharge_layer1') if os.path.isfile(os.path.join('pointcharge_layer1', f)) if f.split('.')[1] != 'xyz'])[-1]+1
    else:
       sublist = ('-d 500 -radii_factor 13 -c 1 -extra 0.1 -o layer1 {} -ofolder pointcharge_layer1'.format(xyz)).split()
       pc_layer1 = place_charge.main(sublist)
    if os.path.isdir('pointcharge_layer2'):
       pc_layer2 = sorted([int(f.split('.')[1]) for f in os.listdir('pointcharge_layer2') if os.path.isfile(os.path.join('pointcharge_layer2', f)) if f.split('.')[1] != 'xyz'])[-1]+1
    else:
       sublist = ('-d 500 -radii_factor 22 -c 1 -extra 0.02 -o layer2 {} -ofolder pointcharge_layer2'.format(xyz)).split()
       pc_layer2 = place_charge.main(sublist)
    if os.path.isdir('pointcharge_layer3'):
       pc_layer3 = sorted([int(f.split('.')[1]) for f in os.listdir('pointcharge_layer3') if os.path.isfile(os.path.join('pointcharge_layer3', f)) if f.split('.')[1] != 'xyz'])[-1]+1
    else:
       sublist = ('-d 500 -radii_factor 40 -c 1 -extra 0.003 -o layer3 {} -ofolder pointcharge_layer3'.format(xyz)).split()
       pc_layer3 = place_charge.main(sublist)
   
    total_pc = pc_layer1 + pc_layer2 + pc_layer3

    print("Successfully place {} perturbations charges, layer1: {}, layer2: {}, layer3: {}".format(total_pc,pc_layer1,pc_layer2,pc_layer3))
   
    os.chdir(ori_dir)

    return total_pc,pc_layer1,pc_layer2,pc_layer3

def place_perturbation_submit(calc_folder,xyz,dplus=0):

    ori_dir = os.getcwd()
    os.chdir(run_dir+'/'+calc_folder)

    # Place perturantion charges
    jobids = []
    if os.path.isdir('pointcharge_layer1') is False:
       d = 700+int(dplus)
       sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
       sublist.append('cd {}/{} \n {}/FF_functions/place_charge.py -d {} -radii_factor 13 -c 1 -extra 0.1 -o layer1 {} -ofolder pointcharge_layer1'.format(run_dir,calc_folder,c['taffi_path'],d,xyz))
       sublist= sublist + ('pc_layer1  -p 1 -t {} -q {}'.format(c["param_fit_wt"],c["param_fit_q"])).split()
       output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
       output = str(output,'utf-8')
       jobids += [ output.split("\n")[-2].split()[-1]]
    else:
       if len(os.listdir('pointcharge_layer1')) == 0:
         shutil.rmtree('pointcharge_layer1')
    if os.path.isdir('pointcharge_layer2') is False:
       d = 500 +int(dplus)
       sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
       sublist.append('cd {}/{} \n {}/FF_functions/place_charge.py -d {} -radii_factor 22 -c 1 -extra 0.02 -o layer2 {} -ofolder pointcharge_layer2'.format(run_dir,calc_folder,c['taffi_path'],d,xyz))
       sublist= sublist + ('pc_layer2  -p 1 -t {} -q {}'.format(c["param_fit_wt"],c["param_fit_q"])).split()
       output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
       output = str(output,'utf-8')
       jobids += [ output.split("\n")[-2].split()[-1]]
    else:
       if len(os.listdir('pointcharge_layer2')) == 0:
         shutil.rmtree('pointcharge_layer2')
    if os.path.isdir('pointcharge_layer3') is False:
       d = 500 +int(dplus)
       sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
       sublist.append('cd {}/{} \n {}/FF_functions/place_charge.py -d {} -radii_factor 40 -c 1 -extra 0.003 -o layer3 {} -ofolder pointcharge_layer3'.format(run_dir,calc_folder,c['taffi_path'],d,xyz)) 
       sublist= sublist + ('pc_layer3  -p 1 -t {} -q {}'.format(c["param_fit_wt"],c["param_fit_q"])).split()
       output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
       output = str(output,'utf-8')
       jobids += [ output.split("\n")[-2].split()[-1]]
    else:
       if len(os.listdir('pointcharge_layer3')) == 0:
         shutil.rmtree('pointcharge_layer3')

    print("Submitting place_charge jobs for {}...".format(calc_folder))
   
    os.chdir(ori_dir)

    return jobids

def get_pc(calc_folder):

    ori_dir = os.getcwd()
    os.chdir(run_dir+'/'+calc_folder)

    # Place perturantion charges
    if os.path.isdir('pointcharge_layer1'):
       if len(os.listdir('pointcharge_layer1')) == 0:
         shutil.rmtree('pointcharge_layer1')
         print("ERROR: charge for layer 1 not placed for {}, try increasing density".format(calc_folder))
       else:
         pc_layer1 = sorted([int(f.split('.')[1]) for f in os.listdir('pointcharge_layer1') if os.path.isfile(os.path.join('pointcharge_layer1', f)) if f.split('.')[1] != 'xyz'])[-1]+1
    else:
       print("ERROR: charge for layer 1 not placed for {}".format(calc_folder))
    if os.path.isdir('pointcharge_layer2'):
       if len(os.listdir('pointcharge_layer2')) == 0:
         shutil.rmtree('pointcharge_layer2')
         print("ERROR: charge for layer 2 not placed for {}, try increasing density".format(calc_folder))
       else:
         pc_layer2 = sorted([int(f.split('.')[1]) for f in os.listdir('pointcharge_layer2') if os.path.isfile(os.path.join('pointcharge_layer2', f)) if f.split('.')[1] != 'xyz'])[-1]+1
    else:
       print("ERROR: charge for layer 2 not placed for {}".format(calc_folder))
    if os.path.isdir('pointcharge_layer3'):
       if len(os.listdir('pointcharge_layer3')) == 0:
         shutil.rmtree('pointcharge_layer3')
         print("ERROR: charge for layer 3 not placed for {}, try increasing density".format(calc_folder))
       else:
         pc_layer3 = sorted([int(f.split('.')[1]) for f in os.listdir('pointcharge_layer3') if os.path.isfile(os.path.join('pointcharge_layer3', f)) if f.split('.')[1] != 'xyz'])[-1]+1
    else:
       print("ERROR: charge for layer 3 not placed for {}".format(calc_folder))
   
    total_pc = pc_layer1 + pc_layer2 + pc_layer3

    print("Successfully place {} perturbations charges for {}, layer1: {}, layer2: {}, layer3: {}".format(total_pc,calc_folder,pc_layer1,pc_layer2,pc_layer3))
   
    os.chdir(ori_dir)

    return total_pc,pc_layer1,pc_layer2,pc_layer3



def change_npert(total_pc,fitfolder):

    inpname = '{}/{}/fitcharge.inp'.format(run_dir,fitfolder)
    # Write # of perturnation charges back to CHARMM fit script
    with open(inpname, 'r') as f:
       # read a list of lines into data
       lines = f.readlines()

    # change npert
    for count_i,i in enumerate(lines):
       fields = i.split() 
       if len(fields) < 2: continue
       if fields[0] == 'NPERT':
         lines[count_i] = '    NPERT {} -             ! number of perturbation ions used in QM ESP\n'.format(total_pc)
         break

    # and write everything back
    with open(inpname, 'w') as f:
       f.writelines( lines )

    return

# Calculating CHELPG potential
def calc_CHELPG(calc_folder,xyz,total_pc,pc_layer1,pc_layer2,pc_layer3,geo,element,Qmag,tag):

    ori_dir = os.getcwd()
    os.chdir(run_dir+'/'+calc_folder)
    config_folder = 'configs{}'.format(tag)

    if os.path.isdir(config_folder) is False:
         os.makedirs(config_folder)
    os.chdir("{}/{}/{}".format(run_dir,calc_folder,config_folder))
    for i in range(0,total_pc):
         if os.path.isdir(str(i)) is False: os.makedirs(str(i))
         else:
            scfpfile='{}/{}_charges.scfp'.format(i,i)
            if  os.path.isfile(scfpfile) is True: 
               continue
         os.chdir(str(i))
         #shutil.copyfile(xyz,'{}.xyz'.format(i))
         if i<pc_layer1:
            shutil.copyfile(run_dir+"/{}/pointcharge_layer1/pointcharges.{}.pc".format(calc_folder,i),'pointcharges.pc')
         elif i< pc_layer1+pc_layer2:
            shutil.copyfile(run_dir+"/{}/pointcharge_layer2/pointcharges.{}.pc".format(calc_folder,i-pc_layer1),'pointcharges.pc')
         else:
            shutil.copyfile(run_dir+"/{}/pointcharge_layer3/pointcharges.{}.pc".format(calc_folder,i-(pc_layer1+pc_layer2)),'pointcharges.pc')
   
         # change charge magnitude
         with open('pointcharges.pc', 'r') as f:
            # read a list of lines into data
            lines = f.readlines()
         fields = lines[1].split()
         fields[0] = '{}'.format(str(Qmag)) 
         lines[1] = " ".join(fields)
         with open('pointcharges.pc', 'w') as f:
           f.writelines( lines )
           
         # Write orca input file   
         name = "{}_charges".format(i)
         Write_charge_in(name,geo,element)
         scfpfile='{}_charges.scfp'.format(i)
         os.chdir("{}/{}/{}".format(run_dir,calc_folder,config_folder))
         

    # Create directory for unperturbed charge potential
    if os.path.isdir("unperturbed") is False: os.makedirs("unperturbed") 
    scfpfile='unperturbed/unperturbed_charges.scfp'
    if  os.path.isfile(scfpfile) is False: 
       os.chdir("unperturbed")
       #shutil.copyfile(xyz,'unperturbed.xyz')
       Write_charge_in('unperturbed_charges',geo,element)
       os.chdir("{}/{}/{}".format(run_dir,calc_folder,config_folder))

    # Submit orca job
    jobids = []
    substring = "python3 {}/Automation_Scripts/orca_submit.py  -p {} -t {} -ppn {} -q {} -sched {} -size {} -o charges -path_to_exe {}  --silent"
    substring = substring.format(c["taffi_path"],c["charges_qc_procs"],c["charges_qc_wt"],c["charges_qc_ppn"],c["charges_qc_q"],c["charges_qc_sched"],c["charges_qc_size"],c["orca_exe"])
    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    output = str(output,'utf-8')
    jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
    #monitor_jobs(jobids,user)

    print("COMPLETE: CHELPG potential submission for {}".format(calc_folder))
    os.chdir(ori_dir)

    return jobids

def subproc_string(s):

    s = s.split()

    # Replace , arguments with spaces
    for count_i,i in enumerate(s):
        if "^" in i:
            s[count_i] = s[count_i].replace('^',' ')
    
    return s

def get_connolly_grid(calc_folder,xyz):

    ori_dir = os.getcwd()
    # Genererating connolly surface grid
    os.chdir(run_dir+'/'+calc_folder)
    if os.path.isfile('map.xyz') is False:
         sublist = ('{} -d 1.41'.format(xyz)).split()
         connolly_grid = generate_connolly_new.main(sublist)
    else:
         with open('map.xyz','r') as f:
            for lc,lines in enumerate(f):
               fields = lines.split()
               if lc == 0:
                  connolly_grid = int(fields[0]) 
                  break
    os.chdir(ori_dir)

    return connolly_grid

# Revaluating potential at connolly grid
def calc_connolly_pot(calc_folder,total_pc,connolly_grid,c,tag):

    ori_dir = os.getcwd()
    os.chdir(run_dir+'/'+calc_folder)
    config_folder = 'configs_connolly{}'.format(tag)
    config_ori_folder = 'configs{}'.format(tag)

    if os.path.isdir(config_folder) is False:
         os.makedirs(config_folder)
    os.chdir('{}/{}/{}'.format(run_dir,calc_folder,config_folder))
    for i in range(0,total_pc):
         if os.path.isdir(str(i)) is False: os.makedirs(str(i))
         else: 
            if os.path.isfile(str(i)+'/tmp.output') is True: 
               # write an out file so that bundle_submit will leave it alone
               with open(str(i)+'/charge.out','w') as f:
                  f.write("complete")
               continue
         os.chdir(str(i))
         
         # copy necessary file from configs folder
         #prepare(run_dir,calc_folder,i,tag)
         if os.path.isfile('tmp.output') is False:
            with open('charge.in','w') as f:
               f.write("#!/bin/bash\n")
               f.write("{} {}".format((c['module_string'].split()[-2]).replace('\\n','\n'),(c['module_string'].split()[-1]).replace('\\n','\n')))
               f.write("/depot/bsavoie/apps/orca_4_1_2/orca_vpot {0:}/{1:}/{2:}/{3:}/{3:}_charges.gbw {0:}/{1:}/{2:}/{3:}/{3:}_charges.scfp {0:}/{1:}/map.xyz tmp.output\n".format(run_dir,calc_folder,config_ori_folder,i))
         os.chdir('{}/{}/{}'.format(run_dir,calc_folder,config_folder))
   
    # Submit orca job
    jobids = []
    substring = "python3 {}/Automation_Scripts/bundle_submit.py  -p {} -t {} -ppn {} -q {} -sched {} -size {} -o charges -path_to_exe {}  --silent"
    substring = substring.format(c["taffi_path"],1,c["charges_qc_wt"],c["charges_qc_ppn"],c["charges_qc_q"],c["charges_qc_sched"],c["charges_qc_ppn"],c["orca_exe"])
    output = subprocess.Popen(subproc_string(substring),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
    output = str(output,'utf-8')
    jobids += [ m.split()[-1] for m in output.split("\n")[:-1]]
    #monitor_jobs(jobids,user)
       
    if os.path.isdir('unperturbed') is False: os.makedirs('unperturbed')
    os.chdir('unperturbed')
    i='unperturbed'
    #prepare(run_dir,calc_folder,'unperturbed',tag)
    if os.path.isfile('tmp.output') is False:
       with open('charge.sh','w') as f:
            f.write("#!/bin/bash\n")
            f.write("{} {}".format((c['module_string'].split()[-2]).replace('\\n','\n'),(c['module_string'].split()[-1]).replace('\\n','\n')))
            f.write("/depot/bsavoie/apps/orca_4_1_2/orca_vpot {0:}/{1:}/{2:}/{3:}/{3:}_charges.gbw {0:}/{1:}/{2:}/{3:}/{3:}_charges.scfp {0:}/{1:}/map.xyz tmp.output\n".format(run_dir,calc_folder,config_ori_folder,i))
       command = 'sh charge.sh'.split()
       process = subprocess.Popen(command)
       status = process.poll()
       while status == None:
          status = process.poll()
          time.sleep(3)
    oldpot = '{}/{}/{}/{}/{}_charges.vpot'.format(run_dir,calc_folder,config_ori_folder,i,i)
    newpot = '{}_charges.vpot'.format(i)
    rewrite_potential(oldpot,newpot,connolly_grid)
    print("COMPLETE: revaluation of potential at connolly grid for {}".format(calc_folder))
    os.chdir(run_dir+'/'+calc_folder)
    os.chdir(ori_dir)

    return jobids


def rewrite_pot_wrap(calc_folder,total_pc,connolly_grid,tag):
    config_folder = 'configs_connolly{}'.format(tag)
    config_ori_folder = 'configs{}'.format(tag)
    os.chdir('{}/{}/{}'.format(run_dir,calc_folder,config_folder))
    for i in range(0,total_pc):
         if os.path.isdir(str(i)) is False: os.makedirs(str(i))
         else: 
            if os.path.isfile(str(i)+'_charges.vpot') is True: continue
         os.chdir(str(i))

         oldpot = '{}/{}/{}/{}/{}_charges.vpot'.format(run_dir,calc_folder,config_ori_folder,i,i)
         newpot = '{}_charges.vpot'.format(i)
         rewrite_potential(oldpot,newpot,connolly_grid)
         os.chdir('{}/{}/{}'.format(run_dir,calc_folder,config_folder))
    return

def fit_exe(logname):
    with open('fit.in','w') as f:
         f.write("#!/bin/bash\n")
         f.write("charmm < fitcharge.inp > {}.log\n".format(logname))
    #sublist = ('sh fit.sh').split()
    #output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]

    return

def get_samples(rundir):
    os.chdir(rundir)
    pc_folder = ['pointcharge_layer1','pointcharge_layer2','pointcharge_layer3']
    bond_lp_count = [0,0,0] 
    for count_i,i in enumerate(pc_folder):
      with open('{}/layer{}.xyz'.format(i,count_i+1),'r') as f:
         for lc,lines in enumerate(f):
            fields = lines.split()
            if len(fields) == 0: continue
            if fields[0] == 'F' or fields[0] == 'Cl': bond_lp_count[count_i] += 1
            if fields[0] == 'Br': break
    return bond_lp_count


# calc_folder: folder that has config_connolly
# fit_folder: folder where the fitting is taking place, the file are written in fit_folder/data/ folder
# start_pc, end_pc : control the number of perturbation charges

def write_charmm_data(geo,fitfolder,calc_folder,samples,tag):
  
    config_folder = 'configs_connolly{}'.format(tag)
    config_ori_folder = 'configs{}'.format(tag) 
    skip_grid = 0
    
    data_dir = '{}/{}/data{}'.format(run_dir,fitfolder,tag)
    # Rewrite everything to CHARMM compatible format
    if os.path.isdir(data_dir) is False:
       os.makedirs(data_dir)

    Ang2Bohr = 1.889725989 
    # ESP unit in CHARMM: e/A
    # ESP unit in ORCA: e/Bohr
    # length unit: CHARMM: Angstrom, ORCA: Bohr
    fitpot = "{}/fitcharge.pot".format(data_dir)
    with open(fitpot,'w') as g:
       for i in samples:
          filename = "{}/{}/{}/{}/{}_charges.vpot".format(run_dir,calc_folder,config_folder,i,i) 
          with open(filename,'r') as f:
              atom_count = 0
              pot_count = 0
              for lc,lines in enumerate(f):
                  fields = lines.split()

                  # Initialize lists/arrays
                  if lc == 0: 
                      n_atom    = int(fields[0])
                      n_grid    = int(fields[1])
                      flag = 1
                      continue

                  # Parse the grid and potential
                  #if flag == 2:
                  if flag == 2 and lc > (n_atom+skip_grid):
                      g.write("{:< 10.7e}\n".format(float(fields[0])*Ang2Bohr))
                      pot_count += 1
                      if pot_count == n_grid:
                          break

                  # Parse the molecule
                  if flag == 1:
                      if len(fields) == 0: continue
                      atom_count += 1
                      if atom_count == n_atom:
                          flag = 2
                          continue
      
    fitpot0 = "{}/fitcharge.pot0".format(data_dir)
    with open(fitpot0,'w') as g:
       i = 'unperturbed'
       filename = "{}/{}/{}/{}/{}_charges.vpot".format(run_dir,calc_folder,config_folder,i,i) 
       with open(filename,'r') as f:
           atom_count = 0
           pot_count = 0
           for lc,lines in enumerate(f):
               fields = lines.split()

               # Initialize lists/arrays
               if lc == 0: 
                   n_atom    = int(fields[0])
                   n_grid    = int(fields[1])
                   flag = 1
                   continue

               # Parse the grid and potential
               #if flag == 2:
               if flag == 2 and lc > (n_atom+skip_grid):
                   g.write("{:< 10.7e} {:< 10.7e} {:< 10.7e} {:< 10.7e}\n".format(float(fields[1])/Ang2Bohr,float(fields[2])/Ang2Bohr,float(fields[3])/Ang2Bohr,float(fields[0])*Ang2Bohr))
                   pot_count += 1
                   if pot_count == n_grid:
                       break

               # Parse the molecule
               if flag == 1:
                   if len(fields) == 0: continue
                   atom_count += 1
                   if atom_count == n_atom:
                       flag = 2
                       continue
    qpos = "{}/fitcharge.qpos".format(data_dir)
    with open(qpos,'w') as g:
       for i in samples:
          for count_j,j in enumerate(geo):
              g.write('{:< 15.6f} {:< 15.6f} {:< 15.6f}\n'.format(j[0],j[1],j[2]))
          filename = "{}/{}/{}/{}/pointcharges.pc".format(run_dir,calc_folder,config_ori_folder,i) 
          with open(filename,'r') as f:
              for lc,lines in enumerate(f):
                  fields = lines.split()

                  # Initialize lists/arrays
                  if lc == 1: 
                      g.write('{:< 15.6f} {:< 15.6f} {:< 15.6f}\n'.format(float(fields[1]),float(fields[2]),float(fields[3])))

    
    return
def plot_custom(x,name='toy_2',folder='.'):

    # Generate fit plot
    color_list = [(0.05,0.35,0.75),(0.05,0.8,0.6),(0.9,0.3,0.05),(0.35,0.7,0.9),(0.9,0.5,0.7),(0.9,0.6,0.05),(0.95,0.9,0.25),(0.05,0.05,0.05)]*10   # Supposedly more CB-friendly
    fig = plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plot_handles = []
    xax = range(1,4)
    for count_i,i in enumerate(x[0].keys()):
       charge = [ j[i] for j in x]
       plot_handel, = ax.plot(xax,charge,marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,color=color_list[count_i],label=i)
       #plot_handel, = ax.plot(xax,charge,marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,alpha=0.3,color=color_list[count_i],label=i)
       #plot_handel, = ax.plot(x,y,marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,alpha=0.3,color=color_list[0])

    y_min,y_max = ax.get_ylim()
    x_min,x_max = ax.get_xlim()
    
    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Set Labels and Save Name
    ax.set_ylabel("xhi2",fontsize=32,labelpad=10,fontweight='bold')
    ax.set_xlabel("iteration",fontsize=32,labelpad=10,fontweight='bold')
    plot_name = folder+'/'+name+'.png'

    # Generate Legend
    ax.legend(loc='best',frameon=False)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1,0.5))

    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=32,pad=10,direction='out',width=4,length=4)
    ax.tick_params(axis='both', which='minor',labelsize=32,pad=10,direction='out',width=4,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]

    # Save the figure
    plt.savefig(plot_name, dpi=300,bbox_extra_artists=(lgd,),bbox_inches='tight')
    plt.close(fig)


    # Save data to file
    #dataname = folder+'/'+name+'.txt'
    #niter = len(x[list(x.keys())[0]])
    #with open(dataname,'w') as  f:
    #     f.write('{:<15s} {}\n'.format('iter',' '.join([ "{:<15s}".format(i) for i in x ])))
    #     for i in range(niter):
    #        f.write('{:<15d} {}\n'.format(i,' '.join([ "{:< 15.8f}".format(y[j][i]) for j in x ])))
    return

def plot_custom2(atomtype,dataname,name='toy_2',folder='.'):

    # Generate fit plot
    color_list = [(0.05,0.35,0.75),(0.05,0.8,0.6),(0.9,0.3,0.05),(0.35,0.7,0.9),(0.9,0.5,0.7),(0.9,0.6,0.05),(0.95,0.9,0.25),(0.05,0.05,0.05)]*10   # Supposedly more CB-friendly
    fig = plt.figure(figsize=(6,10))
    ax = plt.subplot(111)
    plot_handles = []

    with open(dataname,'r') as f:
       data = {}
       for lc,lines in enumerate(f):
         fields  = lines.split()
         if lc == 0:
            for count_i,i in enumerate(fields):
               if i == 'polar':
                  data["_".join([fields[count_i-1],fields[count_i]])] = []
               else:
                  data[i] = []
         else:
            for count_key,key in enumerate(data):
               if count_key != 0:
                  data[key].append(float(fields[count_key]))
               else:
                  data[key].append(fields[count_key])
            
    xax = np.array(range(1,len(data['title'])+1))*20
    xax = range(1,4)
    for count_i,i in enumerate(data.keys()):
       if i == 'title': continue
       plot_handel, = ax.plot(xax,data[i],marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,color=color_list[count_i],label=i)
       #plot_handel, = ax.plot(xax,charge,marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,alpha=0.3,color=color_list[count_i],label=i)
       #plot_handel, = ax.plot(x,y,marker='.',markeredgewidth=0.0,linestyle='-',linewidth=4.0,markersize=20,alpha=0.3,color=color_list[0])

    y_min,y_max = ax.get_ylim()
    x_min,x_max = ax.get_xlim()
    
    ax.set_xlim([x_min,x_max])
    ax.set_ylim([y_min,y_max])

    # Set Labels and Save Name
    ax.set_ylabel("values",fontsize=32,labelpad=10,fontweight='bold')
    ax.set_xlabel("# pt",fontsize=32,labelpad=10,fontweight='bold')
    plot_name = folder+'/'+name+'.png'

    # Generate Legend
    ax.legend(loc='best',frameon=False)
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1,0.5))

    # Format ticks
    ax.tick_params(axis='both', which='major',labelsize=20,pad=10,direction='out',width=4,length=4)
    ax.tick_params(axis='both', which='minor',labelsize=20,pad=10,direction='out',width=4,length=3)
    [j.set_linewidth(2) for j in ax.spines.values()]

    # Save the figure
    plt.savefig(plot_name, dpi=300,bbox_extra_artists=(lgd,),bbox_inches='tight')
    plt.close(fig)


    # Save data to file
    #dataname = folder+'/'+name+'.txt'
    #niter = len(x[list(x.keys())[0]])
    #with open(dataname,'w') as  f:
    #     f.write('{:<15s} {}\n'.format('iter',' '.join([ "{:<15s}".format(i) for i in x ])))
    #     for i in range(niter):
    #        f.write('{:<15d} {}\n'.format(i,' '.join([ "{:< 15.8f}".format(y[j][i]) for j in x ])))
    return


def read_db(dbfile):
    with open(dbfile,'r') as f:
      charge_dict = {}
      polar_dict = {}
      for lc,lines in enumerate(f):
         fields = lines.split()
         if len(fields) == 0: continue
         if fields[0] == 'charge':
            charge_dict[fields[1]] = float(fields[2])
         if fields[0] == 'polar':
            polar_dict[fields[1]] = float(fields[2])
    return (charge_dict,polar_dict)
         

def Drude_charge(FinalFF,localFF,args,c,fitfolder):
    
    ori_dir = os.getcwd()
    os.chdir(run_dir+'/'+fitfolder)
    print("Start parametrization....") 
    if os.path.isfile('fitcharge.log') is False:     
       jobids = []
       sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
       sublist.append('module load openmpi/3.1.4 \n/depot/bsavoie/apps/charmm/exec/gnu/charmm < fitcharge.inp > fitcharge.log')
       sublist= sublist + ('charmm_fit  -p 1 -t {} -q {}'.format(c["param_fit_wt"],c["param_fit_q"])).split()
       output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
       output = str(output,'utf-8')
       jobids += [ output.split("\n")[-2].split()[-1]]
       monitor_jobs(jobids,user)
    if os.path.isfile('fitcharge.log') is False:
       print("ERROR: fitcharge.log not found, something went wrong with executing CHARMM")
       quit()
    else:
       print("COMPLETE: CHARMM parametrization")

    # Converting log file to db file   
    command = ('-o {}'.format(args.outFF)).split()
    log_to_db.main(command)

    os.chdir(run_dir)
    # Merge these Params with the parent
    # XXX can't merge with parent, if the other inchi is also merging then would have error
    option = {}
    #option['master'] = FinalFF
    option['new_params'] = './{}/{}'.format(fitfolder,args.outFF)
    #option['no_ct'] = True
    #option['replace'] = 1
    # Merge these with the local (total_drude.db)
    option['master'] = localFF 
    option['no_ct'] = False
    merge_FF_drude.fun(option)
    print("COMPLETE: Drude charge parametrization")
   
    os.chdir(ori_dir)

    return

def VDW_re(element,geo,FinalFF,localFF,c,args,tag,fitfolder):

    ori_dir = os.getcwd()
    print("Start vdw reparametrization...")
    # move to non_polar folder to relax Drude particles
    nonpolar_dir =  '/home/lin1209/fileserver/halstead/benchmark'
    inchikey = GetInchi(element,geo)
    vdwfolder = nonpolar_dir + '/' + inchikey + '/vdw' 
    if os.path.isdir(vdwfolder) is False:
      print("ERROR: nonpolar vdw folder {} doesn't exist. Exiting...".format(vdw))
    os.chdir(vdwfolder)
    print("move to {} to reparametrize...".format(vdwfolder))
    # converting the new db file to charmm format FF which is needed when relaxing via CHARMM
    if tag == '': FFname = 'FF'
    else: FFname = 'FF_{}'.format(tag)
    if os.path.isdir(vdwfolder+'/'+FFname) is False:
      #sublist = ('{}/geo_opt.xyz  -o FF -q {} -gens {} -mixing_rule wh --polar -FF '.format(run_dir,c['charge'],c['gens'])).split()
      sublist = ('{}/geo_opt.xyz  -o {}  -gens {} -mixing_rule wh --polar -FF '.format(run_dir,FFname,c['gens'])).split()
      sublist.append('{} {}'.format(FinalFF,localFF))
      #gen_md_for_CHARMM.main(sublist)
      db_to_charmm.main(sublist)
    os.chdir('configs')
    if tag == '':
      relax_drude_loop.main([])
    else:
      relax_drude_loop.main('-tag {}'.format(tag).split())

    # initial_param.db serves as extract_vdw's initial guess
    # normally this would be the FF from the cycle before the very last cycle
    # we want to use the fully parametrized FF as our initial guess
    os.chdir('..')
    shutil.copyfile('cycle-2/DFT-AA.db','initial_params.db')
    
    # mv to inchikey folder to execute extract_vdw
    os.chdir('..')

    L2_s = "0.1" 
    L2_e = "0.5"
    jobids = []
    if tag == '':
      vdwdrude_dir = 'vdw_drude' 
    else:
      vdwdrude_dir = 'vdw_drude_{}'.format(tag)
    if os.path.isdir('vdw/'+vdwdrude_dir) is False:
       # we can't use the final db as input or no parametrization would be done
       # so we used FF input of the original parse but switch the charge to Drude charge
       # switching is just precaution, the charges in extract_vdw_drude are from xyz file which relax_drude would write
       FF_path = FinalFF
       charge_FF = run_dir + '/'+fitfolder+'/'+args.outFF 
       sublist = ['{}/Automation_Scripts/shell_submit.sh'.format(c['taffi_path'])]
       if tag == '':
            sublist.append('{}/FF_functions/extract_vdw_drude.py -f vdw -FF_DFT \'{} {}\' -o {} -E_max 0.0 -xhi2_thresh 1E-8 -q_a 0.0 -q_b 0.0 -mixing_rule wh -L2_sigma {} -L2_eps {} -fun {} '.format(c['taffi_path'],FF_path,charge_FF,vdwdrude_dir,L2_s,L2_e,c['functional']))
       else:
            sublist.append('{}/FF_functions/extract_vdw_drude.py -f vdw -FF_DFT \'{} {}\' -o {} -E_max 0.0 -xhi2_thresh 1E-8 -q_a 0.0 -q_b 0.0 -mixing_rule wh -L2_sigma {} -L2_eps {} -fun {} -tag {}'.format(c['taffi_path'],FF_path,charge_FF,vdwdrude_dir,L2_s,L2_e,c['functional'],tag))
         
       sublist= sublist + ('vdw_drude_parse  -p 1 -t {} -q {}'.format(c["param_fit_wt"],c["param_fit_q"])).split()
       output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
       output = str(output,'utf-8')
       jobids += [ output.split("\n")[-2].split()[-1]]
       # Wait until jobs complete
       monitor_jobs(jobids,user)
    if os.path.isfile('vdw/'+vdwdrude_dir+'/DFT-AA.db') is False:
       print("ERROR: vdw/{}/DFT-AA.db not found, something went wrong w/ vdw_parse_drude".format(vdwdrude_dir))
       quit()
    else:
       print("COMPLETE: reparametrize vdw")
    
    # Merge these Params with the Parent
    # XXX can't merge with parent, if the other inchi is also merging then would have error
    option = {}
    option['new_params'] = './vdw/{}/DFT-AA.db'.format(vdwdrude_dir)
    #option['replace'] = 1
    # Merge these Params with the local(total_drude.db)
    option['master'] = localFF
    option['no_ct'] = False
    merge_FF_drude.fun(option)

    os.chdir(ori_dir)

    return
    

def write_datadir(tag):
    with open(run_dir+'/datadir{}.def'.format(tag),'w') as f:
       f.write("* CHARMM testcase data directory assignment\n"+\
               "*\n"+\
               "faster on\n"+\
               "set pnode =\n"+\
               "if ?numnode .gt. 1 set pnode = node 0\n"+\
               "set 0 data{}/     ! input data directory\n".format(tag)+\
               "set 9 scratch/  ! scratch directory\n"+\
               "set testcheck stream @0/test.str\n"+\
               "set qcheck stream @0/qtest.str\n"+
               "set testfail 0\n"+\
               "return\n")
    return


# cp .gbw, .scfp, .xyz file from calc_folder to current folder
def prepare(calc_folder,name,tag):
   config_folder = 'configs{}'.format(tag)
   i = name
   shutil.copyfile('{}/{}/{}/{}/{}_charges.gbw'.format(run_dir,calc_folder,config_folder,i,i),'{}_charges.gbw'.format(i))
   shutil.copyfile('{}/{}/{}/{}/{}_charges.scfp'.format(run_dir,calc_folder,config_folder,i,i),'{}_charges.scfp'.format(i))
   shutil.copyfile('{}/{}/{}/{}/{}.xyz'.format(run_dir,calc_folder,config_folder,i,i),'{}.xyz'.format(i))
   if i != 'unperturbed':
      shutil.copyfile('{}/{}/{}/{}/pointcharges.pc'.format(run_dir,calc_folder,config_folder,i),'pointcharges.pc')
   shutil.copyfile('{}/{}/map.xyz'.format(run_dir,calc_folder),'map')

   return
   

# This function grep FF_path from vdw_parse.submit
# it can also change the original folder path
# this can be handy if the original folder has been moved 
# or when the parametrization happened in a different folder
# it also change total.db to total_drude.db
# This function is only needed right now because the generation structure hasn't been coded

def get_FF_path(filename,changefolder=''):
    with open(filename,'r') as f:
      for lc,lines in enumerate(f):
         fields = lines.split() 
         if len(fields) < 5: continue
         if fields[3] == '-FF_DFT':
            newfields = lines.split('\'')
            FF = newfields[1].split()[1:]
            for count_i,i in enumerate(FF):
               singleFF = i.split('/')
               if singleFF[-1] == 'total.db':
                  singleFF[-1] = 'total_drude.db'
                  FF[count_i] = '/'.join(singleFF)
            FF = ' '.join(FF)
    return FF
            

def Write_charge_in(name,geo,element):
   with open('charge.in','w') as f:
      procs = 8
      D3_option = "D3BJ "
      if(c['functional'] == 'wB97X-D3'):
         D3_option=""
      f.write('#Run DFT-D3 single point of the dimer\n')
      f.write('! {} {} TIGHTSCF {} CHELPG PAL{} KeepDens\n\n'.format(c['functional'],c['basis'],D3_option,procs))
      f.write('%base "{}"\n\n'.format(name))
      if name != 'unperturbed_charges':
         f.write('% pointcharges "pointcharges.pc"\n\n')
      # Increase default max scf iterations
      f.write('\n\n%scf\nMaxIter 1000\nend\n')

      # Add print statements for different charge definitions
      f.write('%output\n')
      f.write('  {}\n  {}\n  {}\n  {}\n  {}\nend\n'.format('Print[ P_Mayer ] 1','Print[ P_NatPop ] 1','Print[ P_Hirshfeld ] 1', 'Print[ P_Mulliken ] 1','Print [ P_Loewdin ] 1'))
      # Write Coordinate header
      f.write('\n* xyz 0 1\n')
      # Write Geometry 
      for count_j,j in enumerate(geo):
          f.write('  {:<20s}   {:< 20.6f} {:< 20.6f} {:< 20.6f}\n'.format(element[count_j],j[0],j[1],j[2]))
      f.write('*\n')
   return

def read_alljson(jsonfile):
    if os.path.isfile(jsonfile) is False:
      print("Error: json file: {} not found".format(jsonfile))
      quit()
    obj_text = codecs.open(jsonfile, 'r').read()
    jload = json.loads(obj_text)
    return jload


def GetInchi(Elements,Geometry): #convert geo2inchi key

    # Open file for writing and write header
    # create a xyz file for openbabel input
    fid = open('inchi.xyz','w')
    fid.write('{}\n\n'.format(len(Elements)))

    for count_i,i in enumerate(Elements):
        fid.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} \n'.format(i,Geometry[count_i,0],Geometry[count_i,1],Geometry[count_i,2]))

    fid.close()
    # create inchikey
    p = subprocess.Popen("obabel -ixyz inchi.xyz -oinchikey".split(),stdin=PIPE, stdout=PIPE, stderr=PIPE)
    inchikey, err = p.communicate()
    inchikey = str(inchikey, 'utf-8').split('\n')[0] ##might be a problem if there are characters utf-8 can't handle
    # remove xyz file
    os.remove("inchi.xyz")
   
    return inchikey
"""
def Get_Inchinew(Elements,Geometry):

    # Open file for writing and write header
    # create a xyz file for openbabel input
    fid = open('inchi.xyz','w')
    fid.write('{}\n\n'.format(len(Elements)))

    for count_i,i in enumerate(Elements):
        fid.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} \n'.format(i,Geometry[count_i,0],Geometry[count_i,1],Geometry[count_i,2]))

    fid.close()
    import openbabel
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "inchikey")

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, "inchi.xyz")

    inchikey=obConversion.WriteString(mol).rstrip()
    # remove xyz file
    os.remove("inchi.xyz")

    return inchikey

def Get_Inchi_xyz(xyz):

    import openbabel
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "inchikey")

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, xyz)

    inchikey=obConversion.WriteString(mol).rstrip()

    return inchikey
"""


# Transform the output potential from orca_vpot to chelpg output potential
# needed the chelpg output for mol's xyz 
#    tmp.output(output from orca_vpot) format: 
#    gridx gridy gridz pot
#    gridx gridy gridz pot
#    ....
#    *_charges.vpot(original from chelpg) format:
#    natoms total_grid_pt
#         xyz of mol
#         ...
#         ...
#    pot  gridx gridy gridz
#    pot  gridx gridy gridz
#         ....
def rewrite_potential(oldpot,newpot,connolly_grid):
   with open(oldpot,'r') as f:
      with open(newpot,'w') as g:
        atom_count = 0
        pot_count = 0
        for lc,lines in enumerate(f):
            fields = lines.split()

            # Initialize lists/arrays
            if lc == 0: 
                n_atom    = int(fields[0])
                cosmo_grid    = int(fields[1])
                g.write("{} {}\n".format(n_atom,connolly_grid))
                flag = 1
                continue

            # Parse the molecule
            if flag == 1:
                if len(fields) == 0: continue
                g.write(lines)
                atom_count += 1
                if atom_count == n_atom:
                    flag = 2
                    break

   with open('tmp.output','r') as f:
        with open(newpot,'a') as g:
           for lc,lines in enumerate(f):
               fields = lines.split()
            # Initialize lists/arrays
               if lc == 0: 
                   if connolly_grid != int(fields[0]):
                        print("ERROR: QM output connolly grid number: {} doesn't match map.xyz:{} ".format(fields[0],connolly_grid))
                        quit()
                   flag = 1
                   continue
               if flag == 1 :
                  g.write("{:< 10.7e}     {:< 10.7e}  {:< 10.7e}  {:< 10.7e}\n".format(float(fields[3]),float(fields[0]),float(fields[1]),float(fields[2])))
   return


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
