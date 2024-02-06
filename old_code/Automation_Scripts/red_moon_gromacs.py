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
    run_dir = current_dir

    """
    # Read all.json (contains model compound and dependency information for the batch 
    status = read_alljson('./all.json')
    # Find inchi keys for this generation
    mc_keys = [ _ for _ in status["mc"].keys() ]
    # Find inchi keys for this generation
    #mc_keys = [ _ for _ in status["mc"].keys() if int(g) in status["mc"][_]["gen"] ]
    gens = sorted(status["gens"].keys())
    with open('tmp_bundle.sh','w') as f:
      f.write("#!/bin/bash\n")
      for g in gens:
         keys = [ _ for _ in status["mc"].keys() if int(g) in status["mc"][_]["gen"] ] 
         f.write("gen{}=(".format(g))
         for i in keys:
            f.write('\'{}\' '.format(i))
         f.write(')\n')
    #for i in mc_keys:
    #  os.makedirs(i)
    quit()
    """
    parser = argparse.ArgumentParser(description='Driver script for taffi benchmark simulations. Expects a configuation file formatted with space delimited keyword arguments')

    #optional arguments                                                                                                                                                        
    parser.add_argument('--restart', dest='restart_flag', default=False, action='store_const', const=True,
                        help = 'When present, will search for the lasted cycle then pick right up') 

    args=parser.parse_args()    
    sys.stdout = Logger(current_dir+'/red_moon.log')
    run_red_moon(restart=args.restart_flag)
    quit()

    
    

    return

def run_red_moon(restart=False):
    # activiation energy for each rxn
    cycle = 0
    rate_const = [2.27e11,3.39e10] # rate constant 
    pre_exp = [1.673,1] # rate constant 
    ER = [ 41.23,44.16] #kcal/mol
    dUo = [ (41.23-44.16), (44.16-41.23)] # E_forward - E_backward
    resi = ['HYD','HYI','IOD']
    kb=0.001987204118                           # Boltzman constant in units of kcal/mol K
    T = 300
    beta = 1/(kb*T)
    box = 250
    tot_cycle = 6000

    if restart:
       line = get_last_line('conc.txt') 
       start_c = int(line.split()[0])
       print("restart from {}...".format(start_c))
       """
       logfile = []
       wild_card_files = [ '*.log']
       for i in wild_card_files:
            path = '/'.join(i.split('/')[:-1])
            if len(path) == 0:
                path = "."
            logfile += [ path+"/"+files for files in os.listdir(path) if fnmatch.fnmatch(files,i.split('/')[-1]) ] 
       logfile=sorted(logfile, key = log_num)   
       start_c = int(logfile[-1].split('_')[0].split('/')[-1])+1
       """
    else:
       atom,id_atom = read_gro('0.gro')
       nmol = get_nmol(id_atom,atom,resi)  
       f = open('traj.pdb','w')
       save_pdb(atom,box,0,f) # save the initial state 
       f.close()
       g = open('conc.txt','w')
       g.write("{:<6s} {:<6s} {:<6s} {:<6s}\n".format('cycle',resi[0],resi[1],resi[2]))
       g.write("{:<6d} {:<6d} {:<6d} {:<6d} {:<6.3f} {:<6.3f} {:<6.3f}\n".format(0,nmol[0],nmol[1],nmol[2],nmol[0]/125,nmol[1]/125,nmol[2]/125))
       g.close()
       start_c = 0
    cycle = start_c
    reject_flag = False
    pair_fail = []
    while cycle < tot_cycle:
      print("Starting cycle: {}".format(cycle))
      atom,id_atom = read_gro('{}.gro'.format(cycle))
      nmol = get_nmol(id_atom,atom,resi)  
      nmol_ori = deepcopy(nmol)
      cand = find_cand()
      os.rename("trj.gro", "prev.gro")
      # if the previous pair rxn is, then we would have the same configuration
      # we want to avoid picking the same pair so that the rxn wouldn't get stucked
      print("failing pair: {}".format(pair_fail))
      if reject_flag:
         for j in pair_fail:
            for rxn_count,rxn_pair in enumerate(cand):
               for count_i,i in enumerate(rxn_pair):
                  if i == j  or (i[1],i[0]) == j:
                     del cand[rxn_count][count_i]
      retry = 0
      while np.amax([len(i) for i in cand]) == 0:
         print("no candidate found, executing NVT...")
         write_top(nmol_ori,resi)
         run_nvt(resi,cycle,T,nmol,box,switch_flag=False,rseed=random.randint(0,10000),min_flag=False)
         cycle += 1
         atom,id_atom = read_gro('{}.gro'.format(cycle))
         f = open('traj.pdb','a')
         save_pdb(atom,box,cycle,f) 
         f.close()
         g = open('conc.txt','a') 
         g.write("{:<6d} {:<6d} {:<6d} {:<6d} {:<6.3f} {:<6.3f} {:<6.3f}\n".format(cycle,nmol_ori[0],nmol_ori[1],nmol_ori[2],nmol_ori[0]/125,nmol_ori[1]/125,nmol_ori[2]/125))
         g.close()
         cand = find_cand()
         os.rename("trj.gro", "prev.gro")
         retry += 1
         pair_fail = []
      print("retry: {}".format(retry))
      # choose rxn
      print(cand)
      weight = [ len(cand[i])*np.exp(-beta*ER[i]) for i in range(0,len(cand))]
      weight[0] = weight[0]*(22.7/3.39)  # preexponential ratio
      #weight[0] = weight[0]*1.673  # preexponential ratio
      #weight = [ len(cand[i])*rate_const[i] for i in range(0,len(cand))]
      [rxn] = random.choices([0,1], weights=weight, k=1) 
      print("rxn: {} chosen!".format(rxn))
      # choose pair(s)
      [pair] = random.sample(range(0,len(cand[rxn])),1)
      pair = cand[rxn][pair]
      # switch pot for chosen pair( change resi name) then write out to n_switch.crd
      switch_pot(atom,id_atom,pair,cycle,rxn,np.sum(np.array(nmol)),box)
      # relax switch structure
      atom,id_atom = read_gro('{}_switch.gro'.format(cycle))
      nmol = get_nmol(id_atom,atom,resi)
      nmol_switch = deepcopy(nmol)
      write_top(nmol_switch,resi)
      run_nvt(resi,cycle,T,nmol,box,switch_flag=True,rseed=random.randint(0,10000),min_flag=True)
      ene_r = get_last_ene('{}.log'.format(cycle))
      ene_s = get_last_ene('{}_switch.log'.format(cycle))
      criteria = np.exp(-beta*(ene_s['E_pot']-ene_r['E_pot']+dUo[rxn]))
      if random.random() > criteria: # reject
        print("P:{}".format(criteria))
        print("delU_MM: {}".format(ene_s['E_pot']-ene_r['E_pot']))
        print("delUo: {}".format(dUo[rxn]))
        print("rxn reject!")
        reject_flag = True
        pair_fail.append(pair)
        f = open('traj.pdb','a')
        atom,id_atom = read_gro('{}.gro'.format(cycle))
        save_pdb(atom,box,cycle+1,f) 
        f.close()
        g = open('conc.txt','a') 
        g.write("{:<6d} {:<6d} {:<6d} {:<6d} {:<6.3f} {:<6.3f} {:<6.3f}\n".format(cycle+1,nmol_ori[0],nmol_ori[1],nmol_ori[2],nmol_ori[0]/125,nmol_ori[1]/125,nmol_ori[2]/125))
        g.close()
        shutil.copyfile('{}.gro'.format(cycle),'{}.gro'.format(cycle+1))
        shutil.copyfile('{}.log'.format(cycle),'{}.log'.format(cycle+1))
        os.rename("prev.gro", "trj.gro")
      else: # accept
        print("delU_MM: {}".format(ene_s['E_pot']-ene_r['E_pot']))
        print("delUo: {}".format(dUo[rxn]))
        print("rxn accept!")
        reject_flag = False
        pair_fail = []
        f = open('traj.pdb','a')
        atom,id_atom = read_gro('{}_switch.gro'.format(cycle))
        save_pdb(atom,box,cycle+1,f) 
        f.close()
        g = open('conc.txt','a') 
        g.write("{:<6d} {:<6d} {:<6d} {:<6d} {:<6.3f} {:<6.3f} {:<6.3f}\n".format(cycle+1,nmol_switch[0],nmol_switch[1],nmol_switch[2],nmol_switch[0]/125,nmol_switch[1]/125,nmol_switch[2]/125))
        g.close()
        shutil.copyfile('{}_switch.gro'.format(cycle),'{}.gro'.format(cycle+1))
        shutil.copyfile('{}_switch.log'.format(cycle),'{}.log'.format(cycle+1))
      delete_wildcard('','#')
      cycle += 1

    return
def write_top(nmol,resi):
    with open('topol.top','w') as f:
      f.write("#include \"all.itp\"\n\n")

      f.write("[system]\n")
      f.write("condense\n\n")

      f.write("[molecules]\n")
      for count_i,i in enumerate(resi):
         f.write("{} {}\n".format(i,nmol[count_i]))
    return

def log_num(x):
    return(int(x.split('_')[0].split('/')[-1]))

def get_last_line(filename):
   with open(filename, 'rb') as f:
       f.seek(-2, os.SEEK_END)
       while f.read(1) != b'\n':
           f.seek(-2, os.SEEK_CUR)
       last_line = f.readline().decode()
   return last_line

def get_last_ene(logfile):
 
   kcal2kJ=4.184
   """
   write_enesh(title) 
   subprocess.call("chmod 777 ener.sh", shell=True)
   sublist = ('sh ener.sh').split()
   output = subprocess.Popen(sublist)
   status = output.poll()
   while status == None:
      status = output.poll()
      time.sleep(0.5)
   with open('energy.xvg','r') as f:
      for line in f:
         pass
      last_line = line
      fields = last_line.split()
      ene = {}
      print(fields[1])
      ene['E_pot'] = float(fields[1])/kcal2kJ
   """
   with open(logfile,'r') as f:
      flag = 0
      for lc,lines in enumerate(f):
         fields = lines.split()
         if len(fields) == 0: continue
         if fields[0] == 'nsteps': 
            final_step = int(fields[2])
            continue
         if flag == 0 and fields[0] == 'Step' and fields[1] == 'Time': 
            flag =1
            continue
         if flag == 1:
            current_step = int(fields[0])
            if current_step == final_step: flag = 2
            else: flag = 0
            continue
         if flag == 2:
            for count_i,i in enumerate(fields):
               if i == 'Potential':
                  flag = 3
                  pot_col = count_i
                  break
            continue
         if flag == 3:
            ene = {}
            if pot_col == 8: pot_col = -1
            ene['E_pot'] = float(fields[pot_col])/kcal2kJ
            flag = 0
         
   return ene
         
   
def write_enesh(title):
   with open('ener.sh','w') as f:
      f.write("#!/bin/bash\n\n") 

      f.write("module load gcc/8.3.0/\n")
      f.write("module load openmpi/3.1.4\n\n")
      f.write("GRO=/home/lin1209/install/gromacs_2/bin/gmx\n")

      f.write("title={}\n".format(title))
      f.write("$GRO energy -f ${title}.edr <<EOF\n")
      f.write("potential\n")
      f.write("EOF\n")
   return

def save_pdb(atom,box,cycle,f):
   f.write("ITEM: TIMESTEP\n")
   f.write("{}\n".format(cycle))
   f.write("ITEM: NUMBER OF ATOMS\n")
   f.write("{}\n".format(len(list(atom.keys()))))
   f.write("ITEM: BOX BOUNDS pp pp pp\n")
   f.write("{:4.3f} {:4.3f}\n".format(-box/2,box/2))
   f.write("{:4.3f} {:4.3f}\n".format(-box/2,box/2))
   f.write("{:4.3f} {:4.3f}\n".format(-box/2,box/2))
   f.write("ITEM: ATOMS id type xs ys zs vx vy vz\n")
   formatting = "{:<5d} {:<5d} {:<11.7f} {:<11.7f} {:<11.7f} {:<5d} {:<5d} {:<5d}\n"
   for a in atom:
      if atom[a]['resi'] == 'HYD': vx = 1
      if atom[a]['resi'] == 'HYI': vx = 2
      if atom[a]['resi'] == 'IOD': vx = 3
      if atom[a]['elem'] == 'H': vy = 1
      if atom[a]['elem'] == 'I': vy = 2
      line = formatting.format(a+1,atom[a]['mol_id']+1,atom[a]['coord'][0]/box,atom[a]['coord'][1]/box,atom[a]['coord'][2]/box,vx,vy,0)
      f.write(line)
   
   """
   f.write("REMARK from red_moon.py\n")
   f.write("TITLE     Ethanol in water t= {:10.5f} step= {:}\n".format(cycle,cycle)) 
   f.write("REMARK    THIS IS A SIMULATION BOX\n")
   f.write("CRYST1   {0:<6.3f}   {0:<6.3f}   {0:<6.3f}  90.00  90.00  90.00 P 1           1 \n".format(box))
   f.write("MODEL        1\n")
   for a in atom:
      f.write("{:<4s}  {:>5d}  {:<4s} {:>3s} {:>4d}  {:>8.3f} {:>8.3f} {:>8.3f} {:>6.2f} {:>6.2f}\n".format('ATOM',a+1,atom[a]['atom_name'],atom[a]['resi'],atom[a]['mol_id']+1,atom[a]['coord'][0],atom[a]['coord'][1],atom[a]['coord'][2],1,0))
   f.write("TER \n")
   f.write("ENDMDL\n")
   """

   return




def rewrite_test(atom):
   box = 250
   g = open('test.gro','w')
   time = 0
   g.write("title line, t= {}\n".format(time))
   g.write(" {}\n".format(len(list(atom.keys()))))
   formatting = "{:5d}{:<5s}{:>5s}{:>5d}{:> 8.3f}{:> 8.3f}{:> 8.3f}\n"
   for a in atom:
      line = formatting.format(atom[a]['mol_id']+1,atom[a]['resi'],atom[a]['atom_name'],a+1,atom[a]['coord'][0]/10,atom[a]['coord'][1]/10,atom[a]['coord'][2]/10,0,0,0)
      g.write(line)
   g.write("{0:<8.3f} {0:<8.3f} {0:<8.3f}\n".format(box/10))

         
   

def opt_structure(choice):
   mol = {}
   mol['HYD'] = [np.array([-9.68307,0.24769,0.0]),np.array([-9.08483,-0.13096,0.0])]
   mol['IOD'] = [np.array([-8.02825,2.05576,0.0]),np.array([-6.18885,0.14534,0.0])]
   mol['HYI'] = [np.array([-16.38615,3.98231,-10.54242]),np.array([-14.73423,3.98231,-10.54242])]

   dist = np.linalg.norm(mol[choice][1]-mol[choice][0])

   return dist

def switch_pot(atom,id_atom,pair,cycle,rxn,tot_mol,box):
   formatting = "{:5d}{:<5s}{:>5s}{:>5d}{:> 8.3f}{:> 8.3f}{:> 8.3f}\n"
   with open('{}_switch.gro'.format(cycle),'w') as f:
      f.write('title\n')
      f.write(' 500\n')
       
      # change label for chosen pair
      # rxn0 : H2+I2 -> 2HI
      # rxn1 : 2HI -> H2+I2
      # atoms that are bonded to the chosen pair atoms
      id_atom[tot_mol] = []
      id_atom[tot_mol+1] = []
      for count_p,p in enumerate(pair):
         # find another atom
         for j in id_atom[atom[p]['mol_id']]:
            if j != p : 
               a = j
         id_atom.pop(atom[p]['mol_id'],None)
         id_atom[tot_mol+count_p].append(p)
         id_atom[tot_mol+count_p].append(a)
         atom[p]['mol_id'] = tot_mol 
         # rxn 0: switch H2 -> HI and I2-> HI
         if atom[p]['resi'] == 'HYD':
            d = unwrap_d(atom[a]['coord']-atom[p]['coord'],box)
            d = d*opt_structure('HYI')/(np.linalg.norm(d))
            atom[a]['coord'] = atom[p]['coord'] + d
            atom[p]['resi'] = 'HYI'
            atom[a]['resi'] = 'HYI'
            atom[p]['atom_name'] = 'H1'
            atom[a]['atom_name'] = 'I1'
            atom[a]['elem'] = 'I'
            continue
         if atom[p]['resi'] == 'IOD':
            d = unwrap_d(atom[a]['coord']-atom[p]['coord'],box)
            d = d*opt_structure('HYI')/(np.linalg.norm(d))
            atom[a]['coord'] = atom[p]['coord'] + d
            atom[p]['resi'] = 'HYI'
            atom[a]['resi'] = 'HYI'
            atom[p]['atom_name'] = 'I1'
            atom[a]['atom_name'] = 'H1'
            atom[a]['elem'] = 'H'
            continue
         if atom[p]['resi'] == 'HYI':
            if atom[p]['elem'] == 'H':
               d = unwrap_d(atom[a]['coord']-atom[p]['coord'],box)
               d = d*opt_structure('HYD')/(np.linalg.norm(d))
               atom[a]['coord'] = atom[p]['coord'] + d
               atom[p]['resi'] = 'HYD'
               atom[a]['resi'] = 'HYD'
               atom[p]['atom_name'] = 'H1'
               atom[a]['atom_name'] = 'H2'
               atom[a]['elem'] = 'H'
            if atom[p]['elem'] == 'I':
               d = unwrap_d(atom[a]['coord']-atom[p]['coord'],box)
               d = d*opt_structure('IOD')/(np.linalg.norm(d))
               atom[a]['coord'] = atom[p]['coord'] + d
               atom[p]['resi'] = 'IOD'
               atom[a]['resi'] = 'IOD'
               atom[p]['atom_name'] = 'I1'
               atom[a]['atom_name'] = 'I2'
               atom[a]['elem'] = 'I'
            continue
      id_resi = [ (atom[id_atom[mol_id][0]]['resi'],mol_id) for mol_id in id_atom]
      id_resi.sort(key = lambda x: x[0])  

      atom_count = 0
      for resi_count,(resi,mol_id) in enumerate(id_resi):
         for a in id_atom[mol_id]:
            line = formatting.format(resi_count+1,atom[a]['resi'],atom[a]['atom_name'],atom_count+1,atom[a]['coord'][0]/10,atom[a]['coord'][1]/10,atom[a]['coord'][2]/10)
            f.write(line)
            atom_count += 1
      f.write("{0:<8.3f} {0:<8.3f} {0:<8.3f}\n".format(box/10))
   
   return 0

def switch_pot_old(atom,id_atom,pair,cycle,rxn,tot_mol,box):
   formatting = "{:5d}{:<5s}{:>5s}{:>5d}{:> 8.3f}{:> 8.3f}{:> 8.3f}\n"
   with open('{}_switch.gro'.format(cycle),'w') as f:
      f.write('title\n')
      f.write(' 500\n')
       
      # change label for chosen pair
      # rxn0 : H2+I2 -> 2HI
      # rxn1 : 2HI -> H2+I2
      # atoms that are bonded to the chosen pair atoms
      other_p = [ a for p in pair for a in id_atom[atom[p]['mol_id']] if a != p] 
      id_atom[tot_mol] = []
      id_atom[tot_mol+1] = []
      if rxn == 0:
         for count_p,p in enumerate(pair):
            id_atom.pop(atom[p]['mol_id'],None)
            id_atom[tot_mol].append(p)
            atom[p]['mol_id'] = tot_mol 
            atom[p]['resi'] = 'HYI'
            if atom[p]['elem'] == 'H': atom[p]['atom_name'] = 'H1'
            if atom[p]['elem'] == 'I': 
               atom[p]['atom_name'] = 'I1'
               d = unwrap_d(atom[pair[1-count_p]]['coord']-atom[p]['coord'],box)
               d = d*1.639/(np.linalg.norm(d))
               atom[pair[1-count_p]]['coord'] = atom[p]['coord'] + d
               #atom[pair[1-count_p]]['coord'] = wrap_coord(atom[p]['coord'],box)
         for count_p,p in enumerate(other_p):
            id_atom[tot_mol+1].append(p)
            atom[p]['mol_id'] = tot_mol+1
            atom[p]['resi'] = 'HYI'
            if atom[p]['elem'] == 'H': atom[p]['atom_name'] = 'H1'
            if atom[p]['elem'] == 'I': 
               atom[p]['atom_name'] = 'I1'
               atom[p]['coord'] = atom[other_p[1-count_p]]['coord'] + np.array([0.946,0.946,0.946])
               #atom[p]['coord'] = wrap_coord(atom[p]['coord'],box)
      if rxn == 1:
         if atom[pair[0]]['elem'] == 'H': # if we pick H- - H
            for count_p,p in enumerate(pair):
               id_atom.pop(atom[p]['mol_id'],None)
               id_atom[tot_mol].append(p)
               atom[p]['mol_id'] = tot_mol
               atom[p]['resi'] = 'HYD'
               atom[p]['atom_name'] = 'H{}'.format(count_p+1)
            for count_p,p in enumerate(other_p):
               id_atom.pop(atom[p]['mol_id'],None)
               id_atom[tot_mol+1].append(p)
               atom[p]['mol_id'] = tot_mol+1
               atom[p]['resi'] = 'IOD'
               atom[p]['atom_name'] = 'I{}'.format(count_p+1)
               if count_p == 1:
                  atom[p]['coord'] = atom[other_p[0]]['coord'] + np.array([2.823,0,0])
                  #atom[p]['coord'] = wrap_coord(atom[p]['coord'],box)
         if atom[pair[0]]['elem'] == 'I': # if we pick I - - I
            for count_p,p in enumerate(pair):
               id_atom.pop(atom[p]['mol_id'],None)
               id_atom[tot_mol].append(p)
               atom[p]['mol_id'] = tot_mol
               atom[p]['resi'] = 'IOD'
               atom[p]['atom_name'] = 'I{}'.format(count_p+1)
            for count_p,p in enumerate(other_p):
               id_atom.pop(atom[p]['mol_id'],None)
               id_atom[tot_mol+1].append(p)
               atom[p]['mol_id'] = tot_mol+1
               atom[p]['resi'] = 'HYD'
               atom[p]['atom_name'] = 'H{}'.format(count_p+1)
               if count_p == 1:
                  atom[p]['coord'] = atom[other_p[0]]['coord'] + np.array([0.742,0,0])
                  #atom[p]['coord'] = wrap_coord(atom[p]['coord'],box)
      id_resi = [ (atom[id_atom[mol_id][0]]['resi'],mol_id) for mol_id in id_atom]
      id_resi.sort(key = lambda x: x[0])  

      atom_count = 0
      for resi_count,(resi,mol_id) in enumerate(id_resi):
         for a in id_atom[mol_id]:
            line = formatting.format(resi_count+1,atom[a]['resi'],atom[a]['atom_name'],atom_count+1,atom[a]['coord'][0]/10,atom[a]['coord'][1]/10,atom[a]['coord'][2]/10)
            f.write(line)
            atom_count += 1
      f.write("{0:<8.3f} {0:<8.3f} {0:<8.3f}\n".format(box/10))
   
   return 0

def wrap_coord(coord,box):
   for count_i,i in enumerate(coord):
      if i > box: coord[count_i] -= box
      if i < 0 : coord[count_i] += box
   return coord

def get_nmol(id_atom,atom,resi):
   nmol = []
   for i in resi:
      num = len([ mol_id  for mol_id in id_atom if atom[id_atom[mol_id][0]]['resi'] == i])
      nmol.append(num)
   return nmol
   

def write_nvt(resi,Temp,nmol,box,rseed,cycle,polar_flag=False,switch=False):
   if switch:
      title = '@simnum_switch'
      run_time = 5000 # 5ps
   else:
      title = '@simnum'
      run_time = 1000 # 1ps
   with open('md-nvt.inp','w') as f:
      f.write("! this is necessart to work with non-interger total charge\n")
      f.write("bomlev -1\n\n")

      f.write("set simnum = {}\n".format(cycle))
      #f.write("calc prevsim = @simnum -1\n")
      f.write("set temp = {}\n\n".format(Temp))

      f.write("! read the psf and coordinate file\n")
      f.write("read rtf card name ff/lin.rtf\n")
      f.write("read para card name ff/lin.prm\n\n")

      for count_i,i in enumerate(resi):
         if nmol[count_i] == 0: continue
         f.write("set residue = {}\n".format(i))
         f.write("READ SEQUENCE @residue {}\n".format(nmol[count_i]))
         f.write("generate @residue first none last none setup warn\n") 
      f.write("OPEN READ CARD UNIT 10 NAME {}.crd\n".format(title))
      f.write("read coor card name {}.crd\n\n".format(title))


      if polar_flag:
         f.write("coor sdrude\n")
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
      f.write("crystal build noper 0 cutoff 12\n\n")

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

      f.write("set NSTEPS = {}\n".format(run_time))
      f.write("set printfreq = 1000\n")
      f.write("OPEN WRITE CARD UNIT 62 NAME {}.nose ! info for temperaature\n".format(title))
      f.write("open unit 35 write form name {}.rst ! restart file\n".format(title))
      f.write("open unit 32 write unfo name {}.trj ! trajectory\n".format(title))
      f.write("open unit 33 write form name {}.kunit ! Temperature\n".format(title))
      f.write("open unit 34 write unfo name {}.vel ! velocity\n\n".format(title))

      f.write("DYNA vv2 start -\n")
      f.write("nstep @NSTEPS time 0.001 -\n")
      f.write("iunread -1 -\n")
      f.write("iunwrite 35 iuncrd 32 iunvel 34 kunit 33 -\n")
      f.write("nsavc @printfreq nsavv @printfreq nprint @printfreq isvfrq @printfreq -\n")
      f.write("iprfrq @NSTEPS ntrfrq @NSTEPS - ! what's used for swm4 water\n")
      f.write("firstt @temp -\n")
      f.write("finalt @temp -\n")
      f.write("iasvel 1 iseed {:7d} -\n".format(rseed))
      f.write("iuno 62 nsnos @printfreq\n\n")

      f.write("open write unit 10 card name {}.crd\n".format(title))
      f.write("write coor card unit 10\n")
      f.write("close unit 10\n\n")

      f.write("stop\n")
   return
   

def run_nvt(resi,cycle,Temp,nmol,box,switch_flag,rseed=123,min_flag=False):
   if switch_flag:
      title = '{}_switch'.format(cycle)
   else:
      title = '{}'.format(cycle)
   write_nvtsh(cycle,switch_flag,min_flag=min_flag)
   subprocess.call("chmod 777 nvt.sh", shell=True)
   sublist = 'sbatch nvt.sh'.split()
   output = subprocess.Popen(sublist,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]
   output = str(output,'utf-8')
   jobids = [ m.split()[-1] for m in output.split("\n")[:-1]]
   monitor_jobs(jobids,user)
   
   #sublist = ('sh nvt.sh').split()
   #output = subprocess.Popen(sublist)
   #status = output.poll()
   #while status == None:
   #   status = output.poll()
   #   time.sleep(5)
   delete_wildcard('','cpt')
   return

def write_nvtsh(cycle,switch_flag,min_flag=False):
   if switch_flag:
      title = '{}_switch'.format(cycle)
   else:
      title = '{}'.format(cycle)
   with open('nvt.sh','w') as f:
      f.write("#!/bin/bash\n\n") 

      f.write("#SBATCH --job-name nvt_re\n")
      f.write("#SBATCH -N 1\n")
      f.write("#SBATCH -n 8\n")
      f.write("#SBATCH -t 04:00:00\n")
      f.write("#SBATCH -A bsavoie\n")
      f.write("#SBATCH -o nvt.out\n")
      f.write("#SBATCH -e nvt.err\n\n")

      f.write("module load gcc/8.3.0/\n")
      f.write("module load openmpi/3.1.4\n\n")
      f.write("GRO=/home/lin1209/install/gromacs_2/bin/gmx\n")

      f.write("title={}\n".format(title))
      if min_flag:
         write_nvtmdp(time=10)
         f.write("$GRO grompp -f em_out.mdp -c ${title}.gro -o em.tpr -maxwarn 1\n")
         f.write("$GRO mdrun -s em.tpr -deffnm em\n")
         f.write("rm mdout.mdp\n")
         f.write("$GRO grompp -f nvt.mdp -c em.gro -o ${title}.tpr -maxwarn 1\n") 
         f.write("$GRO mdrun -s ${title}.tpr -deffnm ${title} -nt 8\n")
         f.write("$GRO trjconv -f ${title}.trr -s ${title}.tpr -o trj.gro <<EOF\n")
         f.write("0\n")
         f.write("EOF\n")
         f.write("rm em.*\n")
         f.write("rm ${title}.xtc\n")
      else:
         write_nvtmdp(time=5)
         f.write("titlenext={}\n".format(cycle+1))
         f.write("$GRO grompp -f nvt.mdp -c ${title}.gro -o ${titlenext}.tpr -maxwarn 1\n") 
         f.write("$GRO mdrun -s ${titlenext}.tpr -deffnm ${titlenext} -nt 8\n")
         f.write("$GRO trjconv -f ${titlenext}.trr -s ${titlenext}.tpr -o trj.gro <<EOF\n")
         f.write("0\n")
         f.write("EOF\n")
         f.write("rm ${titlenext}.xtc\n")
      f.write("rm mdout.mdp\n")
   return

# time: simulation length (ps)
# print_freq: printing freq (ps)
def write_nvtmdp(time=1,print_freq=1):
   time = int(time*1000) # input
   print_freq = int(print_freq*1000) # saving trajectory frequency
   with open('nvt.mdp','w') as f:
      f.write("; RUN CONTROL PARAMETERS\n")
      f.write("integrator               = md\n")
      f.write("; Start time and timestep in ps\n")
      f.write("tinit                    = 0\n")
      f.write("dt                       = 0.001\n")
      f.write("nsteps                   = {}\n".format(time))
      f.write("; For exact run continuation or redoing part of a run\n")
      f.write("init_step                = 0\n")
      f.write("; mode for center of mass motion removal\n")
      f.write("comm-mode                = Linear\n")
      f.write("; number of steps for center of mass motion removal\n")
      f.write("nstcomm                  = 1\n")
      f.write("; group(s) for center of mass motion removal\n")
      f.write("comm-grps                = \n\n")

      f.write("; LANGEVIN DYNAMICS OPTIONS\n")
      f.write("; Friction coefficient (amu/ps) and random seed\n")
      f.write("bd-fric                  = 0\n")
      f.write("ld-seed                  = 1993\n\n")

      f.write("; ENERGY MINIMIZATION OPTIONS\n")
      f.write("; Force tolerance and initial step-size\n")
      f.write("emtol                    = 10\n")
      f.write("emstep                   = 0.01\n")
      f.write("; Max number of iterations in relax_shells\n")
      f.write("niter                    = 20\n")
      f.write("; Step size (ps^2) for minimization of flexible constraints\n")
      f.write("fcstep                   = 0\n")
      f.write("; Frequency of steepest descents steps when doing CG\n")
      f.write("nstcgsteep               = 1000\n")
      f.write("nbfgscorr                = 10\n\n")

      f.write("; OUTPUT CONTROL OPTIONS\n")
      f.write("; Output frequency for coords (x), velocities (v) and forces (f)\n")
      f.write("nstxout                  = {}\n".format(print_freq))
      f.write("nstvout                  = {}\n".format(print_freq))
      f.write("nstfout                  = {}\n".format(print_freq))
      f.write("; Output frequency for energies to log file and energy file\n")
      f.write("nstlog                   = 100\n")
      f.write("nstenergy                = 100\n")
      f.write("; Output frequency and precision for xtc file\n")
      f.write("nstxout-compressed       = 500\n")
      f.write("compressed-x-precision   = 1000\n\n")

      f.write("; NEIGHBORSEARCHING PARAMETERS\n")
      f.write("; cut-off scheme (Verlet: particle based cut-offs, group: using charge groups)\n")
      f.write(";cutoff-scheme            = group ; default in new gromacs is verlet, before 4.6 only support group, that's probably what spoel's group used\n")
      f.write("; nblist update frequency\n")
      f.write("nstlist                  = 5\n")
      f.write("; ns algorithm (simple or grid)\n")
      f.write("ns-type                  = grid\n")
      f.write("; Periodic boundary conditions: xyz (default), no (vacuum)\n")
      f.write("; or full (infinite systems only)\n")
      f.write("pbc                      = xyz\n")
      f.write("; nblist cut-off\n")        
      f.write("rlist                    = 1.2\n\n")

      f.write("; OPTIONS FOR ELECTROSTATICS AND VDW\n")
      f.write("; Method for doing electrostatics\n")
      f.write("coulombtype              = PME\n")
      f.write("rcoulomb-switch          = 0\n")
      f.write("rcoulomb                 = 1.2\n")
      f.write("; Relative dielectric constant for the medium and the reaction field\n")
      f.write("epsilon-r                = 1\n")
      f.write("epsilon-rf               = 1\n")
      f.write("; Method for doing Van der Waals\n")
      f.write("vdwtype                 = cutoff\n")
      f.write("; cut-off lengths       \n")
      f.write("rvdw-switch              = 0\n")
      f.write("rvdw                     = 1.2\n")
      f.write("; Apply long range dispersion corrections for Energy and Pressure\n")
      f.write("DispCorr                 = EnerPres\n")
      f.write("; Extension of the potential lookup tables beyond the cut-off\n")
      f.write("table-extension          = 1\n")
      f.write("; Seperate tables between energy group pairs\n")
      f.write("energygrp-table          = \n")
      f.write("; Spacing for the PME/PPPM FFT grid\n")
      f.write("fourierspacing           = 0.12\n")
      f.write("; FFT grid size, when a value is 0 fourierspacing will be used\n")
      f.write("fourier-nx               = 0\n")
      f.write("fourier-ny               = 0\n")
      f.write("fourier-nz               = 0\n")
      f.write("; EWALD/PME/PPPM parameters\n")
      f.write("pme-order                = 4\n")
      f.write("ewald-rtol               = 1e-05\n")
      f.write("ewald-geometry           = 3d\n")
      f.write("epsilon-surface          = 0\n\n\n")


      f.write("; OPTIONS FOR WEAK COUPLING ALGORITHMS\n")
      f.write("; Temperature coupling  \n")
      f.write("; Groups to couple separately\n")
      f.write("tcoupl                   = nose-hoover\n")
      f.write("tc-grps                  = system\n")
      f.write("tau-t                    = 1\n")
      f.write("ref-t                    = 300\n")
      f.write("; Time constant (ps) and reference temperature (K)\n")
      f.write("; Pressure coupling     \n")
      f.write("; Time constant (ps), compressibility (1/bar) and reference P (bar)\n")
      f.write("pcoupl                   = no \n")
      f.write("pcoupltype               = isotropic\n")
      f.write("tau-p                    = 5\n")
      f.write("compressibility          = 5e-5\n")
      f.write("ref-p                    = 1\n\n\n")


      f.write("; GENERATE VELOCITIES FOR STARTUP RUN\n")
      f.write("gen-vel                  = yes\n")
      f.write("gen-temp                 = 300\n")
      f.write("gen-seed                 = -1\n\n")

      f.write("; OPTIONS FOR BONDS    \n")
      f.write("constraints              = all-bonds\n")
      f.write("; Type of constraint algorithm\n")
      f.write("constraint-algorithm     = Lincs\n")
      f.write("; Highest order in the expansion of the constraint coupling matrix\n")
      f.write("lincs-order              = 8       ; default 4\n")
      f.write("; Number of iterations in the final step of LINCS. 1 is fine for\n")
      f.write("; normal simulations, but use 2 to conserve energy in NVE runs.\n")
      f.write("; For energy minimization with constraints it should be 4 to 8.\n")
      f.write("lincs-iter               = 1\n")
      f.write("; Lincs will write a warning to the stderr if in one step a bond\n")
      f.write("; rotates over more degrees than\n")
      f.write("lincs-warnangle          = 30\n")

   return

def find_cand_old(atom):
    box = 250
    cand = [[],[]]
    # search NR1: H-I < 4.7 A
    for a in atom:
      if atom[a]['elem'] == 'H' and atom[a]['resi'] == 'HYD':
         for b in atom:
            if atom[b]['elem'] == 'I' and  atom[b]['resi'] == 'IOD' and atom[b]['mol_id'] != atom[a]['mol_id']:
               d = unwrap_d(atom[a]['coord']-atom[b]['coord'],box) 
               #if np.linalg.norm(d) <= 4.7:
               if np.linalg.norm(d) <= 8:
                  cand[0].append((a,b)) 

    # search NR2: H-H <1.610 I-I<4.176
    for a in atom:
      if atom[a]['resi'] == 'HYI' and atom[a]['elem'] == 'H':
         for b in atom:
            if atom[b]['resi'] == 'HYI' and atom[b]['elem'] == 'H' and atom[b]['mol_id'] != atom[a]['mol_id']:
               d = unwrap_d(atom[a]['coord']-atom[b]['coord'],box)
               #if np.linalg.norm(d) <= 1.610:
               if np.linalg.norm(d) <= 8:
                  if (a,b) in cand[1] or (b,a) in cand[1]: continue
                  cand[1].append((a,b)) 
    for a in atom:
      if atom[a]['resi'] == 'HYI' and atom[a]['elem'] == 'I':
         for b in atom:
            if atom[b]['resi'] == 'HYI' and atom[b]['elem'] == 'I' and atom[b]['mol_id'] != atom[a]['mol_id']:
               d = unwrap_d(atom[a]['coord']-atom[b]['coord'],box)
               #if np.linalg.norm(d) <= 4.176:
               if np.linalg.norm(d) <= 8:
                  if (a,b) in cand[1] or (b,a) in cand[1]: continue
                  cand[1].append((a,b)) 
               
    return cand



# find candidate pair
def find_cand():
    box = 250
    cand = [[],[]]
    cand_tmp = []
    f = open('trj.gro','r')
    frame_count = -1
    for lc,line in enumerate(f):
       fields = line.split()
       # initialization at every frame
       if len(fields) > 4 and fields[5] == 't=':
         start_lc = lc
         frame_count += 1  
         atom = {}
         continue
       if lc - start_lc == 1:
         tot_atom = int(fields[0])
         continue
       if frame_count == 0: continue
       if (lc-start_lc) == (tot_atom+2):
          # search NR1: H-I < 4.7 A
          for a in atom:
            if atom[a]['elem'] == 'H' and atom[a]['resi'] == 'HYD':
               for b in atom:
                  if atom[b]['elem'] == 'I' and  atom[b]['resi'] == 'IOD' and atom[b]['mol_id'] != atom[a]['mol_id']:
                     d = unwrap_d(atom[a]['coord']-atom[b]['coord'],box) 
                     #if np.linalg.norm(d) <= 4.7:
                     if np.linalg.norm(d) <= 8:
                        if (a,b) in cand[0] or (b,a) in cand[0]: continue
                        cand[0].append((a,b)) 

          # search NR2: H-H <1.610 I-I<4.176
          for a in atom:
            if atom[a]['resi'] == 'HYI' and atom[a]['elem'] == 'H':
               for b in atom:
                  if atom[b]['resi'] == 'HYI' and atom[b]['elem'] == 'H' and atom[b]['mol_id'] != atom[a]['mol_id']:
                     d = unwrap_d(atom[a]['coord']-atom[b]['coord'],box)
                     #if np.linalg.norm(d) <= 1.610:
                     if np.linalg.norm(d) <= 8:
                        if (a,b) in cand[1] or (b,a) in cand[1]: continue
                        cand[1].append((a,b)) 
          for a in atom:
            if atom[a]['resi'] == 'HYI' and atom[a]['elem'] == 'I':
               for b in atom:
                  if atom[b]['resi'] == 'HYI' and atom[b]['elem'] == 'I' and atom[b]['mol_id'] != atom[a]['mol_id']:
                     d = unwrap_d(atom[a]['coord']-atom[b]['coord'],box)
                     #if np.linalg.norm(d) <= 4.176:
                     if np.linalg.norm(d) <= 8:
                        if (a,b) in cand[1] or (b,a) in cand[1]: continue
                        cand[1].append((a,b)) 
          continue
       info = read_gro_format(line)
       atom_count = info['atom_count']-1
       atom[atom_count] = {}
       atom[atom_count]['mol_id'] = info['mol_id']-1
       atom[atom_count]['resi'] = info['resname']
       atom[atom_count]['coord'] = np.array([info['x']*10,info['y']*10,info['z']*10])
       atom[atom_count]['elem'] = info['element']
       atom[atom_count]['atom_name'] = info['atom_name']
               
    return cand

def crd_formatting():
    formatting = '{:>5d}{:>5d} {:<4s} {:<4s}{:>10.5f}{:>10.5f}{:>10.5f} {:<4s} {:<4s}{:10.5f}'
    return

def gro_formatting():
    formatting = "{:5d}{:<5s}{:>5s}{:>5d}{:> 8.3f}{:> 8.3f}{:> 8.3f}\n"


# Description: Performed the periodic boundary unwrap of the geometry
def unwrap_geo(geo,box,total_atoms,total_mol,molids):

    index = [ [] for i in range(0,total_mol)] # do not use [[]]*len, they will have the same reference
    for count_i,i in enumerate(molids):
         index[i-1] += [count_i]

    b_2 = [ box/2,box/2,box/2]

    for i in index:
      for count_j,j in enumerate(i):
         fix = 0
         if count_j == 0: continue
         d = geo[j] - geo[i[0]]
         for count_k,k in enumerate(d):
            if k < -b_2[count_k]: 
               d[count_k] += 2*b_2[count_k]
            if k > b_2[count_k]: 
               d[count_k] -= 2*b_2[count_k]
         geo[j] = geo[i[0]] + d

    return geo

# Description: Performed the periodic boundary unwrap of the geometry
def unwrap_d(d,box):

    b_2 = [ box/2,box/2,box/2]
    for count_i,i in enumerate(d):
      if i < -b_2[count_i]: 
         d[count_i] += 2*b_2[count_i]
      if i > b_2[count_i]: 
         d[count_i] -= 2*b_2[count_i]

    return d
    
def split_string(word): 
    return list(word) 

def read_crd(filename):
    atom = {}
    id_atom = {} # key:mol_id value:atoms
    with open(filename,'r') as f:
       for lc,lines in enumerate(f):
         if lc < 4 : continue
         fields = lines.split()
         atom_count = int(fields[0])-1
         atom[atom_count] = {}
         atom[atom_count]['mol_id'] = int(fields[1])-1
         atom[atom_count]['coord'] = np.array(fix_coord(lines))
         atom[atom_count]['resi'] = fields[2]
         atom[atom_count]['elem'] = split_string(fields[3])[0] 
         atom[atom_count]['atom_name'] = fields[3]
         atom[atom_count]['resi_id'] = fields[-2]
         if atom[atom_count]['mol_id'] not in list(id_atom.keys()): 
            id_atom[atom[atom_count]['mol_id']] = [] 
         id_atom[atom[atom_count]['mol_id']].append(atom_count)
    return atom,id_atom

def read_gro(filename):
    atom = {}
    id_atom = {} # key:mol_id value:atoms
    with open(filename,'r') as f:
       for lc,lines in enumerate(f):
         fields = lines.split()
         if lc == 0 : continue
         if lc == 1: 
            tot_atom = int(fields[0])
            continue
         if lc == tot_atom+2:
            box = float(fields[0])
            continue
         info = read_gro_format(lines)
         atom_count = info['atom_count']-1
         atom[atom_count] = {}
         atom[atom_count]['mol_id'] = info['mol_id']-1
         atom[atom_count]['resi'] = info['resname']
         atom[atom_count]['coord'] = np.array([info['x']*10,info['y']*10,info['z']*10])
         atom[atom_count]['elem'] = info['element']
         atom[atom_count]['atom_name'] = info['atom_name']
         if atom[atom_count]['mol_id'] not in list(id_atom.keys()): 
            id_atom[atom[atom_count]['mol_id']] = [] 
         id_atom[atom[atom_count]['mol_id']].append(atom_count)
    return atom,id_atom

def read_gro_format(line):
   info = {}
   info['mol_id'] = int(line[0:5]) 
   info['resname'] = str(''.join([ i for i in line[5:10] if i != ' ']))
   info['atom_name'] = str(''.join([ i for i in line[10:15] if i != ' ']))
   info['element'] = split_string(info['atom_name'])[0]
   info['atom_count'] = int(line[15:20])
   info['x'] = float(line[20:28])
   info['y'] = float(line[28:36])
   info['z'] = float(line[36:44])

   return info

def resi_gro(num_resi):
   num = []
   resi = []
   for i in split_string(num_resi):
     if is_number(i): num += i 
     else: resi += i
   return int(''.join(num))-1,''.join(resi)
   

# deal with crd file format
# sometimes the coord is too long that there are no space between
def fix_coord(lines):
    fields = lines.split()
    coord = ' '.join(fields[4:-3])
    coord = coord.split('.')
    new_coord = []
    for i in coord: 
      new_coord += i.split()
    final_coord = []
    for i in new_coord:
      if is_number(i) is False: 
         tmp = i.split('-')
         tmp[1] = '-'+tmp[1]
         final_coord += tmp 
      else:
         final_coord += [i]
    final_coord = [float('.'.join(final_coord[i*2:i*2+2])) for i in range(0,3)]
    return final_coord
         

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
         

# Normal Drude fitting, fit all heavy atoms
def fit_Drude_custom(driver_dir,folder,args,tag):

    global run_dir

    run_dir = '{}/{}'.format(driver_dir,folder)
    os.chdir(run_dir) 
 
    nonpolar = '/depot/bsavoie/data/200819_taffi_benchmarks_batch1/' 
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

    print("Starting QM polar calculation...")
    for q in [-0.5,0.5]:
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

    jobids = []
    print("Starting CHELPG potential calculation...")
    for retry in range(0,1):
       for q in [-0.5,0.5]:
          for count_i,i in enumerate(list(range(0,nconf))):
            fitfolder = 'fitcharge_{}'.format(i)
            xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
            element,geo = xyz_parse(xyz)
            total_pc,pc_layer1,pc_layer2,pc_layer3 = get_pc(fitfolder)
            jobids += calc_CHELPG(fitfolder,xyz,total_pc,pc_layer1,pc_layer2,pc_layer3,geo,element,q,'_{}'.format(q))
       monitor_jobs(jobids,user)
       clean_chelpg(nconf)

    jobids = []
    print("Starting reevaluation of potential at Connolly grid...")
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
    quit()
    
    
    """    
    jobid = []
    for i in range(0,100):
      fitfolder = 'fitcharge_{}'.format(i)
      if os.path.isdir(run_dir+'/'+fitfolder) is False:
         os.makedirs(run_dir+'/'+fitfolder)
      xyz = '{}/charges/CHELPG_calcs/configs/{}/{}.xyz'.format(run_dir,i,i)
      total_pc,pc_layer1,pc_layer2,pc_layer3 = place_perturbation(run_dir,fitfolder,xyz)
      element,geo = xyz_parse(xyz)
      tmpjob = calc_CHELPG(run_dir,fitfolder,xyz,total_pc,pc_layer1,pc_layer2,pc_layer3,c,geo,element,user,'_neg')
      jobid += tmpjob
    monitor_jobs(jobid,user)
    quit()
    """ 

    """
    jobid = []
    for i in range(0,100):
      fitfolder = 'fitcharge_{}'.format(i)
      xyz = '{}/charges/CHELPG_calcs/configs/{}/{}.xyz'.format(run_dir,i,i)
      total_pc,pc_layer1,pc_layer2,pc_layer3 = place_perturbation(run_dir,fitfolder,xyz)
      connolly_grid = get_connolly_grid(run_dir,fitfolder,xyz)
      tmpjob = calc_connolly_pot(run_dir,fitfolder,total_pc,connolly_grid,c,user,'_neg')
      jobid += tmpjob
    monitor_jobs(jobid,user)
    quit()
    """
    """ 
    for i in range(0,100):
      fitfolder = 'fitcharge_{}'.format(i)
      xyz = '{}/charges/CHELPG_calcs/configs/{}/{}.xyz'.format(run_dir,i,i)
      total_pc,pc_layer1,pc_layer2,pc_layer3 = place_perturbation(run_dir,fitfolder,xyz)
      connolly_grid = get_connolly_grid(run_dir,fitfolder,xyz)
      rewrite_pot_wrap(run_dir,fitfolder,total_pc,connolly_grid,c,user,'_neg')
    """
    
    quit()

   
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

def place_perturbation(calc_folder,xyz):

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
