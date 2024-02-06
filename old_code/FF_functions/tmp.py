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
import place_charge,generate_connolly
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

    parser = argparse.ArgumentParser(description='takes in the dimer xyz file in configs folder, do drude relaxation, this program assumes there is a FF folder that contains lin.prm and lin.rtf in vdw folder')

    #optional arguments                                                                                                                                                        
    parser.add_argument('-xyz', dest='xyz', default='fitcharge.xyz',
                        help = 'xyz file name, default: fitcharge.xyz')
    parser.add_argument('--silent',dest='silent', default=False, const=True, action='store_const',
                        help = 'Shortcircuits all script related print statements') 
    parser.add_argument('-tag',dest='tag',default='',
                        help = 'tag for FF folder and relax drude, this assumes FF name is FF_tag and output drude: _drude_tag.xyz, default:None')

    args=parser.parse_args(argv)    

    title = args.xyz.split('.')[0]

    if not args.silent: 
      print("PROGRAM CALL: relax_drude.py {}".format(' '.join([ i for i in argv]))) 

    if args.tag != '':
      FFpath = '../../FF_{}'.format(args.tag)
      drude_dir = 'relax_drude_{}'.format(args.tag)
      xyzname = '{}_drude_{}.xyz'.format(title,args.tag)
    else:
      FFpath = '../../FF'
      drude_dir = 'relax_drude'
      xyzname = '{}_drude.xyz'.format(title)

    if not os.path.isdir(FFpath):
         print("FF folder not found. Exiting....")
         quit()

    if not os.path.isdir(drude_dir):  os.makedirs(drude_dir)
    shutil.copy2('{}/lin.prm'.format(FFpath),'{}/lin.prm'.format(drude_dir))
    shutil.copy2('{}/lin.rtf'.format(FFpath),'{}/lin.rtf'.format(drude_dir))


    os.chdir(drude_dir)

    # read in nonpolar config dimer xyz file
    elem_a,elem_b,geo_a,geo_b,types_a,types_b = scrape_xyz('../'+args.xyz,parse_charges=0)
    
    # write CHARMM relaxation inp
    Write_inp(elem_a,elem_b,geo_a,geo_b)

    quit()

    # Run charmm drude relaxation
    Write_exe()
    command = 'sh relax.sh'
    output = subprocess.Popen(command.split(),stdin=PIPE, stdout=PIPE, stderr=PIPE).communicate()

    # Read relax geo and charge and write it in xyz format
    charges = read_charge()
    Write_drude_xyz(charges,'../'+xyzname,elem_a,elem_b,types_a,types_b)
   
    os.chdir('..')


    if not os.path.isfile(xyzname):
         print("{} relaxation failed".format(args.xyz))
         quit()
   
    # if successful, rm relax folder
    else:
         shutil.rmtree(drude_dir)
    
    return
    

def Write_drude_xyz(charges,filename,elem_a,elem_b,types_a,types_b):
    with open('drude.crd','r') as f:
        count = 0
        flag = 0
        for lc,lines in enumerate(f):
            fields = lines.split()
            if lc < 3: continue
            if lc == 3:
               natoms = int(fields[0])
               ele = []
               mol = []
               geos = np.zeros([natoms,3])
               flag = 1
               continue
            if flag == 1:
               geos[count] = np.array([ float(i) for i in fields[4:7] ])
               tmpele = []
               tmp = split_str(fields[3])
               for i in tmp:
                  if i.isnumeric(): break
                  else: tmpele += i
               ele.append("".join(tmpele))
               mol.append(int(fields[1])-1)
               count += 1
               if count == natoms: break
    with open(filename,'w') as f:
         atom_count = 0
         f.write('{} \n\n'.format(natoms))
         for count_i,i in enumerate(ele):
            f.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} {:< 5d} {:< 12.6f} \n'.format(i,geos[count_i,0],geos[count_i,1],geos[count_i,2],mol[count_i],charges[count_i]))
         """ 
         for count_i,i in enumerate(elem_a):
            f.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} {:< 5d} {:< 12.6f} {}\n'.format(i,geos[atom_count,0],geos[atom_count,1],geos[atom_count,2],0,charges[atom_count],types_a[count_i]))
            atom_count += 1
            if i != 'H':
               f.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} {:< 5d} {:< 12.6f} {}\n'.format(i,geos[atom_count,0],geos[atom_count,1],geos[atom_count,2],0,charges[atom_count],types_a[count_i]))
               atom_count += 1
         for count_i,i in enumerate(elem_b):
            f.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} {:< 5d} {:< 12.6f} {}\n'.format(i,geos[atom_count,0],geos[atom_count,1],geos[atom_count,2],1,charges[atom_count],types_b[count_i]))
            atom_count += 1
            if i != 'H':
               f.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} {:< 5d} {:< 12.6f} {}\n'.format(i,geos[atom_count,0],geos[atom_count,1],geos[atom_count,2],1,charges[atom_count],types_b[count_i]))
               atom_count += 1
         """
    # This error could block lone pairs
    #if atom_count != natoms:
    #     print("ERROR: write out number of atoms doesn't match number of atoms in crd file. Exiting....")
    #     quit()
    return
            
def split_str(word): 
    return [char for char in word]



def read_charge():
    with open('relax.log','r') as f:
        flag = 0
        charges = []
        for lc,lines in enumerate(f):
            fields = lines.split()
            if len(fields) < 4: continue
            if fields[1] == 'scalar' and fields[2] == 'charge' and fields[3] == 'show':
               flag = 1
               continue
            if flag == 1 and fields[1] != 'LIN':
               flag = 0
               break 
            if flag == 1:
               charges.append(float(fields[6]))
               continue
    return charges

def Write_exe():
    with open('relax.sh','w') as f:
         f.write('#!/bin/bash\n')
         f.write('export PATH=\"/apps/cent7/xalt/bin:/apps/brown/openmpi.slurm/slurm-19.05/3.1.4_intel-17.0.1.132/bin:/home/lin1209/anaconda3/bin:/home/lin1209/anaconda3/condabin:/apps/cent7/vmd/vmd-1.9.3/bin:/usr/lib64/qt-3.3/bin:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/bin/intel64:/apps/cent7/intel/inspector_2017.1.1.484836/bin64:/apps/cent7/intel/advisor_2017.1.1.486553/bin64:/apps/cent7/intel/vtune_amplifier_xe_2017/bin64:/apps/cent7/intel/opencl-1.2-sdk-6.3.0.1904/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/hpss/bin:/opt/hsi/bin:/opt/ibutils/bin:/opt/thinlinc/bin:/bin:/opt/hpss/bin:/opt/hsi/bin:/opt/intel/mic/bin:/sbin:/usr/bin:/usr/lib64/alliance/bin:/usr/lib64/ccache:/usr/lib64/qt-3.3/bin:/usr/libexec/lightdm:/usr/libexec/sdcc:/usr/local/bin:/usr/local/sbin:/usr/lpp/mmfs/bin:/usr/sbin:/usr/site/rcac/scripts:/usr/site/rcac/sbin:/usr/site/rcac/bin:/usr/site/rcac/scripts:/usr/site/rcac/support_scripts:/home/lin1209/bin:/depot/bsavoie/apps/openbabel/bin:/home/lin1209/install/gnome/bin:/home/lin1209/install/charmm/exec/gnu\"\n\n')
         f.write('export LD_LIBRARY_PATH=\"/apps/cent7/vmd/vmd-1.9.3/lib:/apps/cent7/intel/impi/2017.1.132/lib64:/apps/cent7/intel/itac/2017.1.024/intel64/slib:/apps/cent7/intel/itac/2017.1.024/lib:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/ipp/lib/intel64:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/mkl/lib/intel64:/apps/cent7/intel/compilers_and_libraries_2017.1.132/linux/tbb/lib/intel64/gcc4.7:/apps/cent7/intel/opencl-1.2-sdk-6.3.0.1904/lib64:/apps/cent7/intel/debugger_2017/libipt/intel64/lib:/apps/cent7/intel/debugger_2017/iga/lib\"\n\n')
         f.write('charmm < relax.inp > relax.log\n')
         f.write('wait')
    return


def Write_inp(elem_a,elem_b,geo_a,geo_b):
      
    with open('relax.inp','w') as f:
         f.write('! this is necessart to work with non-interger total charge\n'+\
                 'bomlev -1\n\n'+\
                 'read rtf card name lin.rtf\n'+\
                 'read para card name lin.prm\n'+\
                 'read sequence lin 2\n\n'+
      
                 'generate lin first none last none setup warn drude dmass 0.4\n\n'+\

                 'read coor card free\n'+\
                 '* lin.qpos.1.crd\n'+\
                 '*\n')
         f.write('  {}\n'.format(len(elem_a)*2))
         for count_i,i in enumerate(geo_a):
           atomname = elem_a[count_i]+str(count_i+1)
           f.write("{:< 4d} {:< 4d} LIN {:} {:< 10.6f} {:< 10.6f} {:< 10.6f} LIN {:d} {:7.5f}\n"\
            .format(count_i+1,1,atomname,i[0],i[1],i[2],1,0.00000))
         for count_i,i in enumerate(geo_b):
           atomname = elem_a[count_i]+str(count_i+1)
           f.write("{:< 4d} {:< 4d} LIN {:} {:< 10.6f} {:< 10.6f} {:< 10.6f} LIN {:d} {:7.5f}\n"\
            .format(count_i+1+len(elem_a),2,atomname,i[0],i[1],i[2],1,0.00000))
         f.write('\n\ncoor sdrude\n'+\
                 'coor shake\n\n\n'+\

                 'update inbfrq -1 ihbfrq 0 -'+\
                 '  switch atom vatom vswitch cutnb 990. ctofnb 950. ctonnb 900.\n'+\
                 'cons fix sele .not. type D* end\n'+\
                 'mini ABNR nstep 1000 nprint 20\n'+\
                 'coor print\n\n'+\

                 'scalar charge show\n'+\
                 'open write unit 10 card name drude.crd\n'+\
                 'write coor card unit 10\n')
    return



# Return the coordinates and atom types from an xyz file
# atom types are expected as a fifth column and molecule
# index is expected in the sixth column. These are atypical
# an xyz file, but automatically outputed by the vdw_gen.py 
# program.
def scrape_xyz(name,parse_charges=0):

    atom_count = 0
    with open(name,'r') as f:
        for count_j,j in enumerate(f):
            fields = j.split()

            # Grab the number of atoms in the geometry and initialize the parsing lists
            if count_j == 0:
                elem = ["X"]*int(fields[0])
                geo = np.zeros([int(fields[0]),3]) 
                types = ["X"]*int(fields[0])
                a_ind = []
                b_ind = []
                
                # Only parse charges if flag is enabled
                if parse_charges == 1:
                    charges =np.zeros([int(fields[0])])

            # Within the geometry block parse the geometry, charges, mol_ids, and atom types from the enhanced xyz file
            if count_j > 1 and len(fields) >= 7:
                geo[atom_count,:] = np.array([float(fields[1]),float(fields[2]),float(fields[3])])
                if fields[4] == '0':
                    a_ind += [atom_count]
                else:
                    b_ind += [atom_count]
                elem[atom_count] = fields[0]

                # Only parse charges if flag is enabled
                if parse_charges == 1:
                    charges[atom_count] = float(fields[5])
                types[atom_count] = fields[6]
                atom_count += 1    

    # If flag is enabled return the charges
    if parse_charges == 1:
        return [ i for count_i,i in enumerate(elem) if count_i in a_ind ],\
               [ i for count_i,i in enumerate(elem) if count_i in b_ind ],\
               geo[a_ind],geo[b_ind],\
               [ i for count_i,i in enumerate(types) if count_i in a_ind ],\
               [ i for count_i,i in enumerate(types) if count_i in b_ind ],\
               charges[a_ind],charges[b_ind]

    # Don't return the charges list if the flag isn't enabled
    if parse_charges == 0:
        return [ i for count_i,i in enumerate(elem) if count_i in a_ind ],\
               [ i for count_i,i in enumerate(elem) if count_i in b_ind ],\
               geo[a_ind],geo[b_ind],\
               [ i for count_i,i in enumerate(types) if count_i in a_ind ],\
               [ i for count_i,i in enumerate(types) if count_i in b_ind ]


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
