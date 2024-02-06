#!/bin/env python                                                                                                                                                             
# Author: Brett Savoie (brettsavoie@gmail.com)

import sys,argparse,os,datetime,fnmatch,os,re

# Add TAFFY Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from transify import *
from adjacency import *
from file_parsers import *
from id_types import *
from kekule import *

# NOTE: the order of the imports matters until I resolve the namespace issues (XXX DON'T USE import *)
import random
from numpy import *
import time
import ast
from scipy import *
from scipy.spatial.distance import *
from numpy.linalg import *
from shutil import move,copyfile

def main(argv):

    parser = argparse.ArgumentParser(description='Reads in a geometry and FF.db file (usually outputted from an intramolecular mode parse) '+\
                                                 'and writes CHARMM inp file for fit charge with force field writtern in CHARMM compatible format.')


    #required (positional) arguments                                                                                                  
    parser.add_argument('coord_files', help = 'A quoted list of the *.xyz files to be included in the sim, or a folder holding a previous md cycles. If a folder is supplied then the FF parameters and last configuration '+\
                                              'are used to start the new job along with the most recently fit VDW parameters.')

    #optional arguments    
    parser.add_argument('-FF',dest='FF_files', default="",
                         help = 'A quoted list of the force-field files. Formatting of the force-field files is as produced by the FF_gen.py/FF_extract.py programs)')

    parser.add_argument('-gens', dest='gens', default=2,
                        help = 'When atomtypes are being automatically determined (i.e. when the supplied *.xyz files do not have atomtype information) this variable controls the bond depth used to define unique atomtypes. (default: 2)')

    parser.add_argument('-npert', dest='npert', default=2,
                        help = 'number of perturbationn charges (default: 2)')

    parser.add_argument('-mag_c', dest='mag_c', default=0.5,
                        help = 'magnitude of perturbation charges (default: 0.5)')

    parser.add_argument('-ascale',dest='ascale',default=1,
                        help = 'scale for atomic polarizability (default: 1)')

    parser.add_argument('-o', dest='outputname', default='',
                        help = 'Sets the output filename prefix (default: VDW_selfterms)')

    parser.add_argument('-inpname', dest='inpname', default='fitcharge',
                        help = 'name for inp file (w/o .inp) (default: fitcharge)')

    parser.add_argument('-mixing_rule', dest='mixing_rule', default="wh",
                        help = 'Defines the mixing rule to be used for missing LJ parameters. Waldman-Hagler (wh) and Lorentz-Berthelot (lb) and "none", are valid options. When set to "none", the program will print a message and '+\
                               'exit if there are missing parameters. (default: "none")')

    parser.add_argument('-q', dest='q_list', default="none",
                        help = 'Controls the total charge on the molecule. By default, the rounded integer charge is used for each molecule (round). The program expects a list of integers (or none, or round for the default). '+\
                               'If less integers are supplied than the number of coord files, then the list is automatically expanded using the last supplied (or default) element. (default: None)') 

    parser.add_argument('-eps_scale', dest='eps_scale', default=1.0,
                        help = 'Sets scaling for the default UFF eps parameters. To increase configurational sampling, it can be advantageous '+\
                               'to reduce the stiffness of the VDW potentials. (default: 1.0; full value)')

    parser.add_argument('-sigma_scale', dest='sigma_scale', default=1.0,
                        help = 'Sets scaling for the default UFF sigma parameters. To increase configurational sampling, it can be advantageous '+\
                               'to reduce the LJ minima of the VDW potentials. (default: 1.0; full value)')

    parser.add_argument('-14_scale', dest='14_scale', default=0.0,
                        help = 'scale for 1-4 interaction, default: 0.0(for taffi)')

    parser.add_argument('--UFF', dest='force_UFF', default=0, action='store_const', const=1,
                        help = 'Forces the use of UFF parameters (regardless of the presence of corresponding parameters in the read FF files). (default: off)')

    parser.add_argument('--impropers', dest='improper_flag', default=False, action='store_const', const=True,
                        help = 'When present, the simulation will use improper dihedral terms. (default: False)')

    parser.add_argument('--fit', dest='fit_flag', default=False, action='store_const', const=True,
                        help = 'When present, the simulation will write fitcharge.inp. (default: False)')

    parser.add_argument('--polar', dest='polar_flag', default=False, action='store_const', const=True,
                        help = 'When present, will write polar CHARMM FF files')

    parser.add_argument('--nonC', dest='nonC_flag', default=False, action='store_const', const=True,
                        help = 'When present, only put Drude on non carbon heavy atoms (default: off)') 

    parser.add_argument('--lp', dest='lp_flag', default=False, action='store_const', const=True,
                        help = 'When present, will write LP CHARMM FF fiels')

    parser.add_argument('--overwrite', dest='overwrite', default=False, action='store_const', const=True,
                        help = 'When present, will skip ERROR message when there\'s atomic polarizability missing, this is useful when trying putting Drude particles on partial atoms')

    parser.add_argument('--remove_multi', dest='remove_multi', default=False, action='store_const', const=True,
                        help = 'When present, the simulation will remove bimodal angles from the simulation (default: False)')

    parser.add_argument('--force_read', dest='force_read_opt', default=False, action='store_const', const=True,
                        help = 'When present, the simulation will use the force-field parameters that are discovered, even if the type does not match the expectations of the lewis structure parse. '+\
                               'In particular, sometimes the user wants to use a harmonic dihedral type for an ostensibly flexible dihedral. If only the harmonic type is supplied, the script will '+\
                               'print a missing mode warning and exit. When this flag is supplied the script will proceed to use the dihedral that is discovered in the force-field file. (default: off)')

    # Make parse inputs
    args=parser.parse_args(argv)

    # Make fake input config
    option = dict(args.__dict__.items())
    
    other_opt = {}

    fun(option,other_opt,argv)
    



def fun(config,other_opt,argv):


    option = miss2default(config,other_opt)
    
    coord_files = option['coord_files']
    FF_files = option['FF_files']
    mixing_rule = option['mixing_rule']
    polar_flag = option['polar_flag']
    lp_flag = option['lp_flag']
    q_list = option['q_list']
    gens = option['gens']
    improper_flag = option['improper_flag']
    remove_multi = option['remove_multi']
    force_read_opt = option['force_read_opt']
    outputname = option['outputname']
    onefour_scale = float(option['14_scale'])


    # Converting each list input argument. 
    option['ascale'] = float(option['ascale'])
    option['mag_c'] = float(option['mag_c'])
    gens = int(gens)
    coord_files = coord_files.split()
    FF_files = FF_files.split()
    q_list = [ str(i) for i in q_list.split() ]

    mixing_rule = str(mixing_rule).lower()
    if mixing_rule not in ["wh","lb","none"]:
        print("ERROR: '{}' is not a valid argument for -mixing_rule. Valid options are 'wh' 'lb' and 'none'. Exiting...".format(mixing_rule))
        quit()

    # Parse output filename 
    if outputname != '':

        # Handle the special case where the script is run from without the output directory
        if outputname == ".":
            print("ERROR: the output folder must differ from the current folder.")
            quit()

        # Else, directly assign the output directory based on the -o argument
        else:
            Filename = outputname        

    # If the output folder is not specified then the program defaults to the root of the first supplied *xyz file
    else:
        Filename = coord_files[0].split('.')[0]
    
    # Check that the input is an .xyz file. 
    # Check that the input is an .xyz file. 
    for i in coord_files:
         if i.split('.')[-1] != 'xyz':
             print("ERROR: Check to ensure that the input coordinate file(s) are in .xyz format.")
             quit()
         elif os.path.isfile(i) != True:
             print("ERROR: Specified *.xyz file ({}) does not exit. Please check the path. Exiting...".format(i))
             quit()

    # Check that the supplied FF files exist
    for i in FF_files:
        if os.path.isfile(i) != True:
            print("ERROR: Specified FF_file ({}) does not exist. Please check the path. Exiting...".format(i))
            quit()

    # If the output directory already exists, check for md cycle information to avoid overwriting data    
    if os.path.isdir(Filename) is False:

        os.makedirs(Filename)
    sys.stdout = Logger(Filename)
    print("PROGRAM CALL: python {}\n".format(' '.join([ i for i in sys.argv])))

    # Catenate the FF files and copy them to the output folder
    FF_all = Filename+"/"+Filename.split('/')[-1]+'.db'
    with open(FF_all,'w') as f:
        for i in FF_files:
            with open(i,'r') as ff_file:
                for j in ff_file:
                    f.write(j)
            f.write("\n")

    # Grab data on each molecule being added to the MD simulation.
    Data = get_data(FF_all,coord_files,q_list,gens,Improper_flag=improper_flag,force_read_opt=force_read_opt,remove_multi=remove_multi)
    alpha_dict = {}
    for i in Data:
         print(Data[i]['Charges'])
         for count_j,j in enumerate(Data[i]['Atom_types']):
            alpha_dict[j] = Data[i]['alphas'][count_j]



    # Write the lammps datafile (the function returns a dictionary, Atom_type_dict, holding the mapping between atom_types and the lammps id numbers; this mapping is needed for setting fixes)
    print("Writing CHARMM inpfile ({})...".format(Filename+'.inp'))

    # Generate VDW parameters    
    VDW_params = initialize_VDW(sorted(set(Data[coord_files[0]]['Atom_types'])),sigma_scale=option['sigma_scale'],eps_scale=option['eps_scale'],VDW_FF=Data[coord_files[0]]["VDW_params"],\
                                Force_UFF=option['force_UFF'],mixing_rule=mixing_rule)    

    for i in Data:
      Atom_type_dict,fixed_modes = Write_CHARMM(Filename,Data[i]['Atom_types'],Data[i]['Elements'],Data[i]['Geometry'],Data[i]['Bonds'],Data[i]['Bond_types'],Data[i]['Bond_params'],Data[i]['Angles'],Data[i]['Angle_types'],Data[i]['Angle_params'],\
                                            Data[i]['Dihedrals'],Data[i]['Dihedral_types'],Data[i]['Dihedral_params'],Data[i]['Impropers'],Data[i]['Improper_types'],Data[i]['Improper_params'],Data[i]['Charges'],VDW_params,Data[i]['Masses'],alpha_dict,Data[i]['LP_params'],option['inpname'],Improper_flag = improper_flag,polar_flag=polar_flag,LP_flag=lp_flag,overwrite=option['overwrite'],fit_flag=option['fit_flag'],nonC_flag=option['nonC_flag'],npert=option['npert'],ascale=option['ascale'],mag_c=option['mag_c'],onefour_scale=onefour_scale)

    # Remove concatenated FF file from run folder
    os.remove(FF_all)

    # Print banner
    print("\n{}\n* {:^173s} *\n{}\n".format("*"*177,"Success! Have a Nice Day!","*"*177))

    return

def Write_CHARMM(Filename,Atom_types,Elements,Geometry,Bonds,Bond_types,Bond_params,Angles,Angle_types,Angle_params,Dihedrals,Dihedral_types,Dihedral_params,\
               Impropers,Improper_types,Improper_params,Charges,VDW_params,Masses,alpha_dict,LP_params,inpname,Improper_flag=False,polar_flag=False,LP_flag=False,overwrite=False,fit_flag=True,nonC_flag=False,npert=2,ascale=1,mag_c=0.5,onefour_scale=0.0):

    # Write an xyz for easy viewing
    with open(Filename+'/'+inpname+'.xyz','w') as f:
        f.write('{}\n\n'.format(len(Geometry)))
        for count_i,i in enumerate(Geometry):
            f.write('{:20s} {:< 20.6f} {:< 20.6f} {:< 20.6f}\n'.format(Elements[count_i],i[0],i[1],i[2]))


    # Assign ahp Miller atomic polarizability to heavy atoms
    adj = Table_generator(Elements,Geometry)
    hybrid = Hybridization_finder(Elements,adj,force_opt=True)
    # one bond depth atom type (this makes making neighbor list very easy)
    atomT = id_types(Elements,adj,1)


    # See 2005 Anisimov et. for names' definition
    ahp_miller = {"CTE+4H":2.609, "CTE+3H": 2.222,"CTE+2H": 1.835, "CTE+H":1.448,"CTE":1.061,"CTR":1.352,"CTR+H":1.739,\
                  "CTR+2H": 2.126, "CBR": 1.896, "OTE+H": 1.024, "OTE": 0.637, "OTR4":0.569,\
                  "STE+H": 3.387, "STE": 3.0, "STR4":3.729,\
                  "OTA": 0.858,"NTE+3H":2.125,"NTE+2H":1.738,"NTE+H":1.351,"NTE":0.964,"NPI2+2H":1.864,\
                  "NPI2+H":1.477,"NTR2":1.030,"PTE":2.063,"Cl":2.315,"Br":3.013,'CDI':1.283,\
                  "CDI+H":1.67,'NDI':0.956,'NDI+H':1.343}

    atom_polar = { i:'' for i in Atom_types}
    for count_i,i in enumerate(atomT):
         if Elements[count_i] == 'C':
            # neighbor list (but this also includes element itself)
            nei=[ j.split(']')[0] for j in i.split('[') if j.split(']')[0] != '']
            nei_H = nei[1:].count('1')
            if hybrid[count_i] == 'sp3':
               if nei_H == 4: atom_polar[Atom_types[count_i]] = 'CTE+4H'
               if nei_H == 3: atom_polar[Atom_types[count_i]] = 'CTE+3H'
               if nei_H == 2: atom_polar[Atom_types[count_i]] = 'CTE+2H'
               if nei_H == 1: atom_polar[Atom_types[count_i]] = 'CTE+H'
               if nei_H == 0: atom_polar[Atom_types[count_i]] = 'CTE'
            if hybrid[count_i] == 'sp2':
               if nei_H == 2: atom_polar[Atom_types[count_i]] = 'CTR+2H'
               if nei_H == 1: atom_polar[Atom_types[count_i]] = 'CTR+H'
               if nei_H == 0: atom_polar[Atom_types[count_i]] = 'CTR'
            if hybrid[count_i] == 'sp':
               if nei_H == 1: atom_polar[Atom_types[count_i]] = 'CDI+H'
               if nei_H == 0: atom_polar[Atom_types[count_i]] = 'CDI'
         if Elements[count_i] == 'O':
            # neighbor list (but this also includes element itself)
            nei=[ j.split(']')[0] for j in i.split('[') if j.split(']')[0] != '']
            nei_H = nei[1:].count('1')
            if hybrid[count_i] == 'sp3':
               if nei_H == 1: atom_polar[Atom_types[count_i]] = 'OTE+H'
               else: atom_polar[Atom_types[count_i]] = 'OTE'
            if hybrid[count_i] == 'sp2':
               atom_polar[Atom_types[count_i]] = 'OTR4'
         if Elements[count_i] == 'S':
            # neighbor list (but this also includes element itself)
            nei=[ j.split(']')[0] for j in i.split('[') if j.split(']')[0] != '']
            nei_H = nei[1:].count('1')
            if hybrid[count_i] == 'sp3':
               if nei_H == 1: atom_polar[Atom_types[count_i]] = 'STE+H'
               else: atom_polar[Atom_types[count_i]] = 'STE'
            if hybrid[count_i] == 'sp2':
               atom_polar[Atom_types[count_i]] = 'STR4'
         if Elements[count_i] == 'N':
            # neighbor list (but this also includes element itself)
            nei=[ j.split(']')[0] for j in i.split('[') if j.split(']')[0] != '']
            nei_H = nei[1:].count('1')
            if hybrid[count_i] == 'sp3':
               if nei_H == 3: atom_polar[Atom_types[count_i]] = 'NTE+3H'
               if nei_H == 2: atom_polar[Atom_types[count_i]] = 'NTE+2H'
               if nei_H == 1: atom_polar[Atom_types[count_i]] = 'NTE+H'
               if nei_H == 0: atom_polar[Atom_types[count_i]] = 'NTE'
            if hybrid[count_i] == 'sp2':
               if nei_H == 2: atom_polar[Atom_types[count_i]] = 'NPI2+2H'
               elif nei_H == 1: atom_polar[Atom_types[count_i]] = 'NPI2+H'
               elif nei_H == 0 and (len(nei)-1)==2: atom_polar[Atom_types[count_i]] = 'NTR2'
            if hybrid[count_i] == 'sp':
               if nei_H == 1: atom_polar[Atom_types[count_i]] = 'NDI+H'
               elif nei_H == 0: atom_polar[Atom_types[count_i]] = 'NDI'
         if Elements[count_i] == 'Cl':
               atom_polar[Atom_types[count_i]] = 'Cl'
         if Elements[count_i] == 'Br':
               atom_polar[Atom_types[count_i]] = 'Br'
               
                  
    # Create type dictionaries (needed to convert each atom,bond,angle, and dihedral type to consecutive numbers as per LAMMPS convention)
    # Note: LAMMPS orders atomtypes, bonds, angles, dihedrals, etc as integer types. Each of the following dictionaries holds the mapping between
    #       the true type (held in the various types lists) and the lammps type_id, obtained by enumerated iteration over the respective set(types).
    Atom_type_dict = {}
    for count_i,i in enumerate(sorted(set(Atom_types))):
        for j in Atom_types:
            if i == j:
                Atom_type_dict[i]=count_i+1
            if i in list(Atom_type_dict.keys()):
                break
    Bond_type_dict = {}
    for count_i,i in enumerate(sorted(set(Bond_types))):
        for j in Bond_types:
            if i == j:
                Bond_type_dict[i]=count_i+1
            if i in list(Bond_type_dict.keys()):
                break
    Angle_type_dict = {}
    for count_i,i in enumerate(sorted(set(Angle_types))):
        for j in Angle_types:
            if i == j:
                Angle_type_dict[i]=count_i+1
            if i in list(Angle_type_dict.keys()):
                break
    Dihedral_type_dict = {}
    for count_i,i in enumerate(sorted(set(Dihedral_types))):
        for j in Dihedral_types:
            if i == j:
                Dihedral_type_dict[i]=count_i+1
            if i in list(Dihedral_type_dict.keys()):
                break
    Improper_type_dict = {}
    for count_i,i in enumerate(sorted(set(Improper_types))):
        for j in Improper_types:
            if i == j:
                Improper_type_dict[i]=count_i+1
            if i in list(Improper_type_dict.keys()):
                break

    with open(Filename+'/'+'lin.rtf','w') as g:
        # Write head
        g.write("* Lin topology\n*\n") 
        g.write("41  \n\n") # Version number

        # Atom_type_dict is TAFFI to lammps atomtype
        # Atom_type_dict2 is Lammps to CHARMM atomtype
        Atom_type_dict2 = {} 
        # Write Masses
        for count_i,i in enumerate(sorted(set(Atom_types))):
            for j in set(Atom_types):
                if Atom_type_dict[j] == count_i+1:
                    element_tmp = Num2Element(int(i.split('[')[1]))
                    g.write("MASS {}  {}R{}  {:< 8.6f} {}\n".format(-1,element_tmp,count_i+1,Masses[str(j)],element_tmp)) # count_i+1 bc of LAMMPS 1-indexing
                    Atom_type_dict2[count_i+1] = "{}R{}".format(element_tmp,count_i+1,element_tmp)


        # Add drude particle, drude particle's name is based on atomname not atomtype name
        N_atomtype = len(set(Atom_types))
        Atom_name_dict = {} # Lammps atom name to CHARMM atom name
        
        for count_i,i in enumerate(Atom_types):
            element_tmp = Num2Element(int(i.split('[')[1]))
            Atom_name_dict[count_i+1] = "{}{:<4d}".format(element_tmp,count_i+1)
            #if element_tmp != "H":
            #    g.write("MASS {}  D{}  {:< 8.6f}\n".format(N_atomtype+1,Atom_name_dict[count_i+1],0.1)) # count_i+1 bc of LAMMPS 1-indexing
            #    N_atomtype += 1
        
        g.write("MASS {}  DRUD  {:<8.6f} H ! Drude particle\n".format(-1,0.0))
        g.write("MASS {}  LPD   {:<8.6f} H ! Lone pairs\n".format(-1,0.0))
        g.write("MASS {}  CAL   {:<8.6f} CA\n".format(-1,40.08))
               
        g.write("\n")
        g.write("default first none last none\n")
        if polar_flag or fit_flag:
            g.write("AUTOGENERATE DRUDE\n")
        g.write("RESI LIN 0.00 \n")
        g.write("GROUP\n")
        # Write Atoms
        N_drude = 0
        source_tag = {} # source tag indicate where polarizability comes from, FF means it's from FF.db file, ahp means from ahp_miller
        #if np.sum(np.array(Charges)) > 0.1:
        #    print("Sum of charge larger than 0.1, exiting...")
        #    quit()

        # XXX this should be combined with db_to_charmm.py
        for count_i,i in enumerate(Atom_types):
            element_tmp = Num2Element(int(i.split('[')[1]))
            if polar_flag or fit_flag:
               if element_tmp != "H":
                   N_drude += 1
                   if nonC_flag and element_tmp == 'C':
                      g.write("ATOM {:<6s} {:<6s} {:< 20.16f} ! {}\n"\
                        .format(Atom_name_dict[count_i+1],Atom_type_dict2[Atom_type_dict[i]],Charges[count_i],i))
                      continue
                   if alpha_dict[i] > 0.001: # 0.01 because in FF.db Hydrogen's polarizability is set to 0, just in case             
                        g.write("ATOM {:<6s} {:<6s} {:< 20.16f} ALPHA {:< 3.6f} THOLE 1.3 ! {}\n"\
                        .format(Atom_name_dict[count_i+1],Atom_type_dict2[Atom_type_dict[i]],Charges[count_i],alpha_dict[i]*(-1),i))
                        source_tag[i] = 'FF'
                   else:
                        if overwrite:
                           g.write("ATOM {:<6s} {:<6s} {:< 20.16f} ! {}\n"\
                           .format(Atom_name_dict[count_i+1],Atom_type_dict2[Atom_type_dict[i]],Charges[count_i],i))
                        if fit_flag:
                           if atom_polar[i] =='':
                              print("ERROR: missing ahp miller for {} missing. Exiting...".format(i))
                              quit()
                           else:
                              g.write("ATOM {:<6s} {:<6s} {:< 20.16f} ALPHA {:< 3.6f} THOLE 1.3\n"\
                                 .format(Atom_name_dict[count_i+1],Atom_type_dict2[Atom_type_dict[i]],Charges[count_i],ahp_miller[atom_polar[i]]*(-1)))
                        else:
                           print("ERROR: atomic polarizability for {} missing. Exiting....".format(i))
                           quit()
               else:
                   source_tag[i] = 'H'
                   g.write("ATOM {:<6s} {:<6s} {:< 20.16f} ! {}\n"\
                     .format(Atom_name_dict[count_i+1],Atom_type_dict2[Atom_type_dict[i]],Charges[count_i],i))
            else:
               g.write("ATOM {:<6s} {:<6s} {:< 20.16f} ! {}\n"\
                  .format(Atom_name_dict[count_i+1],Atom_type_dict2[Atom_type_dict[i]],Charges[count_i],i))
        if LP_flag:
            for i in LP_params:
               if len(i.split('_')) == 1:
                    g.write("ATOM {:<6s} {:<6s} {:< 20.16f}\n".format(i,'LPD',LP_params[i]))
        g.write("\n")
        # Write Bonds
        if len(Bonds) > 0:
            for count_i,i in enumerate(Bonds):
                g.write("BOND {} {} \n".format(Atom_name_dict[i[0]+1],Atom_name_dict[i[1]+1]))
            g.write("\n")

        # Write Angles
        if len(Angles) > 0:
            for count_i,i in enumerate(Angles):
                g.write("ANGL {} {} {} \n".format(Atom_name_dict[i[0]+1],Atom_name_dict[i[1]+1],Atom_name_dict[i[2]+1]))
            g.write("\n")

        # Write Dihedrals
        if len(Dihedrals) > 0: 
            for count_i,i in enumerate(Dihedrals):
                g.write("DIHE {} {} {} {} \n".format(Atom_name_dict[i[0]+1],Atom_name_dict[i[1]+1],Atom_name_dict[i[2]+1],Atom_name_dict[i[3]+1]))
            g.write("\n")

        # Write Impropers
        if Improper_flag and len(Impropers) > 0: 
            for count_i,i in enumerate(Impropers):
                g.write("IMPH {} {} {} {} \n".format(Atom_name_dict[i[0]+1],Atom_name_dict[i[1]+1],Atom_name_dict[i[2]+1],Atom_name_dict[i[3]+1]))
            g.write("\n")

        # Write LP
        if LP_flag:
            for i in LP_params:
               if len(i.split('_')) == 1:
                  if LP_params[i+'_dihedral'] < 0: LP_params[i+'_dihedral'] = LP_params[i+'_dihedral']+360
                  g.write("LONEPAIR relative {:<4s} O2 C1 H3 distance {:<3.2f} angle {:<3.2f} dihe {:<3.2f}\n".format(i,LP_params[i+'_dist'],LP_params[i+'_angle'],LP_params[i+'_dihedral']))

        # Write patch
        g.write("PATCH FIRST NONE LAST NONE\n\n")

        # Wrtie perturnbation charge
        g.write("RESI CAL    2.00\n")
        g.write("GROUP \n") 
        g.write("ATOM CAL  CAL 2.00\n")
        g.write("PATCH FIRST NONE LAST NONE\n\n")


        # end of topology file
        g.write("end\n\n")
        

        # Write parameter data
        fixed_modes = {'bonds':[], 'angles':[]}

    with open(Filename+'/'+'lin.prm','w') as g:
        g.write("* Lin parameters\n*\n")
        

        # Write stretching interactions
        # Note: Bond_type_dict was initialized by looping over sorted(set(Bond_types)), so 
        #       iterating over this list resulted in ordered parameters in the in.settings file
        g.write("BONDS\n")
        for i in sorted(set(Bond_types)):
            g.write("{}  {} ".format(Atom_type_dict2[Atom_type_dict[i[0]]],Atom_type_dict2[Atom_type_dict[i[1]]]))
            for j in Bond_params[i]:
                if j == "fixed":
                    continue
                if type(j) is str:
                    g.write("{:20s} ".format(j))
                if type(j) is float:
                    g.write("{:< 20.6f} ".format(j))
            g.write("\n")

            # populate fixed_modes
            if Bond_params[i][0] == "fixed":
                fixed_modes["bonds"] += [Bond_type_dict[i]]

        g.write("X  DRUD  500.0   0\n")
        """
        for count_i,i in enumerate(Atom_types):
            element_tmp = Num2Element(int(i.split('[')[1]))
            if element_tmp != "H":
                g.write("{}  D{}  500.0   0\n".format(Atom_type_dict2[Atom_type_dict[i]],Atom_name_dict[count_i+1])) 
        """
        g.write("\n")
        # Write bending interactions
        # Note: Angle_type_dict was initialized by looping over sorted(set(Angle_types)), so 
        #       iterating over this list resulted in ordered parameters in the in.settings file
        g.write("ANGLES\n")
        for i in sorted(set(Angle_types)):
            g.write("{}  {}  {}  ".format(Atom_type_dict2[Atom_type_dict[i[0]]],Atom_type_dict2[Atom_type_dict[i[1]]],Atom_type_dict2[Atom_type_dict[i[2]]]))
            for j in Angle_params[i]:
                if j == "fixed":
                    continue
                if type(j) is str:
                    g.write("{:20s} ".format(j))
                if type(j) is float:
                    g.write("{:< 20.6f} ".format(j))
            g.write("\n")

            # populate fixed_modes
            if Angle_params[i][0] == "fixed":
                fixed_modes["angles"] += [Angle_type_dict[i]]

        # Write dihedral interactions
        # Note: Dihedral_type_dict was initialized by looping over sorted(set(Dihedral_types)), so 
        #       iterating over this list resulted in ordered parameters in the in.settings file
        # OPLS dihedral style: 1/2*K1*(1+cos(phi)) + 1/2*K2*(1-cos(2phi)) + 1/2*K3*(1+cos(3phi)) + 1/2*K4*(1-cos(4phi))
        # Lammps harmonic dihedral type: K[1+dcos(n phi)]
        # CHARMM dihedrla style: K*(1+cos(n*phi-d))
        # We transform opls to charmm by specifying multliplicity for the same dihedral, and use phase shift to accont +-

        g.write("\nDIHEDRALS\n")
        for i in sorted(set(Dihedral_types)):
            opls_flag = True
            for count_j,j in enumerate(Dihedral_params[i]):
                if count_j == 0  and j == 'harmonic':
                    opls_flag = False
                if opls_flag:
                   if count_j == 1:
                       g.write("{}  {}  {}  {}  ".format(Atom_type_dict2[Atom_type_dict[i[0]]],Atom_type_dict2[Atom_type_dict[i[1]]],Atom_type_dict2[Atom_type_dict[i[2]]],Atom_type_dict2[Atom_type_dict[i[3]]]))
                       g.write("{:< 20.6f} 1 0.00".format(j*0.5))
                       g.write("\n")
                   if count_j == 2:
                       g.write("{}  {}  {}  {}  ".format(Atom_type_dict2[Atom_type_dict[i[0]]],Atom_type_dict2[Atom_type_dict[i[1]]],Atom_type_dict2[Atom_type_dict[i[2]]],Atom_type_dict2[Atom_type_dict[i[3]]]))
                       g.write("{:< 20.6f} 2 180.0".format(j*0.5))
                       g.write("\n")
                   if count_j == 3:
                       g.write("{}  {}  {}  {}  ".format(Atom_type_dict2[Atom_type_dict[i[0]]],Atom_type_dict2[Atom_type_dict[i[1]]],Atom_type_dict2[Atom_type_dict[i[2]]],Atom_type_dict2[Atom_type_dict[i[3]]]))
                       g.write("{:< 20.6f} 3 0.00".format(j*0.5))
                       g.write("\n")
                   if count_j == 4:
                       g.write("{}  {}  {}  {}  ".format(Atom_type_dict2[Atom_type_dict[i[0]]],Atom_type_dict2[Atom_type_dict[i[1]]],Atom_type_dict2[Atom_type_dict[i[2]]],Atom_type_dict2[Atom_type_dict[i[3]]]))
                       g.write("{:< 20.6f} 4 180.0".format(j*0.5))
                       g.write("\n")
                else:
                   if count_j == 1:
                       g.write("{}  {}  {}  {}  ".format(Atom_type_dict2[Atom_type_dict[i[0]]],Atom_type_dict2[Atom_type_dict[i[1]]],Atom_type_dict2[Atom_type_dict[i[2]]],Atom_type_dict2[Atom_type_dict[i[3]]]))
                       g.write("{:< 20.6f} 2 180.0".format(j))
                       g.write("\n")
                   if count_j == 2:
                       if j != -1:
                           print('ERROR: harmonic dihedral not equal K(1-cos(2phi))')
                           quit()
                   if count_j == 3:
                       if j != 2:
                           print('ERROR: harmonic dihedral not equal K(1-cos(2phi))')
                           quit()

            #g.write("{} 3   0.00".format(Dihedral_params[i][3]*0.5))
        g.write("\n")

        g.write("NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vshift -\n")
        g.write("cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac {:3.2f} wmin 1.5\n".format(onefour_scale))

        # Write non-bonded interactions (the complicated form of this loop is owed
        # to desire to form a nicely sorted list in terms of the lammps atom types
        # Note: Atom_type_dict was initialize according to sorted(set(Atom_types) 
        #       so iterating over this list (twice) orders the pairs, too.
        for count_i,i in enumerate(sorted(set(Atom_types))):     
            eps = VDW_params[(i,i)][1] # kcal to kJ
            Rmin = VDW_params[(i,i)][2]*(2**(1.0/6)) # A to nm 
            g.write("{}  0.000000 {:<20.6f}   {:<20.6f} 0.00000  {:<20.6f}   {:<20.6f}\n".format(Atom_type_dict2[Atom_type_dict[i]],eps*(-1),Rmin/2,eps*(-1)*onefour_scale,Rmin/2))
        g.write("DRUD  0.0  -0.0 0.01\n")
        g.write("LPD   0.0  -0.0 0.01\n")
        g.write("CAL   0.0  -0.12 1.71\n")

        # Wrtie crossed terms, in CHARMM NBFIX specified eps and Rmin but in the nonbonded section eps and Rmin/2 are specified
        g.write("\nNBFIX\n")
        print(Atom_types)
        for count_i,i in enumerate(sorted(set(Atom_types))):     
            for count_j,j in enumerate(sorted(set(Atom_types))): 

                # Skip duplicates
                if count_j <= count_i:
                    continue

                # Conform to LAMMPS i <= j formatting
                if Atom_type_dict[i] <= Atom_type_dict[j]:
                    g.write("{:20s} {:20s} ".format(Atom_type_dict2[Atom_type_dict[i]],Atom_type_dict2[Atom_type_dict[j]]))
                else:
                    g.write("{:20s} {:20s} ".format(Atom_type_dict2[Atom_type_dict[j]],Atom_type_dict2[Atom_type_dict[i]]))

                # Determine key (ordered by initialize_VDW function such that i > j)
                if i > j:
                    key = (i,j)
                else:
                    key = (j,i)

                # Write the parameters
                eps = VDW_params[key][1] # kcal to kJ
                Rmin = VDW_params[key][2]*(2**(1.0/6)) # A to nm 
                g.write("{:< 8.6f} {:< 8.6f} ".format(eps*(-1),Rmin))
                # 1-4 interaction
                g.write("{:< 8.6f} {:< 8.6f} \n".format(eps*(-1)*onefour_scale,Rmin))
       
        """
        for count_i,i in enumerate(Atom_types):
            element_tmp = Num2Element(int(i.split('[')[1]))
            if element_tmp != "H":
                g.write("D{}  0.00 -0.000  0.00 \n".format(Atom_name_dict[count_i+1])) 
        """
        g.write("\n\nend\n\n")

    # Write the data file
    with open(Filename+'/'+inpname+'.inp','w') as f:
       
        f.write("* Charge fitting of polarizable molecule\n")
        f.write("* Data Files: fitcharge.pot, fitcharge.pot0, fitcharge.qpos\n\n")
        f.write("stream ../datadir_{}.def\n".format(mag_c))
        f.write("! this is necessart to work with non-interger total charge\n")
        f.write("bomlev -1\n\n")
        f.write("read rtf card name lin.rtf\n")
        f.write("read para card name lin.prm\n")
        f.write("read sequence card\n")
        f.write("* lin\n")
        f.write("*\n")
        f.write("1\n")
        f.write("lin\n\n")
        f.write("generate lin first none last none setup warn drude dmass 0.4\n\n")

        f.write("read sequence card\n") 
        f.write("* Calcium\n")
        f.write("*\n")
        f.write("1\n")
        f.write("CAL\n\n")

        f.write("generate CAl warn setup\n")
        f.write("!assign charge\n")
        f.write("scalar charge set {} select segid CAL end\n\n".format(mag_c))

        f.write("read coor card free\n")
        f.write("* lin.qpos.1.crd\n")
        f.write("*\n")
        f.write("   {}\n".format(len(Atom_types)+1)) 
        atom_per_mol = len(Atom_types)
        for count_i,i in enumerate(Atom_types):
           element_tmp = Num2Element(int(i.split('[')[1]))
           f.write("{:<4d} {:< 4d} LIN {:} {:< 10.6f} {:< 10.6f} {:< 10.6f} LIN {:d} {:7.5f}\n"\
            .format(count_i+1,1,Atom_name_dict[count_i+1],Geometry[count_i,0],Geometry[count_i,1],Geometry[count_i,2],1,0.00000))
        f.write("{:<4d} {:< 4d} CAL {:} {:< 10.6f} {:< 10.6f} {:< 10.6f} CAL {:d} {:7.5f}\n"\
            .format(len(Atom_types)+1,2,'CAL',Geometry[count_i,0],Geometry[count_i,1],Geometry[count_i,2],2,0.00000))

        f.write("\n\n")
        f.write("coor sdrude\n")
        f.write("coor shake\n")

        f.write("! before fitting \n") 
        f.write("scalat charge show\n")
        
        f.write("! print out final coor (drude included)\n")
        f.write("coor print\n")
        f.write("coor dipole oxyz select resn lin end\n\n")

        """
        f.write("! set polarizabilities on non-hydrogen atoms\n")
        # temporary polarizability
       
        for count_i,i in enumerate(Atom_types):
            element_tmp = Num2Element(int(i.split('[')[1]))
            if element_tmp != "H":
                if atom_polar[i] == '':
                     print("WARNING: missing ahp miller atomic prolarizability, use 0.5 A^3 as initial guess for atomtype:{}".format(i))
                     f.write("scalar wmain  set  -{}  select resn lin .and. type {}   END\n".format(0.5,Atom_name_dict[count_i+1]))
                else:     
                     f.write("scalar wmain  set  -{}  select resn lin .and. type {}   END\n".format(ahp_miller[atom_polar[i]],Atom_name_dict[count_i+1]))
        f.write("\n\n")

        f.write("! generate the Drude particles\n")
        f.write("DRUDE thole 2.6 mass 0.1 kdrude 500 -\n")
        f.write("  select .not. ( type H* .or. type CAL ) end\n\n")
        """

        f.write("update inbfrq -1 ihbfrq 0 -\n")
        f.write("  switch atom vatom vswitch cutnb 990. ctofnb 950. ctonnb 900.\n\n")

        f.write("! CHARGE AND DIPOLE BEFORE FITTING\n")
        f.write("scalar charge show\n")
        f.write("coor print\n")
        f.write("coor dipole oxyz select resn lin end\n\n\n")


        f.write("set fname fitcharge\n\n")

        f.write("! potential for unperturbed system\n")
        f.write("open read  card unit 11 name @0@fname.pot0\n\n")
        f.write("! potential for perturbed systems\n")
        f.write("open read  card unit 21 name @0@fname.pot\n\n")
        f.write("! all the positions of the 0.5 charge\n")
        f.write("open read  card unit 31 name @0@fname.qpos\n\n")
        f.write("! scratch file\n")
        f.write("open write card unit 30 name @9@fname.scr\n\n")
        f.write("prnlev 3\n\n")
        f.write("! set restraints\n")
        f.write("scalar wmain set 1.0d-5 select segid lin end  ! put a mild restraint on all atomic charges\n")
        already = [ "type {} .or. type D{}".format(Atom_name_dict[count_j+1],Atom_name_dict[count_j+1]) for count_j,j in enumerate(Atom_types)  if alpha_dict[j] != -1 ]
        if len(already) > 1:
            if len(already)>2:
               #split to multiple lines or CHARMM cannot handle it
               n=2
               split_list = [already[i * n:(i + 1) * n] for i in range((len(already) + n - 1) // n )]
               for count_i,i in enumerate(split_list):
                  if count_i == 0:
                     f.write("define excluded select .not. ( type CAL .or. {} -\n".format(" .or. ".join(i)))
                  elif count_i == len(split_list)-1:
                     f.write(".or. {} ) end\n".format(" .or. ".join(i)))
                  else:
                     f.write(".or. {} -\n".format(" .or. ".join(i)))
            else:
               f.write("define excluded select .not. ( type CAL .or. {} ) end\n".format(" .or. ".join(already)))
        else:
            f.write("define excluded select .not. ( type CAL ) end\n")
        if len(already) == len(Atom_types) : 
            print("No atoms in need of parametrized! Exiting....")
            quit()
            
        f.write("define source select .not. ( type CAL ) end\n")
        f.write("scalar wmain show\n")

        f.write("! enable this unit to have the fitted charges saved in a stream file\n")
        f.write("open write card unit 90 name lin.charge.optimized\n\n")
        f.write("FITCHARGE -\n")


        # define equivalent atom
        for count_i,i in enumerate(set(Atom_types)):
            equiv = [ "type {}".format(Atom_name_dict[count_j+1]) for count_j,j in enumerate(Atom_types) if j==i]
            if(len(equiv) >1):
               #f.write("    equivalent select ( {} ) end - \n".format(" .or. ".join(equiv)))
               n=3
               split_list = [equiv[j * n:(j + 1) * n] for j in range((len(equiv) + n - 1) // n )]
               for count_j,j in enumerate(split_list):
                  if count_j == 0:
                     if len(split_list) == 1:
                        f.write("    equivalent select ( {} ) end - \n".format(" .or. ".join(j)))
                     else:
                        f.write("    equivalent select ( {} - \n".format(" .or. ".join(j)))
                  elif count_j == len(split_list)-1:
                     f.write(".or. {} ) end -\n".format(" .or. ".join(j)))
                  else:
                     f.write(".or. {} -\n".format(" .or. ".join(j)))
        f.write("    select excluded end - ! atoms to fit\n") 
        f.write("    select source end -  ! ESP source\n")
        f.write("    restraint resp para - \n")
        f.write("    flat 0.0  dflat 0.1 -\n")
        f.write("    upot 11 uout 30 -\n")
        f.write("    NITE 50 -\n")
        f.write("    NCONF 1 -\n")
        f.write("    NPERT {} -             ! number of perturbation ions used in QM ESP\n".format(npert))
        f.write("    ascale {} -\n".format(ascale))
        f.write("    uppot 21 -\n")
        f.write("    ucoord 31 altinput -\n")
        f.write("    uprt 90\n")
        f.write("!delete scratch file\n")
        f.write("close unit 30 disposition delete\n\n")
        f.write("! show final fitting charge\n") 
        f.write("scalat charge show\n")
        
        f.write("! print out final coor (drude included)\n")
        f.write("coor print\n")
        f.write("coor dipole oxyz select resn lin end\n\n")

        f.write("! calculate molecular polarizability\n")
        f.write("VIBRAN\n")
        f.write("DIAGONALIZE\n")
        f.write("PRINT NORMAL  DIPOLES SELECT excluded END\n")
        f.write("END\n")


    return Atom_type_dict,fixed_modes



def Num2Element(num):
    # Initialize periodic table

    periodic = {1: 'H', 2: 'HE', 3: 'LI', 4: 'BE', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'NE', 11: 'NA', 12: 'MG', 13: 'AL', 14: 'SI', 15: 'P', 16: 'S', 17: 'CL', 18: 'AR', 19: 'K', 20: 'CA', 21: 'SC', 22: 'TI', 23: 'V', 24: 'CR', 25: 'MN', 26: 'FE', 27: 'CO', 28: 'NI', 29: 'CU', 30: 'ZN', 31: 'GA', 32: 'GE', 33: 'AS', 34: 'SE', 35: 'BR', 36: 'KR', 37: 'RB', 38: 'SR', 39: 'Y', 40: 'ZR', 41: 'NB', 42: 'MO', 43: 'TC', 44: 'RU', 45: 'RH', 46: 'PD', 47: 'AG', 48: 'CD', 49: 'IN', 50: 'SN', 51: 'SB', 52: 'TE', 53: 'I', 54: 'XE', 55: 'CS', 56: 'BA', 72: 'HF', 73: 'TA', 74: 'W', 75: 'RE', 76: 'OS', 77: 'IR', 78: 'PT', 79: 'AU', 80: 'HG', 81: 'TL', 82: 'PB', 83: 'BI', 84: 'PO', 85: 'AT', 86: 'RN'}

    element = periodic[num]
    return element



# A wrapper for the commands to parse the bonds, angles, and dihedrals from the adjacency matrix.
# Returns:   list of atomtypes, bond_types, bond instances, angle_types, angle instances, dihedral_types,
#            diehdral instances, charges, and VDW parameters.
def Find_parameters(Adj_mat,Bond_mats,Geometry,Atom_types,FF_db="FF_file",Improper_flag = False, force_read=False,remove_multi=False):

    # List comprehension to determine bonds from a loop over the adjacency matrix. Iterates over rows (i) and individual elements
    # ( elements A[count_i,count_j] = j ) and stores the bond if the element is "1". The count_i < count_j condition avoids
    # redudant bonds (e.g., (i,j) vs (j,i) ). By convention only the i < j definition is stored.
    print("Parsing bonds...")
    Bonds          = [ (count_i,count_j) for count_i,i in enumerate(Adj_mat) for count_j,j in enumerate(i) if j == 1 ]
    Bond_types     = [ (Atom_types[i[0]],Atom_types[i[1]]) for i in Bonds ]

    # List comprehension to determine angles from a loop over the bonds. Note, since there are two bonds in every angle, there will be
    # redundant angles stored (e.g., (i,j,k) vs (k,j,i) ). By convention only the i < k definition is stored.
    print("Parsing angles...")
    Angles          = [ (count_j,i[0],i[1]) for i in Bonds for count_j,j in enumerate(Adj_mat[i[0]]) if j == 1 and count_j not in i ]
    Angle_types     = [ (Atom_types[i[0]],Atom_types[i[1]],Atom_types[i[2]]) for i in Angles ]

    # List comprehension to determine dihedrals from a loop over the angles. Note, since there are two angles in every dihedral, there will be
    # redundant dihedrals stored (e.g., (i,j,k,m) vs (m,k,j,i) ). By convention only the i < m definition is stored.
    print("Parsing dihedrals...")
    Dihedrals      = [ (count_j,i[0],i[1],i[2]) for i in Angles for count_j,j in enumerate(Adj_mat[i[0]]) if j == 1 and count_j not in i ]
    Dihedral_types = [ (Atom_types[i[0]],Atom_types[i[1]],Atom_types[i[2]],Atom_types[i[3]]) for i in Dihedrals ]

    # List comprehension to determine dihedrals from a loop over the angles. Note, since there are two angles in every dihedral, there will be
    # redundant dihedrals stored (e.g., (i,j,k,m) vs (m,k,j,i) ). By convention only the i < m definition is stored.
    if Improper_flag:    print("Parsing impropers...")
    Impropers      = [ (i[1],i[0],i[2],count_j) for i in Angles for count_j,j in enumerate(Adj_mat[i[1]]) if j == 1 and count_j not in i ]
    Improper_types = [ (Atom_types[i[0]],Atom_types[i[1]],Atom_types[i[2]],Atom_types[i[3]]) for i in Impropers ]

    # Canonicalize the modes
    for i in range(len(Bonds)):
        Bond_types[i],Bonds[i] = canon_bond(Bond_types[i],ind=Bonds[i])
    for i in range(len(Angles)):
        Angle_types[i],Angles[i] = canon_angle(Angle_types[i],ind=Angles[i])
    for i in range(len(Dihedrals)):
        Dihedral_types[i],Dihedrals[i] = canon_dihedral(Dihedral_types[i],ind=Dihedrals[i])
    for i in range(len(Impropers)):
        Improper_types[i],Impropers[i] = canon_improper(Improper_types[i],ind=Impropers[i])        

    # Remove redundancies    
    if len(Bonds) > 0: Bonds,Bond_types = list(map(list, list(zip(*[ (i,Bond_types[count_i]) for count_i,i in enumerate(Bonds) if count_i == [ count_j for count_j,j in enumerate(Bonds) if j == i or j[::-1] == i ][0]  ]))))
    if len(Angles) > 0: Angles,Angle_types = list(map(list, list(zip(*[ (i,Angle_types[count_i]) for count_i,i in enumerate(Angles) if count_i == [ count_j for count_j,j in enumerate(Angles) if j == i or j[::-1] == i ][0]  ]))))
    if len(Dihedrals) > 0: Dihedrals,Dihedral_types = list(map(list, list(zip(*[ (i,Dihedral_types[count_i]) for count_i,i in enumerate(Dihedrals) if count_i == [ count_j for count_j,j in enumerate(Dihedrals) if j == i or j[::-1] == i ][0]  ]))))
    if len(Impropers) > 0: Impropers,Improper_types = list(map(list, list(zip(*[ (i,Improper_types[count_i]) for count_i,i in enumerate(Impropers) if count_i == [ count_j for count_j,j in enumerate(Impropers) if j[0] == i[0] and len(set(i[1:]).intersection(set(j[1:]))) ][0] ]))))

    # Remove bimodal angles (e.g., if the same angle is present in
    # more than one configuration, then only the most abundant is kept)
    if remove_multi is True:
        for i in set(Angle_types):

            # Calculate most populated bin
            bins = [ 20.0*_ for _ in arange(9) ]
            hist = zeros(len(bins))
            angles = []
            for count_j,j in enumerate(Angle_types):
                if j == i:
                    v1 = Geometry[Angles[count_j][0]] - Geometry[Angles[count_j][1]]
                    v2 = Geometry[Angles[count_j][2]] - Geometry[Angles[count_j][1]]
                    v1 = v1/norm(v1)
                    v2 = v2/norm(v2)
                    angles += [ acos(dot(v1,v2)) * 180/pi ]
                    hist[int(angles[-1]/20.0)] += 1            
                else:
                    angles += [0.0]
            angle = bins[where(hist == max(hist))[0][0]]

            # Determine the angles to be removed
            del_inds = []
            for count_j,j in enumerate(Angle_types):                
                if j == i and abs(angles[count_j] - angle) > 20.0:
                    del_inds += [count_j]
            Angle_types = [ j for count_j,j in enumerate(Angle_types) if count_j not in del_inds ]
            Angles = [ j for count_j,j in enumerate(Angles) if count_j not in del_inds ]

    # Add the opls/harmonic type as a fifth element in the Dihedral_types tuples
    for count_i,i in enumerate(Dihedrals):
        if 2 in [ j[i[1],i[2]] for j in Bond_mats ]:
            Dihedral_types[count_i] = tuple(list(Dihedral_types[count_i]) + ["harmonic"])        
        # elif ( "R" in Dihedral_types[count_i][1] or "E" in Dihedral_types[count_i][1] or "Z" in Dihedral_types[count_i][1] ) and \
        #      ( "R" in Dihedral_types[count_i][2] or "E" in Dihedral_types[count_i][2] or "Z" in Dihedral_types[count_i][2] ):
        #     Dihedral_types[count_i] = tuple(list(Dihedral_types[count_i]) + ["harmonic"])        
        else:
            Dihedral_types[count_i] = tuple(list(Dihedral_types[count_i]) + ["opls"])                        

    ##############################################################
    # Read in parameters: Here the stretching, bending, dihedral #
    # and non-bonding interaction parameters are read in from    #
    # parameters file. Mass and charge data is also included.    #
    # The program looks for a simple match for the first entry   #
    # in each line with one of the bond or angle types.          #
    # INPUT: param_file, BOND_TYPES_LIST, ANGLE_TYPES_LIST       #
    #        DIHEDRAL_TYPES_LIST, ELEMENTS                       #
    # OUTPUT: CHARGES, MASSES, BOND_PARAMS, ANGLE_PARAMS,        #
    #         DIHERAL_PARAMS, PW_PARAMS                          #
    ##############################################################

    # Initialize dictionaries

    # Read in masses and charges
    Masses = {}
    with open(FF_db,'r') as f:
        content=f.readlines()
        
    for lines in content:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        if fields[0].lower() == 'atom':
            Masses[fields[1]] = float(fields[3]) 

    # Read in bond parameters
    Bond_params = {}
    with open(FF_db,'r') as f:
        content=f.readlines()

    for lines in content:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        if fields[0] == 'bond':
            if fields[3] == "harmonic":
                Bond_params[(fields[1],fields[2])] = [float(fields[4]),float(fields[5])]
            elif fields[3] == "fixed":
                Bond_params[(fields[1],fields[2])] = ["fixed",0.0,float(fields[4])]
            else:
                print("ERROR: only harmonic bond definitions are currently supported by gen_md_for_sampling.py. Exiting...")
                quit()

    # Read in angle parameters
    Angle_params = {}
    with open(FF_db,'r') as f:
        content=f.readlines()

    for lines in content:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        if fields[0].lower() == 'angle':
            if fields[4] == "harmonic":
                Angle_params[(fields[1],fields[2],fields[3])] = [float(fields[5]),float(fields[6])]
            elif fields[4] == "fixed":
                Angle_params[(fields[1],fields[2],fields[3])] = ["fixed",0.0,float(fields[5])]
            else:
                print("ERROR: only harmonic angle definitions are currently supported by gen_md_for_sampling.py. Exiting...")
                quit()

    # Read in dihedral parameters
    Dihedral_params = {}
    with open(FF_db,'r') as f:
        content=f.readlines()

    for lines in content:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        if fields[0].lower() in ['dihedral','torsion']:            
            if fields[5] == "opls":
                Dihedral_params[(fields[1],fields[2],fields[3],fields[4],fields[5])] = [fields[5]] + [ float(i) for i in fields[6:10] ]
            elif fields[5] == "harmonic":
                Dihedral_params[(fields[1],fields[2],fields[3],fields[4],fields[5])] = [fields[5]] + [ float(fields[6]),int(float(fields[7])),int(float(fields[8])) ]
            else:
                print("ERROR: Only opls and harmonic dihedral types are currently supported by gen_md_for_sampling.py. Exiting...")
                quit()
        
    # Read in improper parameters
    Improper_params = {}
    with open(FF_db,'r') as f:
        content=f.readlines()

    for lines in content:
        fields=lines.split()

        # Skip empty lines
        if len(fields) == 0:
            continue

        if fields[0].lower() in ['improper']:
            if fields[5] == "harmonic":
                Improper_params[(fields[1],fields[2],fields[3],fields[4])] = [fields[5]] + [ float(fields[6]),float(fields[7])]
            else:
                print("ERROR: Only opls type dihedral definitions are currently supported by gen_md_for_vdw.py. Exiting...")
                quit()
                
    # Search for charges based on atom type
    with open(FF_db,'r') as f:
        content=f.readlines()

    Charges = zeros(len(Atom_types))
    for i in range(len(Atom_types)):
        for lines in content:
            fields=lines.split()

            # Skip empty lines
            if len(fields) == 0:
                continue
                    
            if fields[0].lower() in ['charge'] and Atom_types[i] == fields[1]:
                Charges[i] = float(fields[2])

    # Search for polarizability based on atom type
    with open(FF_db,'r') as f:
        content=f.readlines()

    alphas = zeros(len(Atom_types))
    for count_i,i in enumerate(alphas):
        alphas[count_i] = -1
    for i in range(len(Atom_types)):
        for lines in content:
            fields=lines.split()

            # Skip empty lines
            if len(fields) == 0:
                continue
                    
            if fields[0].lower() in ['polar'] and Atom_types[i] == fields[1]:
                alphas[i] = float(fields[2])

    # Search for LP parameters
    LP_params = {}
    with open(FF_db,'r') as f:
        for lines in f:
            fields = lines.split()

            # Skip empty lines
            if len(fields) == 0:
                continue
                    
            if fields[0].lower() in ['charge']:

                tag = split_string(fields[1])
                if tag[0] == 'L' and tag[1] == 'P':
                     LP_params[fields[1]] = float(fields[2])

    # Search for VDW parameters
    VDW_params = {}
    with open(FF_db,'r') as f:
        for lines in f:
            fields = lines.split()

            # Skip empty lines
            if len(fields) == 0:
                continue
                    
            if fields[0].lower() in ['vdw']:

                # Only two parameters are required for lj types
                if fields[3] == "lj":
                    if fields[1] > fields[2]:
                        VDW_params[(fields[1],fields[2])] = [fields[3],float(fields[4]),float(fields[5])]
                    else:
                        VDW_params[(fields[2],fields[1])] = [fields[3],float(fields[4]),float(fields[5])]
                elif fields[3] == "buck":
                    if fields[1] > fields[2]:
                        VDW_params[(fields[1],fields[2])] = [fields[3],float(fields[4]),float(fields[5]),float(fields[6])]
                    else:
                        VDW_params[(fields[2],fields[1])] = [fields[3],float(fields[4]),float(fields[5]),float(fields[6])]

    # Check for missing parameters
    Missing_masses = [ i for i in Atom_types if str(i) not in list(Masses.keys()) ] 
    Missing_charges = [ count_i for count_i,i in enumerate(Charges) if i == -100.0 ]; Missing_charges = [ Atom_types[i] for i in Missing_charges ]
    Missing_bonds = [ i for i in Bond_types if (i[0],i[1]) not in list(Bond_params.keys()) ]
    Missing_angles = [ i for i in Angle_types if (i[0],i[1],i[2]) not in list(Angle_params.keys()) ]
    Missing_dihedrals = [ i for i in Dihedral_types if (i[0],i[1],i[2],i[3],i[4]) not in list(Dihedral_params.keys()) ]
    Missing_impropers = []
    if Improper_flag is True: Missing_impropers = [ i for i in Improper_types if (i[0],i[1],i[2],i[3]) not in list(Improper_params.keys()) ]

    # When the force_read option is set to True, the script will attempt to replace the missing dihedral types with whatever (matching)
    # parameters are in the force-field file. For example, if an opls type is expected but a harmonic type is supplied, the harmonic type
    # will be used IF THE FLAG IS TRUE, otherwise the script will print an error and exit.
    if force_read:

        # Assemble lists of types and keys present in the force field dictionary
        bt,bk = ( [ tuple(i[:2]) for i in list(Bond_params.keys()) ],     [ i for i in list(Bond_params.keys()) ] )
        at,ak = ( [ tuple(i[:3]) for i in list(Angle_params.keys()) ],    [ i for i in list(Angle_params.keys()) ] )
        dt,dk = ( [ tuple(i[:4]) for i in list(Dihedral_params.keys()) ], [ i for i in list(Dihedral_params.keys()) ] )

        # Replace missing dihedral types with whatever matching type is present in the supplied force-field
        for i in set(Missing_dihedrals):
            if i[:4] in dt:
                new_type = dk[dt.index(i[:4])]
                print("using dihedral type {} in place of {}...".format(new_type,i))
                r_inds = [ count_j for count_j,j in enumerate(Dihedral_types) if j == i ]
                for j in r_inds:
                    Dihedral_types[j] = new_type
                Missing_dihedrals = [ j for j in Missing_dihedrals if j != i ]


    # Print diagnostics on missing parameters and quit if the prerequisites are missing.
    if ( len(Missing_masses) + len(Missing_charges) + len(Missing_bonds) + len(Missing_angles) + len(Missing_dihedrals) + len(Missing_impropers) ) > 0:
        print("\nUh Oh! There are missing FF parameters...\n")

        if Missing_masses:
            print("Missing masses for the following atom types: {}".format([ i for i in set(Missing_masses) ]))
        if Missing_charges:
            print("Missing charges for the following atom types: {}".format([ i for i in set(Missing_charges) ]))
        if Missing_bonds:
            print("Missing bond parameters for the following bond types: {}".format([ i for i in set(Missing_bonds) ]))
        if Missing_angles:
            print("Missing angle parameters for the following angle types: {}".format([ i for i in set(Missing_angles) ]))
        if Missing_dihedrals:
            print("Missing dihedral parameters for the following dihedral types: {}".format([ i for i in set(Missing_dihedrals) ]))
        if Improper_flag and Missing_impropers:
            print("Missing improper parameters for the following improper types: {}".format([ i for i in set(Missing_impropers) ]))
        
        print("\nEnsure the specification of the missing parameters. Exiting...")
        quit()

    return list(Bonds),list(Bond_types),Bond_params,list(Angles),list(Angle_types),Angle_params,list(Dihedrals),list(Dihedral_types),Dihedral_params,list(Impropers),list(Improper_types),Improper_params,Charges.tolist(),Masses,VDW_params,alphas.tolist(),LP_params

def split_string(word): 
    return [char for char in word]  

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




# Description: Initialize VDW_dict based on UFF parameters for the initial guess of the fit.
def init_VDW(Data,fit_pairs,read_pairs,Pair_min_dict,fit,FF_db,VDW_dict_prev,UA=0,avoid_read_flag=0,mixing_rule="lb",verbose=True):

    if verbose is True:
        if avoid_read_flag == 0 and mixing_rule in ["none","lb"]:
            print("\nReading VDW parameters or initializing based on UFF values and Lorentz-Berthelot mixing rules\n")
        elif avoid_read_flag == 0 and mixing_rule in ["wh"]:
            print("\nReading VDW parameters or initializing based on UFF values and Waldman-Hagler mixing rules\n")
        elif avoid_read_flag == 1 and mixing_rule in ["none","lb"]:
            print("\nInitializing VDW parameters based on UFF values and Lorentz-Berthelot mixing rules\n")
        elif avoid_read_flag == 1 and mixing_rule in ["wh"]:
            print("\nInitializing VDW parameters based on UFF values and Waldman-Hagler mixing rules\n")

        print("\t{:85s}  {:22s} {:20s}".format("pair_interaction","params","origin"))

    # Initialize UFF parameters (tuple corresponds to eps,sigma pairs for each element)
    # Taken from UFF (Rappe et al. JACS 1992)
    # Note: LJ parameters in the table are specificed in the eps,r_min form rather than eps,sigma
    #       the conversion between r_min and sigma is sigma = r_min/2^(1/6)
    # Note: Units for sigma = angstroms and eps = kcal/mol 
    UFF_dict = { 1:(0.044,2.5711337005530193),  2:(0.056,2.1043027722474816),  3:(0.025,2.183592758161972),  4:(0.085,2.4455169812952313),\
                 5:(0.180,3.6375394661670053),  6:(0.105,3.4308509635584463),  7:(0.069,3.260689308393642),  8:(0.060,3.1181455134911875),\
                 9:(0.050,2.996983287824101),  10:(0.042,2.88918454292912),   11:(0.030,2.657550876212632), 12:(0.111,2.6914050275019648),\
                13:(0.505,4.008153332913386),  14:(0.402,3.82640999441276),   15:(0.305,3.694556984127987), 16:(0.274,3.594776327696269),\
                17:(0.227,3.5163772404999194), 18:(0.185,3.4459962417668324), 19:(0.035,3.396105913550973), 20:(0.238,3.0281647429590133),\
                21:(0.019,2.935511276272418),  22:(0.017,2.828603430095577),  23:(0.016,2.800985569833227), 24:(0.015,2.6931868249382456),\
                25:(0.013,2.6379511044135446), 26:(0.013,2.5942970672246677), 27:(0.014,2.558661118499054), 28:(0.015,2.5248069672097215),\
                29:(0.005,3.113691019900486),  30:(0.124,2.4615531582217574), 31:(0.415,3.904809081609107), 32:(0.379,3.813046513640652),\
                33:(0.309,3.7685015777336357), 34:(0.291,3.746229109780127),  35:(0.251,3.731974730289881), 36:(0.220,3.689211591819145),\
                37:(0.040,3.6651573264293558), 38:(0.235,3.2437622327489755), 39:(0.072,2.980056212179435), 40:(0.069,2.78316759547042),\
                41:(0.059,2.819694442914174),  42:(0.056,2.7190228877643157), 43:(0.048,2.670914356984738), 44:(0.056,2.6397329018498255),\
                45:(0.053,2.6094423454330538), 46:(0.048,2.5827153838888437), 47:(0.036,2.804549164705788), 48:(0.228,2.537279549263686),\
                49:(0.599,3.976080979060334),  50:(0.567,3.9128271700723705), 51:(0.449,3.937772334180300), 52:(0.398,3.982317270087316),\
                53:(0.339,4.009044231631527),  54:(0.332,3.923517954690054),  55:(0.045,4.024189509839913), 56:(0.364,3.2989979532736764),\
                72:(0.072,2.798312873678806),  73:(0.081,2.8241489365048755), 74:(0.067,2.734168165972701), 75:(0.066,2.631714813386562),\
                76:(0.037,2.7796040005978586), 77:(0.073,2.5301523595185635), 78:(0.080,2.453535069758495), 79:(0.039,2.9337294788361374),\
                80:(0.385,2.4098810325696176), 81:(0.680,3.872736727756055),  82:(0.663,3.828191791849038), 83:(0.518,3.893227398273283),\
                84:(0.325,4.195242063722858),  85:(0.284,4.231768911166611),  86:(0.248,4.245132391938716) }

    # Initialize initial_guess dictionary based on the initial_params.db file (if present). 
    # The default is to use these parameters in place of UFF if the file is present.     
    INIT_dict = {}
    if os.path.isfile('initial_params.db') and avoid_read_flag == 0:
        with open('initial_params.db','r') as f:
            for lines in f:
                fields = lines.split()
                if len(fields) > 0 and fields[0] == "vdw":
                    INIT_dict[(fields[1],fields[2])] = [ float(i) for i in fields[4:] ]
                    INIT_dict[(fields[2],fields[1])] = [ float(i) for i in fields[4:] ]

    # Collect types that will be converted into UA
    if UA==1:
        UA_types = {}
        UA_H_types = []
        for i in list(Data.keys()):

            # Iterate over the molecule a adj_mat and identify hydrogens attached to carbons. 
            for count_j,j in enumerate(Data[i]["adj_mat_a"]):
                if int(Data[i]["types_a"][count_j].split('[')[1].split(']')[0]) != 6: continue
                UA_H_types += [ Data[i]["types_a"][count_k] for count_k,k in enumerate(j) if k == 1 and int(Data[i]['types_a'][count_k].split('[')[1].split(']')[0]) == 1 ]

            # Iterate over the molecule b adj_mat and identify hydrogens attached to carbons.
            for count_j,j in enumerate(Data[i]["adj_mat_b"]):
                if int(Data[i]["types_b"][count_j].split('[')[1].split(']')[0]) != 6: continue
                UA_H_types += [ Data[i]["types_b"][count_k] for count_k,k in enumerate(j) if k == 1 and int(Data[i]['types_b'][count_k].split('[')[1].split(']')[0]) == 1 ]

    # Initialize VDW_dict with fit types, first guess based on element types and Lorentz-Berthelot mixing rules
    VDW_dict = {}
    for i in fit_pairs:

        # Grab the base atomic number from the atomtype
        type_1 = int(i[0].split('[')[1].split(']')[0])
        type_2 = int(i[1].split('[')[1].split(']')[0])
        
        # for UA fits, the eps of Hydrogen containing pairs are set to zero
        # so that they do not participate in the fitting. 
        if UA==1 and ( i[0] in UA_H_types or i[1] in UA_H_types ):
            eps    = 0.0
            sigma  = (UFF_dict[type_1][1]+UFF_dict[type_2][1])/2.0
            VDW_dict[i] = (eps,sigma)
            VDW_dict[(i[1],i[0])] = (eps,sigma)

            # Print diagnostic
            if verbose is True:
                print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"Zeroed"))

        # Non-hydrogen UA types read their values from the AA-fit
        elif UA==1:
            eps    = VDW_dict_prev[i][0]
            sigma  = VDW_dict_prev[i][1]
            VDW_dict[i] = (eps,sigma)
            VDW_dict[(i[1],i[0])] = (eps,sigma)

            # Print diagnostic
            if verbose is True:
                print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"AA-Fit"))

        # For all atom fits, initial guesses for eps and sigma are based on scaled UFF parameters with a 
        # check to make sure that the radii aren't too small for the fit set. 
        else:
            if i in INIT_dict:
                eps    = INIT_dict[i][0]
                sigma  = INIT_dict[i][1]
                if sigma < Pair_min_dict[i]*2.**(-1./6.): sigma = Pair_min_dict[i]*2.**(-1./6.)
                VDW_dict[i] = (eps,sigma)
                VDW_dict[(i[1],i[0])] = (eps,sigma)
                if verbose is True:
                    print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"initial_params.db"))
            else:
                if mixing_rule in ["lb",'none']:
                    eps    = (UFF_dict[type_1][0]*UFF_dict[type_2][0])**(0.5) 
                    sigma  = (UFF_dict[type_1][1]+UFF_dict[type_2][1])/2.0
                elif mixing_rule == "wh":
                    sigma  = ((UFF_dict[type_1][1]**(6.0)+UFF_dict[type_2][1]**(6.0))/2.0)**(1.0/6.0)
                    eps    = (UFF_dict[type_1][0]*UFF_dict[type_1][1]**(6.0) * UFF_dict[type_2][0]*UFF_dict[type_2][1]**(6.0) )**(0.5) / sigma**(6.0)
                else: 
                    print("ERROR in init_VDW: {} is not an implemented mixing rule. Exiting...".format(mixing_rule))
                    quit()
                if sigma < Pair_min_dict[i]*2.**(-1./6.): sigma = Pair_min_dict[i]*2.**(-1./6.)
                VDW_dict[i] = (eps,sigma)
                VDW_dict[(i[1],i[0])] = (eps,sigma)
        
                # Print diagnostic
                if verbose is True:
                    print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"UFF"))

    # For diagnostic purposes, the fixed_params.db and user_supplied FF files are reparsed here in order to 
    # print the proper origin information for the user.
    # Read in fixed parameters from the fixed_params.db file
    FIXED_dict = {"vdws":{}}
    if os.path.isfile('fixed_params.db'):         
        FIXED_dict = get_FF_data(['fixed_params.db'])

    # If a FF_db file is supplied then the program avoids fitting any of the
    # pairs that are pairs are in the database
    READ_dict = {"vdws":{}}
    if FF_db != []:        
        READ_dict = get_FF_data(FF_db,mixing_rule=mixing_rule)

    # Initialize read_pairs with entries from the FF_db.
    for i in read_pairs:

        # Grab the base atomic number from the atomtype
        type_1 = int(i[0].split('[')[1].split(']')[0])
        type_2 = int(i[1].split('[')[1].split(']')[0])
        
        # for UA fits, the eps of Hydrogen containing pairs are set to zero
        # so that they do not participate in the fitting. 
        if UA==1 and ( i[0] in UA_H_types or i[1] in UA_H_types ):
            VDW_dict[i] = (0.0,1.0)
            VDW_dict[(i[1],i[0])] = (0.0,1.0)
            continue
        else:        

            # Mixing rule protocols
            if mixing_rule != "none": 

                # Use mixing rule protocol for read parameters
                if (i[0],i[0],'lj') in list(READ_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(READ_dict["vdws"].keys()):

                    eps_1 = float(READ_dict["vdws"][(i[0],i[0],'lj')][3])
                    sigma_1 = float(READ_dict["vdws"][(i[0],i[0],'lj')][4])
                    eps_2 = float(READ_dict["vdws"][(i[1],i[1],'lj')][3])
                    sigma_2 = float(READ_dict["vdws"][(i[1],i[1],'lj')][4])

                # Use mixing rule protocol for fixed parameters
                elif (i[0],i[0],'lj') in list(FIXED_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(FIXED_dict["vdws"].keys()):

                    eps_1 = float(FIXED_dict["vdws"][(i[0],i[0],'lj')][3])
                    sigma_1 = float(FIXED_dict["vdws"][(i[0],i[0],'lj')][4])
                    eps_2 = float(FIXED_dict["vdws"][(i[1],i[1],'lj')][3])
                    sigma_2 = float(FIXED_dict["vdws"][(i[1],i[1],'lj')][4])

                # Use mixing rule protocol for fixed/read parameters
                elif (i[0],i[0],'lj') in list(FIXED_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(READ_dict["vdws"].keys()):

                    eps_1 = float(FIXED_dict["vdws"][(i[0],i[0],'lj')][3])
                    sigma_1 = float(FIXED_dict["vdws"][(i[0],i[0],'lj')][4])
                    eps_2 = float(READ_dict["vdws"][(i[1],i[1],'lj')][3])
                    sigma_2 = float(READ_dict["vdws"][(i[1],i[1],'lj')][4])

                # Use mixing rule protocol for read/fixed parameters
                elif (i[0],i[0],'lj') in list(READ_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(FIXED_dict["vdws"].keys()):

                    eps_1 = float(READ_dict["vdws"][(i[0],i[0],'lj')][3])
                    sigma_1 = float(READ_dict["vdws"][(i[0],i[0],'lj')][4])
                    eps_2 = float(FIXED_dict["vdws"][(i[1],i[1],'lj')][3])
                    sigma_2 = float(FIXED_dict["vdws"][(i[1],i[1],'lj')][4])

                # Apply mixing rule
                if mixing_rule in ["lb"]:
                    eps    = (eps_1*eps_2)**(0.5)
                    sigma  = (sigma_1+sigma_2)/2.0
                elif mixing_rule == "wh":
                    sigma  =  ( ( sigma_1**(6.0) + sigma_2**(6.0) ) / 2.0 )**(1.0/6.0)
                    eps    =  ( eps_1*sigma_1**(6.0) * eps_2*sigma_2**(6.0) )**(0.5) / sigma**(6.0)    
                else: 
                    print("ERROR in init_VDW: {} is not an implemented mixing rule. Exiting...".format(mixing_rule))
                    quit()
                VDW_dict[i] = (eps,sigma)
                VDW_dict[(i[1],i[0])] = (eps,sigma)

                # Print commands
                if verbose is True:
                    if (i[0],i[0],'lj') in list(READ_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(READ_dict["vdws"].keys()):
                        print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"supplied_file"))
                    elif (i[0],i[0],'lj') in list(FIXED_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(FIXED_dict["vdws"].keys()):
                        print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"fixed_params.db"))
                    elif (i[0],i[0],'lj') in list(FIXED_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(READ_dict["vdws"].keys()):
                        print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"fixed/supplied mixing"))
                    elif (i[0],i[0],'lj') in list(READ_dict["vdws"].keys()) and (i[1],i[1],'lj') in list(FIXED_dict["vdws"].keys()):
                        print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"fixed/supplied mixing"))



            # Add the pair type (along with the reverse combination) using the parameters found in the READ_dict
            elif (i[0],i[1],'lj') in list(READ_dict["vdws"].keys()) or (i[1],i[0],'lj') in list(READ_dict["vdws"].keys()):
                if (i[0],i[1],'lj') in list(READ_dict["vdws"].keys()):
                    VDW_dict[i] = (float(READ_dict["vdws"][(i[0],i[1],'lj')][3]),float(READ_dict["vdws"][(i[0],i[1],'lj')][4]))
                    VDW_dict[(i[1],i[0])] = (float(READ_dict["vdws"][(i[0],i[1],'lj')][3]),float(READ_dict["vdws"][(i[0],i[1],'lj')][4]))            

                else:
                    VDW_dict[i] = (float(READ_dict["vdws"][(i[1],i[0],'lj')][3]),float(READ_dict["vdws"][(i[1],i[0],'lj')][4]))
                    VDW_dict[(i[1],i[0])] = (float(READ_dict["vdws"][(i[1],i[0],'lj')][3]),float(READ_dict["vdws"][(i[1],i[0],'lj')][4]))

                if verbose is True:
                    print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"supplied_file"))

            # Add the pair type (along with the reverse combination) using the parameters found in the FIXED_dict
            elif (i[0],i[1],'lj') in list(FIXED_dict["vdws"].keys()) or (i[1],i[0],'lj') in list(FIXED_dict["vdws"].keys()):

                if (i[0],i[1],'lj') in list(FIXED_dict["vdws"].keys()):
                    VDW_dict[i] = (float(FIXED_dict["vdws"][(i[0],i[1],'lj')][3]),float(FIXED_dict["vdws"][(i[0],i[1],'lj')][4]))
                    VDW_dict[(i[1],i[0])] = (float(FIXED_dict["vdws"][(i[0],i[1],'lj')][3]),float(FIXED_dict["vdws"][(i[0],i[1],'lj')][4]))            

                else:
                    VDW_dict[i] = (float(FIXED_dict["vdws"][(i[1],i[0],'lj')][3]),float(FIXED_dict["vdws"][(i[1],i[0],'lj')][4]))
                    VDW_dict[(i[1],i[0])] = (float(FIXED_dict["vdws"][(i[1],i[0],'lj')][3]),float(FIXED_dict["vdws"][(i[1],i[0],'lj')][4]))

                if verbose is True:
                    print("\tpair {:80s}: {:20s} {:20s}".format(i,", ".join([ "{:<10.4f}".format(j) for j in VDW_dict[i]]),"fixed_params.db"))

            else:
                print("ERROR in init_VDW: Information for pair type {} was expected in the the fixed_params.db".format((i[0],i[1],'lj')))
                print("                   and/or the supplied FF files, but it could not be found. Please check that the parameters exist. Exiting...")
                quit()


    # Perform conversions to buckingham style. Conversions follow the approach of Rappe et al. 
    # with a "form factor" of 12.0 which reproduces the long-range LJ behavior
    if fit == 'buck':
        for i in list(VDW_dict.keys()):
            A = VDW_dict[i][0]*0.5*exp(12.0)
            B = VDW_dict[i][1]/12.0
            C = VDW_dict[i][0] * 2.0 * VDW_dict[i][1]**(6.0)
            VDW_dict[i] = (A,B,C)

    return VDW_dict

# Description: Initialize VDW_dict based on UFF parameters for the initial guess of the fit.
def initialize_VDW(atomtypes,sigma_scale=1.0,eps_scale=1.0,VDW_FF={},Force_UFF=0,mixing_rule='lb'):

    # Initialize UFF parameters (tuple corresponds to eps,sigma pairs for each element)
    # Taken from UFF (Rappe et al. JACS 1992)
    # Note: LJ parameters in the table are specificed in the eps,r_min form rather than eps,sigma
    #       the conversion between r_min and sigma is sigma = r_min/2^(1/6)
    # Note: Units for sigma = angstroms and eps = kcal/mol 
    UFF_dict = { 1:(0.044,2.5711337005530193),  2:(0.056,2.1043027722474816),  3:(0.025,2.183592758161972),  4:(0.085,2.4455169812952313),\
                 5:(0.180,3.6375394661670053),  6:(0.105,3.4308509635584463),  7:(0.069,3.260689308393642),  8:(0.060,3.1181455134911875),\
                 9:(0.050,2.996983287824101),  10:(0.042,2.88918454292912),   11:(0.030,2.657550876212632), 12:(0.111,2.6914050275019648),\
                13:(0.505,4.008153332913386),  14:(0.402,3.82640999441276),   15:(0.305,3.694556984127987), 16:(0.274,3.594776327696269),\
                17:(0.227,3.5163772404999194), 18:(0.185,3.4459962417668324), 19:(0.035,3.396105913550973), 20:(0.238,3.0281647429590133),\
                21:(0.019,2.935511276272418),  22:(0.017,2.828603430095577),  23:(0.016,2.800985569833227), 24:(0.015,2.6931868249382456),\
                25:(0.013,2.6379511044135446), 26:(0.013,2.5942970672246677), 27:(0.014,2.558661118499054), 28:(0.015,2.5248069672097215),\
                29:(0.005,3.113691019900486),  30:(0.124,2.4615531582217574), 31:(0.415,3.904809081609107), 32:(0.379,3.813046513640652),\
                33:(0.309,3.7685015777336357), 34:(0.291,3.746229109780127),  35:(0.251,3.731974730289881), 36:(0.220,3.689211591819145),\
                37:(0.040,3.6651573264293558), 38:(0.235,3.2437622327489755), 39:(0.072,2.980056212179435), 40:(0.069,2.78316759547042),\
                41:(0.059,2.819694442914174),  42:(0.056,2.7190228877643157), 43:(0.048,2.670914356984738), 44:(0.056,2.6397329018498255),\
                45:(0.053,2.6094423454330538), 46:(0.048,2.5827153838888437), 47:(0.036,2.804549164705788), 48:(0.228,2.537279549263686),\
                49:(0.599,3.976080979060334),  50:(0.567,3.9128271700723705), 51:(0.449,3.937772334180300), 52:(0.398,3.982317270087316),\
                53:(0.339,4.009044231631527),  54:(0.332,3.923517954690054),  55:(0.045,4.024189509839913), 56:(0.364,3.2989979532736764),\
                72:(0.072,2.798312873678806),  73:(0.081,2.8241489365048755), 74:(0.067,2.734168165972701), 75:(0.066,2.631714813386562),\
                76:(0.037,2.7796040005978586), 77:(0.073,2.5301523595185635), 78:(0.080,2.453535069758495), 79:(0.039,2.9337294788361374),\
                80:(0.385,2.4098810325696176), 81:(0.680,3.872736727756055),  82:(0.663,3.828191791849038), 83:(0.518,3.893227398273283),\
                84:(0.325,4.195242063722858),  85:(0.284,4.231768911166611),  86:(0.248,4.245132391938716) }

    # NEW: Initialize VDW_dict first guess based on element types and Lorentz-Berthelot mixing rules
    # Order of operations: (1) If the parameters are in the supplied FF database then they are used as is
    # (2) If the self-terms are in the supplied FF datebase then they are used to generate the mixed
    # interactions (3) UFF parameters are used. 
    VDW_dict = {}
    origin = {}
    VDW_styles = []
    for count_i,i in enumerate(atomtypes):
        for count_j,j in enumerate(atomtypes):
            if count_i < count_j:
                continue

            # Check for parameters in the database
            if (i,j) in VDW_FF and Force_UFF != 1:

                # Determine appropriate lammps style
                if VDW_FF[(i,j)][0] == "lj":
                    VDW_type = "lj/cut/coul/long"
                elif VDW_FF[(i,j)][0] == "buck":
                    VDW_type = "buck/coul/long"
                else:
                    print("ERROR in initialize_VDW: only lj and buck pair types are supported. Exiting...")
                    quit()

                # Assign style
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type] + VDW_FF[(i,j)][1:]
                else:
                    VDW_dict[(j,i)] = [VDW_type] + VDW_FF[(i,j)][1:]
                origin[(i,j)] = origin[(j,i)] = "read"

            # Check for reverse combination
            elif (j,i) in VDW_FF and Force_UFF != 1:

                # Determine appropriate lammps style
                if VDW_FF[(j,i)][0] == "lj":
                    VDW_type = "lj/cut/coul/long"
                elif VDW_FF[(j,i)][0] == "buck":
                    VDW_type = "buck/coul/long"
                else:
                    print("ERROR in initialize_VDW: only lj and buck pair types are supported. Exiting...")
                    quit()

                # Assign style
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type] + VDW_FF[(j,i)][1:]
                else:
                    VDW_dict[(j,i)] = [VDW_type] + VDW_FF[(j,i)][1:]
                origin[(i,j)] = origin[(j,i)] = "read"

            # Check if the database has the self-terms necessary for applying mixing rules
            elif (i,i) in VDW_FF and (j,j) in VDW_FF and Force_UFF != 1 and mixing_rule == 'lb':

                # Check compatibility with mixing rules
                if VDW_FF[(i,i)][0] != "lj" or VDW_FF[(j,j)][0] != "lj":
                    print("ERROR in initialize_VDW: only lj styles support mixing rules. Exiting...")
                    quit()

                # Apply mixing rules and assign
                VDW_type = "lj/cut/coul/long"
                eps    = (VDW_FF[(i,i)][1]*VDW_FF[(j,j)][1])**(0.5) * eps_scale
                sigma  = (VDW_FF[(i,i)][2]+VDW_FF[(j,j)][2])/2.0 * sigma_scale
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type,eps,sigma]
                else:
                    VDW_dict[(j,i)] = [VDW_type,eps,sigma]
                origin[(i,j)] = origin[(j,i)] = "lb"

            # Check if the database has the self-terms necessary for applying mixing rules
            elif (i,i) in VDW_FF and (j,j) in VDW_FF and Force_UFF != 1 and mixing_rule == 'wh':

                # Check compatibility with mixing rules
                if VDW_FF[(i,i)][0] != "lj" or VDW_FF[(j,j)][0] != "lj":
                    print("ERROR in initialize_VDW: only lj styles support mixing rules. Exiting...")
                    quit()

                # Apply mixing rules and assign
                VDW_type = "lj/cut/coul/long"
                sigma  = ((VDW_FF[(i,i)][2]**(6.0)+VDW_FF[(j,j)][2]**(6.0))/2.0)**(1.0/6.0)
                eps    = (VDW_FF[(i,i)][1]*VDW_FF[(i,i)][2]**(6.0) * VDW_FF[(j,j)][1]*VDW_FF[(j,j)][2]**(6.0) )**(0.5) / sigma**(6.0)
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type,eps,sigma]
                else:
                    VDW_dict[(j,i)] = [VDW_type,eps,sigma]
                origin[(i,j)] = origin[(j,i)] = "wh"

            # Last resort: Use UFF parameters.
            else:
                VDW_type = "lj/cut/coul/long"
                type_1 = int(i.split('[')[1].split(']')[0])
                type_2 = int(j.split('[')[1].split(']')[0])
                eps    = (UFF_dict[type_1][0]*UFF_dict[type_2][0])**(0.5) * eps_scale
                sigma  = (UFF_dict[type_1][1]+UFF_dict[type_2][1])/2.0 * sigma_scale
                if i > j:
                    VDW_dict[(i,j)] = [VDW_type,eps,sigma]
                else:
                    VDW_dict[(j,i)] = [VDW_type,eps,sigma]
                origin[(i,j)] = origin[(j,i)] = "UFF"
                
            # Collect a list of the LAMMPS styles used in the simulation
            VDW_styles += [VDW_type]

    # Print summary
    print("\n{}".format("*"*177))
    print("* {:^173s} *".format("Initializing VDW parameters for the simulation (those with * were read from the FF file(s))"))
    print("*{}*".format("-"*175))
    print("* {:<50s} {:<50s} {:<20s}  {:<18s} {:<18s} {:<8s}    *".format("Type","Type","VDW_type","eps (kcal/mol)","sigma (angstroms)","origin"))
    print("{}".format("*"*177))
    for j in list(VDW_dict.keys()):
        print("  {:<50s} {:<50s} {:<20s} {:< 18.4f} {:< 18.4f}  {:<18s}".format(j[0],j[1],VDW_dict[j][0],VDW_dict[j][1],VDW_dict[j][2],origin[j]))
    print("")


#    return VDW_dict,sorted(set(VDW_styles))
    return VDW_dict

# Description: Initialize VDW_dict based on UFF parameters for the initial guess of the fit.
def initialize_VDW_old(atomtypes,sigma_scale=1.0,eps_scale=1.0,VDW_type="lj/cut/coul/long",VDW_FF={},Force_UFF=0):

    # Initialize UFF parameters (tuple corresponds to eps,sigma pairs for each element)
    # Taken from UFF (Rappe et al. JACS 1992)
    # Note: LJ parameters in the table are specificed in the eps,r_min form rather than eps,sigma
    #       the conversion between r_min and sigma is sigma = r_min/2^(1/6)
    # Note: Units for sigma = angstroms and eps = kcal/mol 
    UFF_dict = { 1:(0.044,2.5711337005530193),  2:(0.056,2.1043027722474816),  3:(0.025,2.183592758161972),  4:(0.085,2.4455169812952313),\
                 5:(0.180,3.6375394661670053),  6:(0.105,3.4308509635584463),  7:(0.069,3.260689308393642),  8:(0.060,3.1181455134911875),\
                 9:(0.050,2.996983287824101),  10:(0.042,2.88918454292912),   11:(0.030,2.657550876212632), 12:(0.111,2.6914050275019648),\
                13:(0.505,4.008153332913386),  14:(0.402,3.82640999441276),   15:(0.305,3.694556984127987), 16:(0.274,3.594776327696269),\
                17:(0.227,3.5163772404999194), 18:(0.185,3.4459962417668324), 19:(0.035,3.396105913550973), 20:(0.238,3.0281647429590133),\
                21:(0.019,2.935511276272418),  22:(0.017,2.828603430095577),  23:(0.016,2.800985569833227), 24:(0.015,2.6931868249382456),\
                25:(0.013,2.6379511044135446), 26:(0.013,2.5942970672246677), 27:(0.014,2.558661118499054), 28:(0.015,2.5248069672097215),\
                29:(0.005,3.113691019900486),  30:(0.124,2.4615531582217574), 31:(0.415,3.904809081609107), 32:(0.379,3.813046513640652),\
                33:(0.309,3.7685015777336357), 34:(0.291,3.746229109780127),  35:(0.251,3.731974730289881), 36:(0.220,3.689211591819145),\
                37:(0.040,3.6651573264293558), 38:(0.235,3.2437622327489755), 39:(0.072,2.980056212179435), 40:(0.069,2.78316759547042),\
                41:(0.059,2.819694442914174),  42:(0.056,2.7190228877643157), 43:(0.048,2.670914356984738), 44:(0.056,2.6397329018498255),\
                45:(0.053,2.6094423454330538), 46:(0.048,2.5827153838888437), 47:(0.036,2.804549164705788), 48:(0.228,2.537279549263686),\
                49:(0.599,3.976080979060334),  50:(0.567,3.9128271700723705), 51:(0.449,3.937772334180300), 52:(0.398,3.982317270087316),\
                53:(0.339,4.009044231631527),  54:(0.332,3.923517954690054),  55:(0.045,4.024189509839913), 56:(0.364,3.2989979532736764),\
                72:(0.072,2.798312873678806),  73:(0.081,2.8241489365048755), 74:(0.067,2.734168165972701), 75:(0.066,2.631714813386562),\
                76:(0.037,2.7796040005978586), 77:(0.073,2.5301523595185635), 78:(0.080,2.453535069758495), 79:(0.039,2.9337294788361374),\
                80:(0.385,2.4098810325696176), 81:(0.680,3.872736727756055),  82:(0.663,3.828191791849038), 83:(0.518,3.893227398273283),\
                84:(0.325,4.195242063722858),  85:(0.284,4.231768911166611),  86:(0.248,4.245132391938716) }

    # NEW: Initialize VDW_dict first guess based on element types and Lorentz-Berthelot mixing rules
    # Order of operations: (1) If the parameters are in the supplied FF database then they are used as is
    # (2) If the self-terms are in the supplied FF datebase then they are used to generate the mixed
    # interactions (3) UFF parameters are used. 
    VDW_dict = {}
    for count_i,i in enumerate(atomtypes):
        for count_j,j in enumerate(atomtypes):
            if count_i < count_j:
                continue

            # Check for parameters in the database
            if (i,j) in VDW_FF and Force_UFF != 1:
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type,VDW_FF[(i,j)][1],VDW_FF[(i,j)][2]]
                else:
                    VDW_dict[(j,i)] = [VDW_type,VDW_FF[(i,j)][1],VDW_FF[(i,j)][2]]
            elif (j,i) in VDW_FF and Force_UFF != 1:
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type,VDW_FF[(j,i)][1],VDW_FF[(j,i)][2]]
                else:
                    VDW_dict[(j,i)] = [VDW_type,VDW_FF[(j,i)][1],VDW_FF[(j,i)][2]]

            # Check if the database has the self-terms necessary for applying mixing rules
            elif (i,i) in VDW_FF and (j,j) in VDW_FF and Force_UFF != 1:
                eps    = (VDW_FF[(i,i)][1]*VDW_FF[(j,j)][1])**(0.5) * eps_scale
                sigma  = (VDW_FF[(i,i)][2]+VDW_FF[(j,j)][2])/2.0 * sigma_scale
                if i > j: 
                    VDW_dict[(i,j)] = [VDW_type,eps,sigma]
                else:
                    VDW_dict[(j,i)] = [VDW_type,eps,sigma]

            # Last resort: Use UFF parameters.
            else:
                type_1 = int(i.split('[')[1].split(']')[0])
                type_2 = int(j.split('[')[1].split(']')[0])
                eps    = (UFF_dict[type_1][0]*UFF_dict[type_2][0])**(0.5) * eps_scale
                sigma  = (UFF_dict[type_1][1]+UFF_dict[type_2][1])/2.0 * sigma_scale
                if i > j:
                    VDW_dict[(i,j)] = [VDW_type,eps,sigma]
                else:
                    VDW_dict[(j,i)] = [VDW_type,eps,sigma]

    # Print summary
    print("\n{}".format("*"*167))
    print("* {:^163s} *".format("Initializing VDW parameters for the simulation (those with * were read from the FF file(s))"))
    print("*{}*".format("-"*165))
    print("* {:<50s} {:<50s} {:<20s}  {:<18s} {:<18s}   *".format("Type","Type","VDW_type","eps (kcal/mol)","sigma (angstroms)"))
    print("{}".format("*"*167))
    for j in list(VDW_dict.keys()):
        if j in VDW_FF and Force_UFF != 1:
            print("  {:<50s} {:<50s} {:<20s} {:< 18.4f} {:< 18.4f} *".format(j[0],j[1],VDW_dict[j][0],VDW_dict[j][1],VDW_dict[j][2]))
        else:
            print("  {:<50s} {:<50s} {:<20s} {:< 18.4f} {:< 18.4f}".format(j[0],j[1],VDW_dict[j][0],VDW_dict[j][1],VDW_dict[j][2]))
    print("")


    return VDW_dict

# Description:   A wrapper for the commands to create a data dictionary holding the FF, mode, and geometry information for each molecule
# Inputs:  Filename:
# Returns: Data: a dictionary keyed to the list of coord_files
#                a keyed entry is itself a dictionary holding the geometry, 
#                elements, charges, bonds, bond_types, angles, angle_types,
#                dihedrals, dihedral_types and masses.
#                subdictionary keys (e.g., Data[*]["Geometry"])
#                Geometry: array of the molecule's geometry
#                Elements: list of the molecule's elements (indexed to the rows of "Geometry")
#                Atom_types: list of the molecule's atom_types (indexed to the rows of "Geometry")
#                Charges: list of the molecule's charges (indexed to the rows of "Geometry")
#                Bonds: a list of the molecule's bonds, each element being a tuple of the atom indices that are bonded
#                Angles: a list of the molecule's angles, each element being a tuple of the atom indices that form an angle mode
#                Dihedrals: a list of the molecule's dihedrals, each element being a tuple of the atom indices that form a dihedral
#                Bond_types: a list of the molecule's bond_types, indexed to Bonds
#                Angle_types: a list of the molecule's angle_types, indexed to Angles
#                Dihedral_types: a list of the molecule's dihedral_types, indexed to Dihedrals
#                Impropers: a list of the molecule's improper dihedrals by atom tindex
#                Improper_types: a list of the molecule's improper_types, indexed to Impropers
#
def get_data(FF_all,coord_files,q_list,gens,Improper_flag=False,force_read_opt=False,remove_multi=False):

    # Initialize dictionary to hold the FF, mode, and geometry information about each molecule
    Data = {}
    N_mol = [1]

    # Iterate over all molecules being simulated and collect their FF, mode, and geometry information.
    for count_i,i in enumerate(coord_files):

        # Initialize dictionary for this geometry and set number of molecules to place in the simulated system
        Data[i] = {}      
        Data[i]["N_mol"] = N_mol[count_i]

        # Extract Element list and Coord list from the file
        Data[i]["Elements"],Data[i]["Geometry"] = xyz_parse(i)

        # Generate adjacency table
        Data[i]["Adj_mat"] = Table_generator(Data[i]["Elements"],Data[i]["Geometry"])

        # Automatically parse the atom types based on the adjacency matrix and gens argument. 
        # Find Hybridizations (needed to identify EZ steroisomers in id_types call)
        Hybridizations = Hybridization_finder(Data[i]["Elements"],Data[i]["Adj_mat"],force_opt=True)
    
        # Find atom types
        print("Determining atom types based on a {}-bond deep search...".format(i,gens))
        Data[i]["Atom_types"] = id_types(Data[i]["Elements"],Data[i]["Adj_mat"],gens)

        # Determine bonding matrix for the compound
        if q_list[count_i] == "none":
            q_tmp = 0.0
        else:
            q_tmp = q_list[count_i]
        Data[i]["Bond_mats"] = find_lewis(Data[i]["Elements"], Data[i]["Adj_mat"], q_tot=int(q_tmp),b_mat_only=True,verbose=False)        

        # Check the number of molecules
        mol_in_in = mol_count(Data[i]["Adj_mat"])
        if mol_in_in > 1:
            print("ERROR: {} molecules were discovered in geometry {}. Check the geometry of the input file. Exiting...".format(mol_in_in,i))
            quit()

        # Generate list of bonds angles and dihedrals    
        print("\n{}\n* {:^163s} *\n{}\n".format("*"*167,"Parsing Modes and FF Information for Molecule {}".format(i),"*"*167))
        Data[i]["Bonds"],Data[i]["Bond_types"],Data[i]["Bond_params"],Data[i]["Angles"],Data[i]["Angle_types"],Data[i]["Angle_params"],\
        Data[i]["Dihedrals"],Data[i]["Dihedral_types"],Data[i]["Dihedral_params"],Data[i]["Impropers"],Data[i]["Improper_types"],Data[i]["Improper_params"],Data[i]["Charges"],Data[i]["Masses"],Data[i]["VDW_params"],Data[i]['alphas'],Data[i]['LP_params'] =\
            Find_parameters(Data[i]["Adj_mat"],Data[i]["Bond_mats"],Data[i]["Geometry"],Data[i]["Atom_types"],FF_db=FF_all,Improper_flag = Improper_flag, force_read=force_read_opt,remove_multi=remove_multi)

        # Print System characteristics
        print("\n{}\n* {:^163s} *\n{}\n".format("*"*167,"Mode Summary for Molecule {}".format(i),"*"*167))
        print("\nAtom_types ({}):\n".format(len(set(Data[i]["Atom_types"]))))
        for j in sorted(set(Data[i]["Atom_types"])):
            print("\t{}".format(j))
        print("\nBond types ({}):\n".format(len(set(Data[i]["Bond_types"]))))
        for j in sorted(set(Data[i]["Bond_types"])):
            print("\t{}".format(j))
        print("\nAngle types ({}):\n".format(len(set(Data[i]["Angle_types"]))))
        for j in sorted(set(Data[i]["Angle_types"])):
            print("\t{}".format(j))
        print("\nDihedral types ({}):\n".format(len(set(Data[i]["Dihedral_types"]))))
        for j in sorted(set(Data[i]["Dihedral_types"])):
            print("\t{}".format(j))
        if Improper_flag:
            print("\nImproper types ({}):\n".format(len(set(Data[i]["Improper_types"]))))
            for j in sorted(set(Data[i]["Improper_types"])):
                print("\t{}".format(j))
        
        # Subtract off residual
        print("\n{:40s} {}".format("Residual Charge:",sum(Data[i]["Charges"])))
        if q_list[count_i] == "round":
            q_tot = int(round(sum(Data[i]["Charges"])))
        elif q_list[count_i] != "none":
            q_tot = int(q_list[count_i])

        # Avoid subtracting off the residual if q_list is none
        if q_list[count_i] != "none":
            correction = (float(q_tot)-sum(Data[i]["Charges"]))/float(len(Data[i]["Atom_types"]))
            print("{:40s} {}".format("Charge added to each atom:",correction))
            for j in range(len(Data[i]["Atom_types"])):
                Data[i]["Charges"][j] += correction

        print("{:40s} {}".format("Final Total Charge:",sum(Data[i]["Charges"])))
        print("\n{}".format("*"*192))
        print("* {:^188s} *".format("System Characteristics"))
        print("*{}*".format("-"*190))
        print("* {:<87s} {:<24s} {:<24s} {:<25s} {:<25s}*".format("Type","Element","Mass","Charge","Polarizability"))
        print("{}".format("*"*192))
        for j in range(len(Data[i]["Atom_types"])):
            print(" {:<88s} {:<23s} {:< 24.6f} {:< 24.6f} {:< 24.6f}".format(Data[i]["Atom_types"][j],Data[i]["Elements"][j],Data[i]["Masses"][str(Data[i]["Atom_types"][j])],Data[i]["Charges"][j],Data[i]["alphas"][j]))

    return Data

# A simple wrapper function for generating a fix statement for scaled charges
def scale_charges(charge_scale,Atom_type_dict,Atom_types,Charges):

    fixes = []
    if charge_scale != 1.0:
        fixes = '# Setting up fixes to scale the charges by {}\n'.format(charge_scale)
        for i in sorted([ Atom_type_dict[j] for j in list(Atom_type_dict.keys()) ]):
            fixes += 'group charge_{} type {}\n'.format(i,i)        
        for i in sorted([ Atom_type_dict[j] for j in list(Atom_type_dict.keys()) ]):
            index = Atom_types.index([ k for k in list(Atom_type_dict.keys()) if Atom_type_dict[k] == i ][0])
            fixes += 'variable charge_{} equal {}\n'.format(i,Charges[index]*charge_scale)
        for i in sorted([ Atom_type_dict[j] for j in list(Atom_type_dict.keys()) ]):        
            index = Atom_types.index([ k for k in list(Atom_type_dict.keys()) if Atom_type_dict[k] == i ][0])
            fixes += 'fix {} charge_{} adapt 0 atom charge v_charge_{} reset yes\n'.format(i,i,i)
        fixes += '\n'
    return fixes

                    

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder+"/gen_md_for_CHARMM.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass

    def flush(self):
        pass

def miss2default(config,other_opt):

 
    # Create default dict
    options = ['coord_files','FF_files','N_mol','T_equil','t_equil','T_anneal','t_anneal','T_dilatometry','N_dilatometry',\
               't_ext','pressure','density','N_density','fixes','frequency','q_list','outputname','timestep','gens','onefourscale_coul'\
               ,'onefourscale_lj','eps_scale','sigma_scale','charge_scale','force_UFF','mixing_rule','pair_styles','tail_opt','improper_flag'\
               ,'velocity_flag','molecule_flag','remove_multi','force_read_opt']

    defaults = ["", "25", 400, 1E6, 400, 1E6, None, 20, 1E6, "1", 0, 0, "", 1000, "none", '', 1.0, 2, 0.0, 0.0, 1.0, 1.0, 1.0, 0, "none",\
                'lj/cut/coul/long 14.0 14.0', False, False, False, False, False, False]  

    N_position = int (len(options) - len(defaults))

    default = {}
    for count_i,i in enumerate(defaults):
      default[options[N_position + count_i]] = i


    # Combine config and other_opt
    option = {}
    for key in config:
      option[key] = config[key]

    #other_opt's priority > config
    for key in other_opt:
      option[key] = other_opt[key]

    missing = [ i in option for i in options]
    
    
    # set missing option to default
    for count_i,i in enumerate(missing):
      if i is False:
         option[options[count_i]] = default[options[count_i]]

    return option

if __name__ == "__main__":
   main(sys.argv[1:])
