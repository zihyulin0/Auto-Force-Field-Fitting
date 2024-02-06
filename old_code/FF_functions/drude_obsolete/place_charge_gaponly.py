#!/bin/env python                                                                                                                                                             
# Author: Lin (lin1209@purdue.edu)

import sys,argparse,subprocess,os,time,math
from subprocess import PIPE
#from numpy import *
from math import sqrt,sin,cos,tan,factorial
import math
#from scipy import *
#from numpy.linalg import *
import shutil
from shutil import move,copyfile
from fnmatch import fnmatch
from copy import deepcopy
from itertools import permutations

# Add TAFFY Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from kekule import *
from transify import *
from adjacency import *
from file_parsers import *
from id_types import *
import json, codecs
#from __builtin__ import any  # any has a namespace conflict with numpy. XXX WE REALLY NEED TO CLEAN UP THESE NAMESPACE ISSUES
from builtins import any
import numpy as np

# Lib for Grid_gen
sys.path.append('/home/lin1209/Grid_gen/bin/')
import jenerate_resp_points as GRC # Generate Resp Charge 
from scipy.spatial.distance import cdist

def main(argv):

    parser = argparse.ArgumentParser(description='place perturbation charges for input xyz file according to Roux\'s group method')

    #required (positional) arguments                                                                                                  
    parser.add_argument('coord_files', help = 'Space delimited string holding the names of the .xyz input geometries of one or several molecules in need of parametrization. Wildcards are supported. ')

    parser.add_argument('-d', dest='density', default=6,
                        help = 'Sets density for generating the grid points (default: 6)')

    parser.add_argument('-radii_factor', dest='size_factor', default='10',
                        help = 'radii factor for VDW radii default: 10 (10 is VDW_radii*1.0)')

    parser.add_argument('-d_factor', dest='density_scale_factor', default=1,
                        help = 'density scale factor (default: 1)')
    
    parser.add_argument('-o', dest='file_name', default='charge',
                        help = 'output xyz name (default: charge.xyz)')

    parser.add_argument('-ofolder', dest='folder_name', default='pointcharge',
                        help = 'output folder name (default: pointcharge)')

    parser.add_argument('-c', dest='criteria', default=5.0,
                        help = 'criteria for overlap intersection vector (default: 5 degree)')
    parser.add_argument('-extra', dest='extra', default=0.1,
                        help = 'extra density factor for placing charges in gap (default: 0.1)')

    parser.add_argument('--test', dest='test', default=False, action='store_const', const=True,
                        help = 'when this falg is on, test mode is True, will report grid pts generated only, this is for getting extra_scale for placing gap charges')


    # Make relevant inputs lowercase
    args=parser.parse_args(argv)    

    # Create fake option for function input
    option = dict(args.__dict__.items())

    # Find files
    args.coord_files = args.coord_files.split()
    wild_card_files  = [ i for i in args.coord_files if "*" in i ]
    args.coord_files = [ i for i in args.coord_files if "*" not in i ]
    for i in wild_card_files:
        path = '/'.join(i.split('/')[:-1])
        if len(path) == 0:
          path = "."
        args.coord_files += [ path+"/"+files for files in os.listdir(path) if fnmatch(files,i) ] 

    # Handle "./" condition
    for count_i,i in enumerate(args.coord_files):
        if i[:2] == "./": args.coord_files[count_i] = args.coord_files[count_i][2:]    
    args.coord_files = list(set(args.coord_files))

    # Perform some consistency checks
    if False in [ i.split('.')[-1] == 'xyz' for i in args.coord_files ]:
        print( "ERROR: Check to ensure that the input file(s) are in .xyz format.")
        return
    if False in [ os.path.isfile(i) for i in args.coord_files ]:
        print( "ERROR: Could not find file {}. Exiting...".format(next([ i for i in args.coord_files if os.path.isfile(i) == False ])))


    # Create geo_dict for function input
    geo_dict_in = {}
    for file in args.coord_files:
    
      elements,geo = xyz_parse(file)
      inchikey = GetInchi(elements,geo)
      geo_dict_in[inchikey] = {}
      geo_dict_in[inchikey]["elements"] = elements
      geo_dict_in[inchikey]["geo"] = geo
    
    total_pc = fun(option,geo_dict_in)
    
    return total_pc



def fun(option,geo_dict):
    gens =2 
    q_tot = 0
    arg_density = float(option['density'])
    size_factor = float(option['size_factor'])
    density_scale_factor = float(option['density_scale_factor']) 
    file_name = option['file_name']
    folder_name = option['folder_name']
    criteria = float(option['criteria'])
    extra = float(option['extra'])
    test = option['test']
    
    # Initialize dictionary and lists for holding modes and geometry information
    #geo_dict    = {}
    bonds       = []
    bond_types  = []
    bond_files  = []
    angles      = []
    angle_types = []
    angle_files = []
    dihedrals   = []
    dihedral_types = []
    dihedral_files = []

    for file in geo_dict:        
        inchikey_in = file

        # Print the name
        print( "  {}".format(file))


        # Generate adjacency table
        geo_dict[file]["adj_mat"] = Table_generator(geo_dict[file]["elements"],geo_dict[file]["geo"])
        

        # Find Hybridizations (needed to identify EZ steroisomers in id_types call)
        hybridizations = Hybridization_finder(geo_dict[file]["elements"],geo_dict[file]["adj_mat"])

        # Find atom types
        geo_dict[file]["atom_types"] = id_types(geo_dict[file]["elements"],geo_dict[file]["adj_mat"],gens)


        # Determine bonding matrix for the compound
        geo_dict[file]["bond_mat"] = find_lewis(geo_dict[file]["elements"], geo_dict[file]["adj_mat"], q_tot=q_tot, b_mat_only=True,verbose=False)

        # Find modes (returns all modes in canonicalized format)
        bonds_tmp,angles_tmp,dihedrals_tmp,one_fives_tmp,bond_types_tmp,angle_types_tmp,dihedral_types_tmp,one_five_types_tmp = Find_modes(geo_dict[file]["adj_mat"],geo_dict[file]["atom_types"],geo_dict[file]["bond_mat"],return_all=1)

        if(test):
            test_pt = grid_gen(geo_dict[file]["elements"],geo_dict[file]["geo"],6,size_factor,extra)
            return

        grid_points = grid_gen(geo_dict[file]["elements"],geo_dict[file]["geo"],arg_density,size_factor,density_scale_factor)
        # additional density scale factor, this is by trial and error to get an approximate grid pts that are similar to the final # of charges placed reported by Roux's group
        total_charges=place_all(criteria,grid_points,geo_dict[file]["geo"],bonds_tmp,hybridizations,geo_dict[file]["elements"],arg_density,size_factor,extra,file_name,folder_name)


    return total_charges



def remove_contact(grid_points,intersec):
    shape1 = np.array(grid_points).shape
    shape2 = np.array(intersec).shape
    if (len(shape1) ==1 and shape1[0] == 0) or (len(shape2) == 1 and shape2[0] == 0):
      return grid_points
    if len(shape1) == 1 and shape1[0] != 0:
      grid_points = np.array(grid_points).reshape(-1,1)
    if len(shape2) == 1 and shape2[0] != 0:
      intersec = np.array(intersec).reshape(-1,1)
    Dist_Mat = np.triu(cdist(grid_points,intersec))
    rm_ind,y_ind = np.where( (Dist_Mat < 1.5) & (Dist_Mat > 0.0))
    rm_ind = list(dict.fromkeys(rm_ind))
    
    print("{} charges removed due to close contact w/ already placed charge".format(len(rm_ind)))
    grid_points = [ i for count_i,i in enumerate(grid_points) if count_i not in rm_ind]


    return grid_points

def place_all(criteria,points,geo,bonds,hybridizations,element,density,size_factor,extra_d_scale,name,folder_name):

    if(os.path.isdir(folder_name)):
         print("ERROR: output folder {}  exist, exiting avoid overwriting....".format(folder_name))
         quit()
    os.makedirs(folder_name) 
    os.chdir(folder_name)
    # place charges in the gap
    intersec3 = place_gap(element,geo,size_factor,extra_d_scale)

    total_charges = len(intersec3)
      
    print("Total {} charges placed!".format(total_charges))

    name = name + '.xyz'
    with open(name,'w') as f:
         f.write("{} \n\n".format(len(geo)+total_charges))
         for count_i,i in enumerate(geo):
            f.write("{} {} {} {}\n".format(element[count_i],i[0],i[1],i[2]))
         for i in intersec3:
            f.write("Br {} {} {} \n".format(i[0],i[1],i[2]))
 
    count = 0 
    for i in intersec3:
         filename = 'pointcharges.{}.pc'.format(count)
         with open(filename,'w') as f:
            f.write("1\n")
            f.write("0.5 {} {} {} \n".format(i[0],i[1],i[2]))
         count += 1
    os.chdir("..")

    return total_charges
    
            
   #with open('pointcharges.pc','w') as f:
   #      f.write("{} \n".format(total_charges))
   #      for i in intersec1:
   #         f.write("0.5 {} {} {} \n".format(i[0],i[1],i[2]))
   #      for i in intersec2:
   #         f.write("0.5 {} {} {} \n".format(i[0],i[1],i[2]))
   #      for i in intersec3:
   #         f.write("0.5 {} {} {} \n".format(i[0],i[1],i[2]))
   
    


def place_gap(element,geo,size_factor,extra_d_scale):
   
    density = 6 # this doesn't matter cuz extra_d_scale controls it
    grid_points = grid_gen(element,geo,density,size_factor,extra_d_scale)
    print("Generting initial gap grids... {} pts".format(len(grid_points)))

    #with open('gap.xyz','w') as f:
    #     f.write("{} \n\n".format(len(geo)+len(grid_points)))
    #     for count_i,i in enumerate(geo):
    #        f.write("{} {} {} {}\n".format(elements[count_i],i[0],i[1],i[2]))
    #     for i in grid_points:
    #        f.write("F {} {} {} \n".format(i[0],i[1],i[2]))


    print("Placing {} charges in gap...".format(len(grid_points)))
    return grid_points
    

def place_bonds(criteria,points,geo,bonds,element):


    #name = name + '.xyz'
    # Get bonds vectors
    intersec = []
    for pair in bonds:
         # Get bonds vector and bond central coordinate
         v = np.zeros(3)
         bc = np.zeros(3)
         for i in range(0,3):
            v[i] = geo[pair[0]][i] - geo[pair[1]][i] 
            bc[i] = (geo[pair[0]][i] + geo[pair[1]][i]) / 2.0 
         # Get the intersection of bondvector and connolly surface
         intersec_tmp = {}
         for p in points: 
            v2 = np.zeros(3)
            for i in range(0,3):
               v2[i] = p[i] - bc [i]
            angle = angle_between(v,v2)
            if angle < criteria :
               intersec_tmp[angle] = p
         # if more than 2 grid points fits the criteria for one bond, delete the one that's off most
         if len(intersec_tmp) > 2:
            del intersec_tmp[max(intersec_tmp.keys())] 
         elif len(intersec_tmp) == 1:
            print("WARNING: Only 1 charge placed for the bond, try increasing density")
         elif len(intersec_tmp) == 0:
            print("WARNING: 0 charge placed for the bond, try increasing density")
         for key in intersec_tmp:
            intersec.append(intersec_tmp[key])
    # if distance between two pts <1.5 A then remove one
    shape = np.array(intersec).shape
    if len(shape) == 1 and shape[0] == 0 :
      return intersec 
    elif len(shape) == 1 and shape[0] != 0:
      intersec = np.array(intersec).reshape(-1,1)
    Dist_Mat = np.triu(cdist(intersec,intersec))
    rm_ind,y_ind = np.where( (Dist_Mat < 1.5) & (Dist_Mat > 0.0))
    # remove pts that are repeatedly chosen ( because it sneak through Dist_Mat > 0.0)
    tmp_x,tmp_y = np.where( (Dist_Mat == 0.0) )
    for i in range(0,len(tmp_x)):
      if tmp_y[i] > tmp_x[i]:
         rm_ind = np.append(rm_ind,tmp_y[i])
         
    rm_ind = list(dict.fromkeys(rm_ind))
    intersec = [ i for count_i,i in enumerate(intersec) if count_i not in rm_ind] 
      
             
      
    #with open(name,'w') as f:
    #     f.write("{} \n\n".format(len(geo)+len(intersec)))
    #     for count_i,i in enumerate(geo):
    #        f.write("{} {} {} {}\n".format(element[count_i],i[0],i[1],i[2]))
    #     for i in intersec:
    #        f.write("F {} {} {} \n".format(i[0],i[1],i[2]))
    print("Placing {} charges along {} bonds....".format(len(intersec),len(bonds)))  


    return intersec
    
def place_lone_pair(criteria,points,geo,bonds,hybridizations,element):


    #name = name + '.xyz'
    # Get bonds vectors
    intersec = []

    # atoms that have special contraints (sp3 O and sp2 N)
    special = {}
    for count_i,hybrid in enumerate(hybridizations): 
         if hybrid == 'sp3' and element[count_i] == 'O':
            link = []
            for pair in bonds:
               if pair[0] == count_i:
                  link.append(pair[1])
               elif pair[1] == count_i:
                  link.append(pair[0])
            if len(link) == 2:
               special[count_i] = {}
               special[count_i]['element'] = 'N'
               special[count_i]['atom_num'] = count_i 
               special[count_i]['link'] = link
         elif hybrid == 'sp2' and element[count_i] == 'N':
            link = []
            for pair in bonds:
               if pair[0] == count_i:
                  link.append(pair[1])
               elif pair[1] == count_i:
                  link.append(pair[0])
            if len(link) == 2:
               special[count_i] = {}
               special[count_i]['element'] = 'N'
               special[count_i]['atom_num'] = count_i 
               special[count_i]['link'] = link

    for atom in special:
         v1 = np.zeros(3)
         v2 = np.zeros(3)
         for i in range(0,3):
            v1[i] = geo[special[atom]['link'][0]][i] - geo[atom][i] 
            v2[i] = geo[special[atom]['link'][1]][i] - geo[atom][i] 

         # A-P directions for P1, P2', P2'', P3', P3'' (see paper)
         P_d = []
         # points (for visual check if the direction is correct)
         P_d_p = []
         # A-P1 's direction (see paper)
         P_d.append(-bisector(v1,v2))
         # A-P2' 's direction
         P_d.append(np.cross(v1,v2))
         # A-P2'' 's directiono
         P_d.append(np.cross(v2,v1))
         # A-P3' 's direction
         P_d.append(bisector(P_d[1],P_d[0]))
         # A-P3 '' 's direction
         P_d.append(bisector(P_d[2],P_d[0]))
         for i in range(0,5):
            P_d_p.append(geo[atom]+P_d[i])
         intersec_tmp = {}
         count = [0] * 5
         for p in points: 
            v3 = p - geo[atom]
            for i in range(0,5):
               angle = angle_between_normal(P_d[i],v3)
               if angle < criteria :
                  intersec_tmp[angle] = p
                  count[i] += int(1)

         intersec = []
         # if more than 1 grid points fits the criteria for line pair, only keep one that has min angle
         for i in range(0,5):
            begin = int(np.sum(count[:i]))
            end = int(np.sum(count[:i+1])) 
            if count[i] >1:
               key = [min(list(intersec_tmp.keys())[begin:end])]

            elif count[i] == 0:
               print("WARNING: 0 charge placed for one of the direction, try increasing density")
               key = list(intersec_tmp.keys())[begin:end]
            else:
               key = list(intersec_tmp.keys())[begin:end]

            if key != []:
               intersec.append(intersec_tmp[key[0]])
             
    #with open(name,'a') as f:
    #     #f.write("{} \n\n".format(len(geo)+len(intersec)))
    #     #for count_i,i in enumerate(geo):
    #     #   f.write("{} {} {} {}\n".format(element[count_i],i[0],i[1],i[2]))
    #     for i in intersec:
    #        f.write("Cl {} {} {} \n".format(i[0],i[1],i[2]))
    print("Placing {} charges along {} special atoms....".format(len(intersec),len(special.keys())))  

    #with open(name,'w') as f:
    #     f.write("{} \n\n".format(8))
    #     for atom in special:
    #        f.write("{} {} {} {}\n".format(special[atom]['element'],geo[atom][0],geo[atom][1],geo[atom][2]))
    #        for i in special[atom]['link']:
    #           f.write("{} {} {} {}\n".format('O',geo[i][0],geo[i][1],geo[i][2]))
    #        
    #     
    #     for i in range(0,5):  
    #        f.write("{} {} {} {}\n".format('C',P_d_p[i][0],P_d_p[i][1],P_d_p[i][2]))


    return intersec

    
    
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
# Return acute angle
    """ Returns the angle in degrees between vectors 'v1' and 'v2'(will turn >90 to <90)::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))/math.pi*180.0
    if angle>90.0:
      angle = 180.0 - angle
    return angle

def angle_between_normal(v1,v2):
# Return acute and obstuse angle
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))/math.pi*180.0
    
    return angle

def bisector(v1,v2):

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    bi = v1_u + v2_u

    return bi

def grid_gen(elements,geo,arg_density,size_factor,density_scale):
 
    conf = [] 
    for count_i,i in enumerate(elements):
        conf.append(dict(element=i,xyz=tuple(geo[count_i]))) 
    
    # Use Atom class in jenerate_resp.py 
    atoms = [GRC.Atom(**a) for a in conf]
    x, y, z = np.array([c['xyz'] for c in conf]).T
    #elements = [c['element'] for c in conf]
    # array of points for each scale of molecular surface
    #try different density, 6 is a good number
    points = []
    
    # Size factors (multiplication constants) multiplied by the VDW radii of corresponding atoms
    # sizefactor and density scale are set to be the same as 2005 Anisimov et al.
    # doi:10.1021/ct049930p
    for atom in atoms:
        points.extend(atom.radius_points(density=float(arg_density)*density_scale, scale=size_factor/10, check_against=atoms))

    points = np.array(points) 
    rows = points.shape[0]
    print("Total of {} grid pts generated for Radii scale at {}".format(rows,float(size_factor/10)))

    
    #with open('grid.xyz','w') as f:
    #     f.write("{} \n\n".format(len(geo)+rows))
    #     for count_i,i in enumerate(geo):
    #        f.write("{} {} {} {}\n".format(elements[count_i],i[0],i[1],i[2]))
    #     for i in points:
    #        f.write("F {} {} {} \n".format(i[0],i[1],i[2]))
    #np.savetxt('ttt', points, header=str(len(points)), comments='')
    return points


def write_xyz(Output,Elements,Geometry,Atom_types):
    
    # Open file for writing and write header
    fid = open(Output+'.xyz','w')
    fid.write('{}\n\n'.format(len(Elements)))

    for count_i,i in enumerate(Elements):
        fid.write('{: <4} {:< 12.6f} {:< 12.6f} {:< 12.6f} {}\n'.format(i,Geometry[count_i,0],Geometry[count_i,1],Geometry[count_i,2],Atom_types[count_i]))

    fid.close()

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

def write_pdb(Output,Elements,Geometry,Atom_types):

    # Open file for writing and write header
    f = open(Output+'.pdb','w')
    f.write('COMPND    SYSTEM\nAUTHOR    Generated by polygen.py\n')
    count=0
    for count_i,i in enumerate(Elements):
        count=count+1
        f.write('HETATM {:>4} {:<3}  {:>3} 1       {:> 8.3f}{:> 8.3f}{:> 8.3f}  1.00 {: >6}           {:<3}\n'\
                .format(count,i," ",float(Geometry[count_i,0]),float(Geometry[count_i,1]),float(Geometry[count_i,2]),count_i,Atom_types[count_i]))     
    f.close()


# Shortcut for normalizing a vector
def normalize(x):
    return x/sum(x**(2.0))**(0.5)

# Logger object redirects standard output to a file.
class Logger(object):
    def __init__(self,folder):
        self.terminal = sys.stdout
        self.log = open(folder+"/place_charge.log", "a",buffering=1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  
    def flush(self):
        pass
   
if __name__ == "__main__":
   main(sys.argv[1:])



