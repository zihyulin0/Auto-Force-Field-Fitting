#!/bin/env python                                                                                                                                                             
# Author: Zih-Yu Lin (lin1209@purdue.edu)

import sys,argparse,subprocess,os,time,math
# Add TAFFY Lib to path
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from adjacency import Table_generator, find_lewis, Find_modes
from id_types import Hybridization_finder, id_types

sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
from kekule import *
from transify import *
from adjacency import *
from file_parsers import *
from id_types import *

# Lib for Grid_gen
sys.path.append('/home/lin1209/Grid_gen/bin/')
import jenerate_resp_points as GRC # Generate Resp Charge 
from scipy.spatial.distance import cdist


#### output map.xyz : in Bohr (for ORCA)
#### output connoly_vis.xyz (in Anstrom)


def main(argv):

    parser = argparse.ArgumentParser(description='Reads in xyz files and generate connolly grid for ORCA to reevaluate and a xyz file for visualization')

    #required (positional) arguments                                                                                                  
    parser.add_argument('coord_files', help = 'Space delimited string holding the names of the .xyz input geometries of one or several molecules in need of parametrization. Wildcards are supported. ')

    parser.add_argument('-d', dest='density',type=float, default=6,
                        help = 'Sets density for generating the grid points (default: 6)')

    parser.add_argument('-radii_factor', dest='size_factor',type=float, default=10,
                        help = 'radii factor for VDW radii, this is only used for test default: 10 (10 is VDW_radii*1.0)')

    parser.add_argument('-gen', dest='gens',type=int, default=2,
                        help = 'taffi atom-type generation, default: 2')

    parser.add_argument('-qtot', dest='qtot',type=int, default=0,
                        help = 'total charge on molecule, default: 0 (neutral species)')

    parser.add_argument('-o', dest='file_name', default='charge',
                        help = 'output xyz name (default: charge.xyz)')

    parser.add_argument('-extra', dest='extra',type=float, default=0.1,
                        help = 'extra density factor for placing charges in gap this is only used for test(default: 0.1)')

    parser.add_argument('--test', dest='test', default=False, action='store_const', const=True,
                        help = 'when this falg is on, test mode is True, will report grid pts generated only, this is for getting extra_scale for placing gap charges')


    global args
    # Make relevant inputs lowercase
    args=parser.parse_args(argv)    

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
    for f in args.coord_files:
    
      elements,geo = xyz_parse(f)
      inchikey = GetInchi(elements,geo)
      geo_dict_in[inchikey] = {}
      geo_dict_in[inchikey]["elements"] = elements
      geo_dict_in[inchikey]["geo"] = geo
    
    total_pts = fun(option,geo_dict_in)
    
    return total_pts 



def fun(geo_dict):
    
    # Initialize dictionary and lists for holding modes and geometry information
    bonds       = []
    bond_types  = []
    bond_files  = []
    angles      = []
    angle_types = []
    angle_files = []
    dihedrals   = []
    dihedral_types = []
    dihedral_files = []

    for xyz_file in geo_dict:        
        inchikey_in = xyz_file

        # Print the name
        print( "  {}".format(xyz_file))


        # Generate adjacency table
        geo_dict[xyz_file]["adj_mat"] = Table_generator(geo_dict[xyz_file]["elements"],geo_dict[xyz_file]["geo"])
        

        # Find Hybridizations (needed to identify EZ steroisomers in id_types call)
        hybridizations = Hybridization_finder(geo_dict[xyz_file]["elements"],geo_dict[xyz_file]["adj_mat"])

        # Find atom types
        geo_dict[xyz_file]["atom_types"] = id_types(geo_dict[xyz_file]["elements"],geo_dict[xyz_file]["adj_mat"],args.gens)


        # Determine bonding matrix for the compound
        geo_dict[xyz_file]["bond_mat"] = find_lewis(geo_dict[xyz_file]["elements"], geo_dict[xyz_file]["adj_mat"], q_tot=args.qtot, b_mat_only=True,verbose=False)

        # Find modes (returns all modes in canonicalized format)
        bonds_tmp,angles_tmp,dihedrals_tmp,one_fives_tmp,bond_types_tmp,angle_types_tmp,dihedral_types_tmp,one_five_types_tmp = Find_modes(geo_dict[xyz_file]["adj_mat"],geo_dict[xyz_file]["atom_types"],geo_dict[xyz_file]["bond_mat"],return_all=1)

        if(args.test):
            test_pt = grid_gen(geo_dict[xyz_file]["elements"],geo_dict[xyz_file]["geo"],6,args.size_factor,args.extra)
            return

        scale_factors = [ 13, 22, 30, 50, 60 ]
        density_scales = [ 6.0, 1.1, 1.3, 0.6, 0.2]
        grid_points = {}
        total_points = 0
        for count_i,i in enumerate(scale_factors):
            grid_points[i] = grid_gen(geo_dict[xyz_file]["elements"],geo_dict[xyz_file]["geo"],args.density,scale_factors[count_i],density_scales[count_i])
            total_points += len(grid_points[i])

        # additional density scale factor, this is by trial and error to get an approximate grid pts that are similar to the final # of charges placed reported by Roux's group
        #convert to Bohr
        xyzscale = 1.889725989
        with open('map.xyz','w') as f:
            f.write("{}\n".format(total_points))
            for key in grid_points:
               for i in grid_points[key]:
                  f.write("{} {} {}\n".format(i[0]*xyzscale,i[1]*xyzscale,i[2]*xyzscale))

        with open('connolly_vis.xyz','w') as f:
            f.write("{} \n\n".format(total_points+len(geo_dict[xyz_file]["elements"])))
            for count_i,i in enumerate(geo_dict[xyz_file]["geo"]):
               f.write("{} {} {} {}\n".format(geo_dict[xyz_file]["elements"][count_i],i[0],i[1],i[2]))
            tmp = ['Li','Na','F', 'Cl', 'Br']
            for count_k,key in enumerate(grid_points):
               for i in grid_points[key]:
                  f.write("{} {} {} {}\n".format(tmp[count_k],i[0],i[1],i[2]))
            


    return total_points 



def remove_contact(grid_points,intersec):
    if len(grid_points) == 0 or len(intersec) == 0:
      return grid_points
    Dist_Mat = np.triu(cdist(grid_points,intersec))
    rm_ind,y_ind = np.where( (Dist_Mat < 1.5) & (Dist_Mat > 0.0))
    rm_ind = list(dict.fromkeys(rm_ind))
    
    print("{} charges removed due to close contact w/ already placed charge".format(len(rm_ind)))
    grid_points = [ i for count_i,i in enumerate(grid_points) if count_i not in rm_ind]


    return grid_points

def place_all(criteria,points,geo,bonds,hybridizations,element,density,size_factor,extra_d_scale,name):

    if(os.path.isdir("pointcharge")):
         shutil.move("pointcharge","pointcharge.backup")   
    os.makedirs("pointcharge") 
    os.chdir("pointcharge")
    # place charges along the bonds
    print("Placing charges along the bonds...")
    intersec1 = place_bonds(criteria,points,geo,bonds,element)
    # place charges for lone pairs
    print("Placing charges along lone pairs...")
    intersec2 = place_lone_pair(criteria,points,geo,bonds,hybridizations,element)
    # rm charges placed for lone_pair that are too close to charges placed for bonds
    intersec2 = remove_contact(intersec2,intersec1)
    # place charges in the gap
    intersec3 = place_gap(element,geo,size_factor,extra_d_scale,intersec1,intersec2)

    total_charges = len(intersec1)+len(intersec2)+len(intersec3)
      
    print("Total {} charges placed!".format(total_charges))

    name = name + '.xyz'
    with open(name,'w') as f:
         f.write("{} \n\n".format(len(geo)+total_charges))
         for count_i,i in enumerate(geo):
            f.write("{} {} {} {}\n".format(element[count_i],i[0],i[1],i[2]))
         for i in intersec1:
            f.write("F {} {} {} \n".format(i[0],i[1],i[2]))
         for i in intersec2:
            f.write("Cl {} {} {} \n".format(i[0],i[1],i[2]))
         for i in intersec3:
            f.write("Br {} {} {} \n".format(i[0],i[1],i[2]))
 
    count = 0 
    for i in intersec1:
         filename = 'pointcharges.{}.pc'.format(count)
         with open(filename,'w') as f:
            f.write("0.5 {} {} {} \n".format(i[0],i[1],i[2]))
         count += 1
    for i in intersec2:
         filename = 'pointcharges.{}.pc'.format(count)
         with open(filename,'w') as f:
            f.write("0.5 {} {} {} \n".format(i[0],i[1],i[2]))
         count += 1
    for i in intersec3:
         filename = 'pointcharges.{}.pc'.format(count)
         with open(filename,'w') as f:
            f.write("0.5 {} {} {} \n".format(i[0],i[1],i[2]))
         count += 1
    os.chdir("..")
    
            
   #with open('pointcharges.pc','w') as f:
   #      f.write("{} \n".format(total_charges))
   #      for i in intersec1:
   #         f.write("0.5 {} {} {} \n".format(i[0],i[1],i[2]))
   #      for i in intersec2:
   #         f.write("0.5 {} {} {} \n".format(i[0],i[1],i[2]))
   #      for i in intersec3:
   #         f.write("0.5 {} {} {} \n".format(i[0],i[1],i[2]))
   
    


def place_gap(element,geo,size_factor,extra_d_scale,intersec1,intersec2):
   
    density = 6 # this doesn't matter cuz extra_d_scale controls it
    grid_points = grid_gen(element,geo,density,size_factor,extra_d_scale)
    print("Generting initial gap grids... {} pts".format(len(grid_points)))
    grid_points = remove_contact(grid_points,intersec1)
    grid_points = remove_contact(grid_points,intersec2)

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

# Wrapper function for write commands for *modes files
def write_modelist(name,modes,bond_mat=[]):
    with open(name,'w') as f:
        for count_i,i in enumerate(modes.keys()):
            if len(i) == 2:
                modetype='bond'
            elif len(i) == 3:
                modetype='angle'
            elif len(modes[i]["modes"][0]) == 4 and 2 in [ j[modes[i]["modes"][0][1],modes[i]["modes"][0][2]] for j in bond_mat ]:
                modetype='harmonic_dihedral'
            elif len(modes[i]["modes"][0]) == 4:
                modetype='dihedral'                
            f.write("\n{} start\n".format(modetype))
            f.write('{}\n'.format(len(modes[i]["modes"])))
            f.write("{}\n".format(" ".join([ "{:<60s}".format("_".join([ str(k) for k in j])) for j in modes[i]["modes"] ])))
            for count_j,j in enumerate(modes[i]["atomtypes"][0]):
                f.write('{}\n'.format(" ".join([ "{:<60s}".format(modes[i]["atomtypes"][k][count_j]) for k in range(len(modes[i]["atomtypes"])) ])))
            f.write("{} end\n".format(modetype))
    return 

# Returns the geometry and relevant property lists corresponding to the
# smallest fragment that is consistent with the mode being parametrized.
# If the 1 or 4 atom is a ring, then the entire ring is included. 

# M is mode
def mode_frag(M,Geometry,Adj_mat,Bond_mat,Elements,Atom_types,gens,q_tot=0,keep_terminal=False,keep_rings=False,keep_types=True):

    # Initialize mass_dict (used for identifying the dihedral among a coincident set that will be explicitly scanned)
    if hasattr(mode_frag,"mass_dict") is False:
        mode_frag.mass_dict = {'H':1.00794,'He':4.002602,'Li':6.941,'Be':9.012182,'B':10.811,'C':12.011,'N':14.00674,'O':15.9994,'F':18.9984032,'Ne':20.1797,\
                               'Na':22.989768,'Mg':24.3050,'Al':26.981539,'Si':28.0855,'P':30.973762,'S':32.066,'Cl':35.4527,'Ar':39.948,\
                               'K':39.0983,'Ca':40.078,'Sc':44.955910,'Ti':47.867,'V':50.9415,'Cr':51.9961,'Mn':54.938049,'Fe':55.845,'Co':58.933200,'Ni':58.6934,'Cu':63.546,'Zn':65.39,\
                               'Ga':69.723,'Ge':72.61,'As':74.92159,'Se':78.96,'Br':79.904,'Kr':83.80,\
                               'Rb':85.4678,'Sr':87.62,'Y':88.90585,'Zr':91.224,'Nb':92.90638,'Mo':95.94,'Tc':98.0,'Ru':101.07,'Rh':102.90550,'Pd':106.42,'Ag':107.8682,'Cd':112.411,\
                               'In':114.818,'Sn':118.710,'Sb':121.760,'Te':127.60,'I':126.90447,'Xe':131.29,\
                               'Cs':132.90545,'Ba':137.327,'La':138.9055,'Hf':178.49,'Ta':180.9479,'W':183.84,'Re':186.207,'Os':190.23,'Ir':192.217,'Pt':195.078,'Au':196.96655,'Hg':200.59,\
                               'Tl':204.3833,'Pb':207.2,'Bi':208.98038,'Po':209.0,'At':210.0,'Rn':222.0}

    # Check the consistency of the supplied arguments
    if gens < 1: print( "ERROR in mode_frag: gens variable must be an integer >= 1. Exiting..."); quit();
    if keep_terminal not in [ True, False]: print( "ERROR in mode_frag: keep_terminal argument must be set to a bool. Exiting..."); quit();
    if keep_rings not in [ True, False]: print( "ERROR in mode_frag: keep_rings argument must be set to a bool. Exiting..."); quit();

    # Generate the fragment geometry/attribute lists corresponding to the mode
    N_M,N_Geometry,N_Adj_mat,tmp = mode_geo(M,Geometry,Adj_mat,gens,dup=[Elements,Atom_types])
    N_Elements,N_Atom_types = tmp

    # Preservation protocol depends on the mode type the fragment was generated for.
    fixed_bonds = []
    if len(N_M) == 1:
        loop_list = N_M
    elif len(N_M) == 2:
        loop_list = N_M

        # Updated fixed bonds based on the resonance structure with the largest bond order
        bmat_ind = [ i[M[0],M[1]] for i in Bond_mat ]
        bmat_ind = bmat_ind.index(max(bmat_ind))
        fixed_bonds += [(N_M[0],N_M[1],int(Bond_mat[bmat_ind][M[0],M[1]]))]
    elif len(N_M) == 3:
        loop_list = [N_M[1]]

        # Updated fixed bonds based on the resonance structure with the largest bond order of the bonds involved in the bend
        bmat_ind = [ i[M[0],M[1]]+i[M[1],M[2]] for i in Bond_mat ]
        bmat_ind = bmat_ind.index(max(bmat_ind))
        fixed_bonds += [(N_M[0],N_M[1],int(Bond_mat[bmat_ind][M[0],M[1]]))]
        fixed_bonds += [(N_M[1],N_M[2],int(Bond_mat[bmat_ind][M[1],M[2]]))]
    elif len(N_M) == 4:
        if sum([ N_Adj_mat[N_M[0],j] for j in N_M[1:] ]) == 3: 
            print( "WARNING: USING IMPROPER CRITERIA ON {} CHECK THE RESULT".format([ N_Atom_types[i] for i in N_M ]))
            loop_list = N_M[0]
        else:
            loop_list = N_M[1:3]

            # Updated fixed bonds based on the resonance structure with the largest bond order
            bmat_ind = [ i[M[1],M[2]] for i in Bond_mat ]
            bmat_ind = bmat_ind.index(max(bmat_ind))            
            fixed_bonds += [(N_M[1],N_M[2],int(Bond_mat[bmat_ind][M[1],M[2]]))]
    else:
        print( "ERROR in mode_frag: Protocol doesn't exist for modes involving more than 4 atoms.")
        quit()

    # Include the atoms in the mode and connected atoms within the preserve list.
    preserve = []
    gs = graph_seps(N_Adj_mat)
    for i in loop_list:
        preserve += [ count_i for count_i,i in enumerate(gs[i]) if i < gens ]
    preserve = set(preserve)

    # Perform Hydrogenation    
    N_Geometry,tmp_Atom_types,N_Elements,N_Adj_mat,added_idx = add_hydrogens(N_Geometry,N_Adj_mat,deepcopy(N_Atom_types),N_Elements,q_tot=q_tot,preserve=preserve,fixed_bonds=fixed_bonds)

    # Update the link types
    N_Atom_types += [ 'link-'+tmp_Atom_types[i] for i in added_idx ]

    # Canonicalize by sorting the elements based on hashing (NOTE: range(len(N_atom_types)) is used here rather than "atoms" as in the keep_terminal is True option. 
    Masses = [ mode_frag.mass_dict[N_Elements[i]] for i in range(len(N_Atom_types)) ]
    hash_list,atoms = [ list(j) for j in zip(*sorted([ (atom_hash(count_i,N_Adj_mat,Masses),i) for count_i,i in enumerate(range(len(N_Atom_types))) ],reverse=True)) ]

    # Update lists/arrays based on the sorted atoms
    N_M          = tuple([ atoms.index(i) for i in N_M ])
    N_Geometry   = N_Geometry[atoms]
    N_Adj_mat    = N_Adj_mat[atoms]
    N_Adj_mat    = N_Adj_mat[:,atoms]
    N_Elements   = [ N_Elements[i] for i in atoms ]
    N_Atom_types = [ N_Atom_types[i] for i in atoms ]

    # If keep_types is False then the atomtypes are recalculated here
    if keep_types is False:
        N_Atom_types = id_types(N_Elements,N_Adj_mat,gens=gens)

    return N_M,N_Geometry,N_Adj_mat,N_Elements,N_Atom_types,hash_list

# Description: This is a simple implementation of the Dijkstra algorithm for 
#              finding the backbone of a polymer 
def Dijkstra(Adj_mat):

    # Remove terminal sites (sites with only a single length 2
    # self walk). Continue until all of the terminal structure 
    # has been removed from the topology.
    Adj_trimmed = copy(Adj_mat)

    # Initialize Distances, Previous, and Visited lists    
    Distances = np.array([0] + [100000]*(len(Adj_mat)-1)) # Holds shortest distance to origin from each site
    Previous = np.array([0]*len(Adj_mat)) # Holds the previous site on the short distance to origin
    Visited = [0]*len(Adj_mat) # Holds which sites have been visited

    # Initialize current site and neighbors list
    i = 0 # current site
    neighbors = []

    # Iterate through sites. At each step find the shortest distance between all of hte 
    # current sites neighbors and the START. Update the shortest distances of all sites
    # and choose the next site to iterate on based on which has the shortest distance of
    # among the UNVISITED set of sites. Once the terminal site is identified the shortest
    # path has been found
    while( 0 in Visited):

        # If the current site is the terminal site, then the algorithm is finished
        if i == len(Adj_trimmed)-1:
            break

        # Add new neighbors to the list
        neighbors = [ count_j for count_j,j in enumerate(Adj_trimmed[i]) if j == 1 ]

        # Remove the current site from the list of unvisited
        Visited[i] = 1

        # Iterate over neighbors and update shortest paths
        for j in neighbors:

            # Update distances for current generation of connections
            if Distances[i] + Adj_trimmed[j,i] < Distances[j]:
                Distances[j] = Distances[i] + Adj_trimmed[j,i]
                Previous[j] = i
        
        # Find new site 
        tmp = min([ j for count_j,j in enumerate(Distances) if Visited[count_j] == 0])
        i = [ count_j for count_j,j in enumerate(Distances) if j == tmp and Visited[count_j] == 0 ][0]
        
    # Find shortest path by iterating backwards
    # starting with the end site.
    Shortest_path = [len(Adj_trimmed)-1]    
    i=len(Adj_trimmed)-1

    while( i != 0):
        Shortest_path = Shortest_path + [Previous[i]]    
        i = Previous[i]

    # Reverse order of the list to go from start to finish
    Shortest_path = Shortest_path[::-1]
    return Shortest_path

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

# Adds parameters from the FF file(s) to the FF_dict.
def parse_FF_params(FF_files,FF_dict={"masses":{},"charges":{},"bonds":{},"angles":{},"dihedrals":{},"dihedrals_harmonic":{},"vdw":{}}):
                   
    modes_from_FF = []
    for i in FF_files:
        with open(i,'r') as f:
            for lines in f:
                fields = lines.split()
                if len(fields) == 0: continue
                if fields[0].lower() == "atom":   FF_dict["masses"][fields[1]] = float(fields[3])
                if fields[0].lower() == "charge": FF_dict["charges"][fields[1]] = float(fields[2])
                if fields[0].lower() == "bond":   
                    modes_from_FF += [(fields[1],fields[2])]
                    modes_from_FF += [(fields[2],fields[1])]
                    FF_dict["bonds"][(fields[1],fields[2])] = [fields[3],float(fields[4]),float(fields[5])]
                if fields[0].lower() == "angle":
                    modes_from_FF += [(fields[1],fields[2],fields[3])]
                    modes_from_FF += [(fields[3],fields[2],fields[1])]
                    FF_dict["angles"][(fields[1],fields[2],fields[3])] = [fields[4],float(fields[5]),float(fields[6])]
                if fields[0].lower() in ["dihedral","torsion"]: 
                    modes_from_FF += [(fields[1],fields[2],fields[3],fields[4])]
                    modes_from_FF += [(fields[4],fields[3],fields[2],fields[1])]
                    if fields[5] == "opls":       
                        FF_dict["dihedrals"][(fields[1],fields[2],fields[3],fields[4],fields[5])] = [fields[5]] + [ float(i) for i in fields[6:10] ]
                    elif fields[5] == "harmonic":
                        FF_dict["dihedrals_harmonic"][(fields[1],fields[2],fields[3],fields[4],fields[5])] = [fields[5]] + [ float(fields[6]),int(float(fields[7])),int(float(fields[8])) ] 
                    elif fields[5] == "quadratic":
                        FF_dict["dihedrals_harmonic"][(fields[1],fields[2],fields[3],fields[4],fields[5])] = [fields[5]] + [ float(fields[6]),float(fields[7]) ]
                if fields[0].lower() == "vdw":    
                    FF_dict["vdw"][(fields[1],fields[2])] = [fields[3],float(fields[4]),float(fields[5])]
                    FF_dict["vdw"][(fields[2],fields[1])] = [fields[3],float(fields[4]),float(fields[5])]

    return FF_dict,modes_from_FF

# # This function handles remapping the atomtypes of one fragment on to another such that the modes (bonds
# # angles, dihedrals, etc) in the remapped fragment are the same after remapping. 
def map_frag(fixed_atomtypes,fixed_adj_mat,fixed_bond_mats,fixed_hash,remapped_atomtypes,remapped_adj_mat,remapped_bond_mats,remapped_hash):
    
    # The equivalency of the two topologies is based upon the mode types and number each topology possesses
    bond_types_0,angle_types_0,dihedral_types_0,one_five_types_0 = Find_modes(remapped_adj_mat,remapped_atomtypes,remapped_bond_mats,return_all=1)[4:]
    comp_obj = ( sorted(bond_types_0),sorted(angle_types_0),sorted(dihedral_types_0) )

    # Seed the mapping by placing the first atom in the mappable set onto the first instance 
    # of its hash type in the fixed set.
    #first = len(remapped_hash)-1
    first = 0
    mapping   = {first:fixed_hash.index(remapped_hash[first])}  # fixed indices keyed to remapped indices
    R_mapping = {fixed_hash.index(remapped_hash[first]):first}  # remapped indices keyed to fixed indices
    remapped_place_me = set(range(len(remapped_atomtypes)))     # remapped atom indices that need to be placed
    fixed_place_me = set(range(len(remapped_atomtypes)))        # fixed atom indices that haven't found a match in the remapped set
    remapped_place_me.remove(first)                             # remove the starting atom index
    fixed_place_me.remove(mapping[first])                       # remove the starting atom index

    # Seed the connections list (idx in fixed, connection restraint in remapped) and iterate over the connections until no new connections are found
    cons=[ (count_i,first) for count_i,i in enumerate(fixed_adj_mat[mapping[first]]) if i == 1 and count_i in fixed_place_me ]
    for i in cons:        
        if i[0] not in fixed_place_me: continue
        remapped_match = next( j for j in remapped_place_me if remapped_hash[j] == fixed_hash[i[0]] and remapped_adj_mat[i[1],j] == 1 )
        mapping[remapped_match] = i[0]
        R_mapping[i[0]] = remapped_match
        fixed_place_me.remove(i[0])
        remapped_place_me.remove(remapped_match)
        cons += [ (count_j,remapped_match) for count_j,j in enumerate(fixed_adj_mat[i[0]]) if j == 1 and count_j in fixed_place_me ] 

    # Remap the atomtypes
    new_atomtypes = [ "X" ]*len(remapped_atomtypes)
    for i in mapping: new_atomtypes[mapping[i]] = remapped_atomtypes[i]

    # Check if the set of types and number of each mode match the original fragment
    # NOTE: the fixed_adj_mat is used for determining the connectivity
    # NOTE: the fixed_bond_mats is used for determining dihedral types
    bond_types_1,angle_types_1,dihedral_types_1,one_five_types_1 = Find_modes(fixed_adj_mat,new_atomtypes,fixed_bond_mats,return_all=1)[4:]
    if comp_obj != ( sorted(bond_types_1),sorted(angle_types_1),sorted(dihedral_types_1)):
        print( "ERROR in map_frag: No mapping was discovered that preserved the topology of the original fragment. Check that the two fragments are indeed identical. Exiting...")
        print( "{:60s} {:60s}".format("Original","New"))
        for count_i,i in enumerate(remapped_atomtypes):
            print( "{:<3d} {:60s} {:60s}".format(count_i,i,new_atomtypes[count_i]))
        return new_atomtypes,mapping,1
    
    return new_atomtypes,mapping,0        

# Add hydrogens based upon the supplied atom types. 
# This function is only compatible with TAFFI atom types
# NOTE: Hydrogenation heuristics for geometry assume carbon behavior. This isn't usually a problem when the results are refined with transify, but more specific rules should be implemented in the future
def add_hydrogens(geo,adj_mat,atomtypes,elements,q_tot=0,preserve=[],saturate=True,retype=True,fixed_bonds=[]):

    # Initialize the saturation dictionary the first time this function is called
    if not hasattr(add_hydrogens, "sat_dict"):
        add_hydrogens.sat_dict = {  'H':1, 'He':1,\
                                   'Li':1, 'Be':2,                                                                                                                'B':3,     'C':4,     'N':3,     'O':2,     'F':1,    'Ne':1,\
                                   'Na':1, 'Mg':2,                                                                                                               'Al':3,    'Si':4,     'P':3,     'S':2,    'Cl':1,    'Ar':1,\
                                    'K':1, 'Ca':2, 'Sc':None, 'Ti':None,  'V':None, 'Cr':None, 'Mn':None, 'Fe':None, 'Co':None, 'Ni':None, 'Cu':None, 'Zn':None, 'Ga':None, 'Ge':None, 'As':None, 'Se':None, 'Br':1,    'Kr':None,\
                                   'Rb':1, 'Sr':2,  'Y':None, 'Zr':None, 'Nb':None, 'Mo':None, 'Tc':None, 'Ru':None, 'Rh':None, 'Pd':None, 'Ag':None, 'Cd':None, 'In':None, 'Sn':None, 'Sb':None, 'Te':None,  'I':1,    'Xe':None,\
                                   'Cs':1, 'Ba':2, 'La':None, 'Hf':None, 'Ta':None,  'W':None, 'Re':None, 'Os':None, 'Ir':None, 'Pt':None, 'Au':None, 'Hg':None, 'Tl':None, 'Pb':None, 'Bi':None, 'Po':None, 'At':None, 'Rn':None  }

        add_hydrogens.lone_e = {    'H':0, 'He':2,\
                                   'Li':0, 'Be':2,                                                                                                                'B':0,     'C':0,     'N':2,     'O':4,     'F':6,    'Ne':8,\
                                   'Na':0, 'Mg':2,                                                                                                               'Al':0,    'Si':0,     'P':2,     'S':4,    'Cl':6,    'Ar':8,\
                                    'K':0, 'Ca':2, 'Sc':None, 'Ti':None,  'V':None, 'Cr':None, 'Mn':None, 'Fe':None, 'Co':None, 'Ni':None, 'Cu':None, 'Zn':None, 'Ga':None, 'Ge':0,    'As':3,    'Se':4,    'Br':6,    'Kr':None,\
                                   'Rb':0, 'Sr':2,  'Y':None, 'Zr':None, 'Nb':None, 'Mo':None, 'Tc':None, 'Ru':None, 'Rh':None, 'Pd':None, 'Ag':None, 'Cd':None, 'In':None, 'Sn':None, 'Sb':None, 'Te':None,  'I':6,    'Xe':None,\
                                   'Cs':0, 'Ba':2, 'La':None, 'Hf':None, 'Ta':None,  'W':None, 'Re':None, 'Os':None, 'Ir':None, 'Pt':None, 'Au':None, 'Hg':None, 'Tl':None, 'Pb':None, 'Bi':None, 'Po':None, 'At':None, 'Rn':None  }

        add_hydrogens.frag = 0

    # Intermediate scalars
    H_length = 1.1
    N_atoms  = len(geo)
    init_len = len(geo)

    # If the user specifies a set of atoms to preserve as is, then
    # then bonding_pref entry is set to full saturation.
    if preserve != []:
        bonding_pref = [ (i,add_hydrogens.sat_dict[elements[i]]) for i in preserve ]
    else:
        bonding_pref = None

    # Get the lewis structure
    lone_electrons,bonding_electrons,core_electrons,bonding_pref = check_lewis(atomtypes,adj_mat,q_tot=q_tot,bonding_pref=bonding_pref,return_pref=True,fixed_bonds=fixed_bonds)

    # Update the preserved atoms (check_lewis will extend this list if there are special groups (e.g., nitro moieties) that need to be conserved for the sake of the lewis structure)
    preserve = set([ i[0] for i in bonding_pref ])

    # Loop over the atoms in the geometry
    for count_i,i in enumerate(geo):
        
        # ID undercoordinated atoms
        if count_i in preserve:
            continue
        elif add_hydrogens.sat_dict[elements[count_i]] is not None:
            B_expected = add_hydrogens.sat_dict[elements[count_i]]
        else:
            print( "ERROR in add_hydrogens: could not determine the number of hydrogens to add to {}. Exiting...").format(elements[count_i])
            quit()
        B_current  = bonding_electrons[count_i]

        # Determine the number of nuclei that are attached and expected.
        N_current   = sum(adj_mat[count_i])
        N_expected = N_current + (B_expected - B_current)

        # Add hydrogens to undercoordinated atoms
        if N_expected > N_current:
            
            old_inds = [ count_j for count_j,j in enumerate(adj_mat[count_i]) if j == 1 ]

            # Protocols for 1 missing hydrogen
            if N_expected - N_current == 1:
                if N_expected == 1:
                    new = i + np.array([H_length,0.0,0.0])
                elif N_expected == 2:
                    new = -1.0 * normalize(geo[old_inds[0]] - i) * H_length + i + np.array([random.random(),random.random(),random.random()])*0.01 #random factor added for non-carbon types to relax during FF-opt
                elif N_expected == 3:
                    new = -1.0 * normalize( normalize(geo[old_inds[0]] - i) + normalize(geo[old_inds[1]] - i) ) * H_length + i
                elif N_expected == 4:
                    print( "WARNING in add_hydrogens: 1:4 (unchecked)")
                    new = -1.0 * normalize( normalize(geo[old_inds[0]] - i) + normalize(geo[old_inds[1]] - i) + normalize(geo[old_inds[2]] - i) ) * H_length + i                

                # Update geometry, adj_mat, elements, and atomtypes with one new atoms
                geo = np.vstack([geo,new])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])]
                elements += ["H"]
                tmp = np.zeros([N_atoms+1,N_atoms+1])
                tmp[:N_atoms,:N_atoms] = adj_mat
                tmp[-1,count_i] = 1
                tmp[count_i,-1] = 1
                adj_mat = tmp                
                N_atoms += 1

            # Protocols for 2 missing hydrogens
            # ISSUE, NEW ALGORITHM IS BASED ON BONDED ATOMS NOT BONDED CENTERS
            if N_expected - N_current == 2:
                if N_expected == 2:
                    new_1 = i + np.array([H_length,0.0,0.0])
                    new_2 = i - np.array([H_length,0.0,0.0])
                elif N_expected == 3:
                    rot_vec = normalize(np.cross( geo[old_inds[0]] - i, np.array([random.random(),random.random(),random.random()]) ))
                    new_1 = normalize(axis_rot(geo[old_inds[0]],rot_vec,i,120.0) - i)*H_length + i
                    new_2 = normalize(axis_rot(geo[old_inds[0]],rot_vec,i,240.0) - i)*H_length + i
                elif N_expected == 4:
                    bisector = normalize(geo[old_inds[0]] - i + geo[old_inds[1]] - i) 
                    new_1    = axis_rot(geo[old_inds[0]],bisector,i,90.0)
                    new_2    = axis_rot(geo[old_inds[1]],bisector,i,90.0) 
                    rot_vec  = normalize(np.cross(new_1-i,new_2-i))
                    angle    = ( 109.5 - acos(np.dot(normalize(new_1-i),normalize(new_2-i)))*180.0/np.pi ) / 2.0
                    new_1    = axis_rot(new_1,rot_vec,i,-angle)
                    new_2    = axis_rot(new_2,rot_vec,i,angle)
                    new_1    = -1*H_length*normalize(new_1-i) + i
                    new_2    = -1*H_length*normalize(new_2-i) + i
                    
                # Update geometry, adj_mat, elements, and atomtypes with two new atoms
                geo = np.vstack([geo,new_1])
                geo = np.vstack([geo,new_2])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])]*2
                elements += ["H","H"]
                tmp = np.zeros([N_atoms+2,N_atoms+2])
                tmp[:N_atoms,:N_atoms] = adj_mat
                tmp[[-1,-2],count_i] = 1
                tmp[count_i,[-1,-2]] = 1
                adj_mat = tmp
                N_atoms += 2

            # Protocols for 3 missing hydrogens
            if N_expected - N_current == 3:
                if N_expected == 3:
                    rot_vec = np.array([0.0,1.0,0.0])
                    new_1 = i + np.array([H_length,0.0,0.0])
                    new_2 = axis_rot(new_1,rot_vec,i,120.0)
                    new_3 = axis_rot(new_1,rot_vec,i,240.0)
                if N_expected == 4:
                    rot_vec = normalize(np.cross( geo[old_inds[0]] - i, np.array([random.random(),random.random(),random.random()]) ))
                    new_1 = H_length*normalize(axis_rot(geo[old_inds[0]],rot_vec,i,109.5)-i) + i
                    new_2 = axis_rot(new_1,normalize(i-geo[old_inds[0]]),i,120.0)
                    new_3 = axis_rot(new_2,normalize(i-geo[old_inds[0]]),i,120.0)

                # Update geometry, adj_mat, elements, and atomtypes with three new atoms
                geo = np.vstack([geo,new_1])
                geo = np.vstack([geo,new_2])
                geo = np.vstack([geo,new_3])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])]*3
                elements += ["H","H","H"]
                tmp = np.zeros([N_atoms+3,N_atoms+3])
                tmp[:N_atoms,:N_atoms] = adj_mat
                tmp[[-1,-2,-3],count_i] = 1
                tmp[count_i,[-1,-2,-3]] = 1
                adj_mat = tmp
                N_atoms += 3

            # Protocols for 4 missing hydrogens
            if N_expected - N_current == 4:
                if N_expected == 4:
                    new_1 = i + np.array([H_length,0.0,0.0])
                    rot_vec = normalize(np.cross( new_1 - i, np.array([random.random(),random.random(),random.random()]) ))
                    new_2 = H_length*normalize(axis_rot(new_1,rot_vec,i,109.5)-i) + i
                    new_3 = axis_rot(new_2,normalize(i-new_1),i,120.0)
                    new_4 = axis_rot(new_3,normalize(i-new_1),i,120.0)
                    
                # Update geometry, adj_mat, elements, and atomtypes with three new atoms
                geo = np.vstack([geo,new_1])
                geo = np.vstack([geo,new_2])
                geo = np.vstack([geo,new_3])
                geo = np.vstack([geo,new_4])
                atomtypes += ["[1[{}]]".format(atomtypes[count_i].split(']')[0].split('[')[1])]*4
                elements += ["H","H","H","H"]
                tmp = np.zeros([N_atoms+4,N_atoms+4])
                tmp[:N_atoms,:N_atoms] = adj_mat
                tmp[[-1,-2,-3,-4],count_i] = 1
                tmp[count_i,[-1,-2,-3,-4]] = 1
                adj_mat = tmp
                N_atoms += 4

    if retype is True:
        return geo,id_types(elements,adj_mat,gens=2),elements,adj_mat,range(init_len,len(geo))
    else:
        return geo,atomtypes,elements,adj_mat,range(init_len,len(geo))

# Shortcut for normalizing a vector
def normalize(x):
    return x/sum(x**(2.0))**(0.5)

# Recursively generates the permutations of indices corresponding to identical values in types
def gen_perms(types,start=0):
    ind = [ count_i for count_i,i in enumerate(types) if i==types[0] ]
    not_ind = [ count_i for count_i,i in enumerate(types) if i !=types[0] ]
    for i in permutations(range(start,start+len(ind)),len(ind)):
        if len(not_ind) > 0:
            for j in gen_perms([ types[j] for j in not_ind ],start=len(ind)+start):
                yield i+j
        else:
            yield i

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
   
def miss2default(config,other_opt):

 
    # Create default dict
    options = ['FF_db','outputname','gens','recursive_opt','q_tot','keep_terminal_opt','avoid_frags_opt']

    defaults = ['', '.', 2, False, 0, False, False]  

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



