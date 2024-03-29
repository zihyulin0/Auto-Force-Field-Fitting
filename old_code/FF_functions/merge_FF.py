#!/bin/env python                                                                                                                                                             
# Author: Brett Savoie (brettsavoie@gmail.com)

import sys,argparse,os,datetime
from numpy import *
import time
import ast
from scipy import *
from scipy.spatial.distance import *
from numpy.linalg import *
import shutil
import random

def main(argv):


    parser = argparse.ArgumentParser(description='Reads in a geometry and FF.db file (usually outputted from an intramolecular mode parse) '+\
                                                 'and writes inputs for a LAMMPS job to generate configurations for the parsing VDW parameters.')


    #required (positional) arguments                                                                                                  
    parser.add_argument('master', help = 'A text file holding pre-existing FF data, or the name of the new FF database file to be created. '+\
                                         'This file has minimal formatting requirements other than that it should conform to the general formatting '+\
                                         'scheme for param files that are compatible with polygen.py. NOTE: the program expects to find a *.molecules file in the same location as this file. '+\
                                         'The "molecules" file is used to link the FF parameters with the sample molecules that were used to generate the parameters. If "master" does not exist, then' +\
                                         'it is created along with its *.molecules file.')

    parser.add_argument('new_params', help = 'A text file containing FF parameters to be merged with the master database. NOTE: the program expects to find a molecules.db file in the '+\
                                          'same folder as this file. This is a safeguard to ensure that all parameter files are properly generated by the proper machinery, which '+\
                                          'automatically generate this file alongside the parameters files. The molecules.db file is used to link each parameter in the FF database with the '+\
                                          'sample molecule(s) used to generate it.')

    #optional arguments    
    parser.add_argument('-only', dest='only', default='',
                        help = 'A space delimitted string containing the types of FF parameters the user wants to merge with the master force-field database. '+\
                               'For example, "bond angle torsion" would only parse the bonds angles and dihedrals in the input FF file. The arguments should match FF mode keys '+\
                               '(vdw,atom,bond,angle,dihedral,torsion,charge). Arguments supplied to this option override those of -except. When left empty, the arguments supplied to '+\
                               '-except are used. When neither is set, all discovered parameters are parsed. (default: none)')

    parser.add_argument('-no_ct', dest='no_ct', default=False, action='store_const', const=True,
                        help = 'If this flag is present then the cross terms in the vdw entries are not saved. (default: off)')

    parser.add_argument('-except', dest='besides', default='',
                        help = 'A space delimitted string containing the types of FF parameters the user wants to avoid merging with the master force-field database. '+\
                               'For example, "vdw" would only parse all parameters except the vdw discovered in the input FF file. The arguments should match FF mode keys '+\
                               '(vdw,atom,bond,angle,dihedral,torsion,charge). Arguments supplied to this option overriden by those of -only. When neither -except or -only are set,'+\
                               'all discovered parameters are parsed. (default: none)')

    parser.add_argument('-num_bkp', dest='num_bkp', default=5,
                        help = 'Everytime the database is updated an image of the unupdated database and accompanying molfile are copied to .bkp1.*master_name*. In turn, the previous '+\
                               'backup, if it exists, is copied to .bkp2.*master_name* and so on for however many backups the user wants to keep. This is a precaution that allowed updates '+\
                               'to be reversed. By default, 5 updated are kept (default: 5)')

    parser.add_argument('-R', dest='replace', default=0, action='store_const', const=1,
                        help = 'By default, the script will not replace existing database entries, regardless of whether new parameters for the same modes/charges are discovered. When -R is '+\
                               'present, however, the script will replace existing database entries with any new entries discovered.')

    # Make parse inputs
    args=parser.parse_args(argv)


    # Create fake option for function input
    option = dict(args.__dict__.items())

    fun(option)
   


def fun(option):

    option = miss2default(option)

    master = option['master']
    new_params = option['new_params']
    only = option['only']
    no_ct = option['no_ct']
    besides = option['besides']
    num_bkp = option['num_bkp']
    replace = option['replace']

    # Parse input arguments
    if only != '': only = only.split()
    if besides != '': besides = besides.split()
    
    # Initialize list of all supported mode types
    types = [ "atom", "vdw", "bond", "angle", "torsion", "dihedral", "charge" ]    

    if only != '':
        types = [ i for i in types if i in only ]
    elif besides != '':
        types = [ i for i in types if i not in besides ]

    # generate *.molecules.db names
    if len([ i for i in master.split('/')[:-1]]) > 1: master_path = '/'.join([ i for i in master.split('/')[:-1]])
    else: master_path = '.'
    if len([ i for i in new_params.split('/')[:-1]]) > 1: new_path    = '/'.join([ i for i in new_params.split('/')[:-1]])
    else: new_path = '.'

    # Check that the supplied FF files and contingent *.molecules files exist
    if os.path.isfile(new_params) != True:
        print("ERROR: Specified database file ({}) doesn't exist. Please check the path. Exiting...".format(new_params))
        quit()

    if os.path.isfile(master) != True:
        print("WARNING: Specified database file ({}) doesn't exist. Creating a new database...".format(master))
        with open(master,'w') as f:
            f.write("# Database file for FF parameters\n")
            f.write("# formatting notes:\n")
            f.write("#     vdw-lj units:          kcal/mol angstrom (eps,sigma respectively)\n")
            f.write("#     harmonic-bond units:   kcal/mol angstroms (k,r_0 respectively)\n")
            f.write("#     harmonic-angle units:  kcal/mol degrees (k,theta_0 respectively)\n")
            f.write("#     opls-torsion units:    kcal/mol for all (fourier coefficients)\n")
            f.write("#     charge units:          fraction of elementary charge\n")

    # Create backups if they don't already exist
    # NOTE: convoluted "/".join call is to properly parse the path and ".".join call to
    #       to handle names with multiple periods.
    for i in range(num_bkp):        
        if os.path.isfile(master_path+'/.bkp'+str(num_bkp-i)+'.'+master.split('/')[-1]) != True:
            shutil.copyfile(master,master_path+'/.bkp'+str(num_bkp-i)+'.'+master.split('/')[-1])            

    # Read in the master list data
    Master_Data = get_data(master,keep_types = [ "atom", "vdw", "bond", "angle", "torsion", "dihedral", "charge" ],no_ct=no_ct)
    New_Data = get_data(new_params,keep_types=types,no_ct=no_ct)
    
    # Merge new parameters with master    
    added_elements = 0
    updated_elements = 0
    for i in list(New_Data.keys()):
        for j in list(New_Data[i].keys()):

            # Add new parameters that don't exist in current master database.
            if j not in list(Master_Data[i].keys()):
                
                # Add new parameter and append mol_id to the info
                Master_Data[i][j] = New_Data[i][j]
                added_elements += 1

            # Update parameters that exist in current master database (if -R is toggled in command line call).
            elif j in list(Master_Data[i].keys()) and replace == 1:

                # Add new parameter and append mol_id to the info
                Master_Data[i][j] = New_Data[i][j]
                updated_elements += 1
                
    # Write updated master FF 
    if added_elements > 0 or updated_elements > 0:
        if added_elements > 0:
            print("Updating {} with {} new parameters from {}...".format(master,added_elements,new_params))
        if updated_elements > 0:
            print("Updating {} with {} updated parameters from {}...".format(master,updated_elements,new_params))
        write_updated(Master_Data,master,num_bkp)
    else:
        print("No new parameters discovered in {}. Exiting without modifying {}...".format(new_params,master))

    return
    
# Description: Reads in the FF
def get_data(master_file,keep_types=[ "atom", "vdw", "bond", "angle", "torsion", "dihedral", "charge" ],no_ct=False):

    Data = {"atoms":{},"bonds":{},"angles":{},"dihedrals":{},"vdws":{},"charges":{}}

    # Read in Atom data
    with open(master_file,'r') as f:
        for lines in f:
            fields = lines.split()

            if len(fields) == 0: continue
            if fields[0] == "atom" and "atom" in keep_types:
                Data["atoms"][fields[1]] = fields[1:]
            if fields[0] == "vdw" and "vdw" in keep_types:
                if no_ct is True and fields[1] == fields[2]:
                    Data["vdws"][(fields[1],fields[2],fields[3])] = fields[1:]
                elif no_ct is False:
                    Data["vdws"][(fields[1],fields[2],fields[3])] = fields[1:]
            if fields[0] == "bond" and "bond" in keep_types:
                Data["bonds"][(fields[1],fields[2],fields[3])] = fields[1:]
            if fields[0] == "angle" and "angle" in keep_types:
                Data["angles"][(fields[1],fields[2],fields[3],fields[4])] = fields[1:]
            if fields[0] == "dihedral" or fields[0] == "torsion" and ( "dihedral" in keep_types or "torsion" in keep_types ):
                Data["dihedrals"][(fields[1],fields[2],fields[3],fields[4],fields[5])] = fields[1:]
            if fields[0] == "charge" and "charge" in keep_types:
                Data["charges"][fields[1]] = fields[1:]
    
    return Data

def write_updated(Data,master_FF,num_bkp):

    # generate master_path and *.molecules.db name
    if len([ i for i in master_FF.split('/')[:-1]]) > 1: master_path = '/'.join([ i for i in master_FF.split('/')[:-1]])
    else: master_path = './'

    #######################
    #  Write .tmp FF file #
    #######################

    # Copy header    
    with open(master_FF,'r') as f:
        with open(master_FF+'.tmp', 'w') as h:
            for i in f:
                if i[0] != "#": break
                else: h.write(i)

    # The formatting of this database file is compatible 
    # with the input for the polygen.py program.
    with open(master_FF+'.tmp', 'a') as f:

        # Write atom type definitions
        f.write("\n\n# Atom type definitions\n#\n{:<10s} {:<60s} {:<59s} {:<20s} {:<6s}\n"\
                .format("#","Atom_type","Label","Mass", "Mol_ID"))        
        for i in sorted(Data["atoms"].keys()):
            if len(Data["atoms"][i]) > 3:
                f.write("{:<10s} {:<60s} {:<59s} {:<20s} {:<s}\n"\
                        .format("atom",Data["atoms"][i][0],Data["atoms"][i][1],Data["atoms"][i][2],Data["atoms"][i][3]))
            elif len(Data["atoms"][i]) == 3:
                f.write("{:<10s} {:<60s} {:<59s} {:<20s}\n"\
                        .format("atom",Data["atoms"][i][0],Data["atoms"][i][1],Data["atoms"][i][2]))                

        # Write VDW definitions
        f.write("\n# VDW definitions\n#\n{:<10s} {:<60s} {:<60s} {:<15s} {:<41s} {:<6s}\n"\
                .format("#","Atom_type","Atom_type","Potential","params (style determines #args)","Mol_ID"))
        for i in sorted(Data["vdws"].keys()):
            f.write("{:<10s} {:<60s} {:<60s} {:<15s} {} {:<s}\n"\
                    .format("vdw",Data["vdws"][i][0],Data["vdws"][i][1],Data["vdws"][i][2]," ".join([ "{:<20s}".format(j) for j in Data["vdws"][i][3:-1] ]),Data["vdws"][i][-1])) 

        # Write bond definitions
        f.write("\n# Bond type definitions\n#\n{:<10s} {:<40s} {:<40s} {:<15s} {:<41s} {:<6s}\n".format("#","Atom_type","Atom_type","style","params (style determines #args)","Mol_ID"))
        for i in sorted(Data["bonds"].keys()):
            f.write("{:<10s} {:<40s} {:<40s} {:<15s} {} {:<s}\n"\
                    .format("bond",Data["bonds"][i][0],Data["bonds"][i][1],Data["bonds"][i][2]," ".join([ "{:<20s}".format(j) for j in Data["bonds"][i][3:-1] ]),Data["bonds"][i][-1]))

        # Write angle definitions
        f.write("\n# Angle type definitions\n#\n{:<10s} {:<40s} {:<40s} {:<40s} {:<15s} {:<41s} {:<6s}\n"\
                .format("#","Atom_type","Atom_type","Atom_type","style","params (style determines #args)","Mol_ID"))
        for i in sorted(Data["angles"].keys()):
            f.write("{:<10s} {:<40s} {:<40s} {:<40s} {:<15s} {} {:<s}\n"\
                    .format("angle",Data["angles"][i][0],Data["angles"][i][1],Data["angles"][i][2],Data["angles"][i][3],\
                            " ".join([ "{:<20s}".format(j) for j in Data["angles"][i][4:-1] ]),Data["angles"][i][-1]))

        # Write dihedral definitions
        f.write("\n# Dihedral/Torsional type definitions\n#\n{:<10s} {:<40s} {:<40s} {:<40s} {:<40s} {:<15s} {:<83s} {:<6s}\n"\
                .format("#","Atom_type","Atom_type","Atom_type","Atom_type","style","params (style determines #args)","Mol_ID"))
        for i in sorted(Data["dihedrals"].keys()):
            f.write("{:<10s} {:<40s} {:<40s} {:<40s} {:<40s} {:<15s} {} {:<s}\n"\
                    .format("torsion",Data["dihedrals"][i][0],Data["dihedrals"][i][1],Data["dihedrals"][i][2],Data["dihedrals"][i][3],Data["dihedrals"][i][4],\
                            " ".join([ "{:<20s}".format(j) for j in Data["dihedrals"][i][5:-1] ]),Data["dihedrals"][i][-1]))

        # Write charge definitions
        f.write("\n# Charge definitions\n#\n{:<10s} {:<61s} {:<6s}\n".format("#","Atom_type","Charge","Mol_ID"))
        for i in sorted(Data["charges"].keys()):
            if len(Data["charges"][i]) > 2:
                f.write("{:<10s} {:<60s} {:<20s} {:<20s}\n".format("charge",Data["charges"][i][0],Data["charges"][i][1],Data["charges"][i][2]))
            elif len(Data["charges"][i]) == 2:
                f.write("{:<10s} {:<60s} {:<20s}\n".format("charge",Data["charges"][i][0],Data["charges"][i][1]))

    # copy backups
    # NOTE: convoluted "/".join call is to properly parse the path and ".".join call to
    #       to handle names with multiple periods.
    for i in range(num_bkp):        
        if i == num_bkp-1:
            current_FF  = master_FF
            old_FF      = master_path+'/.bkp'+str(num_bkp-i)+'.'+master_FF.split('/')[-1]
        else:
            current_FF  = master_path+'/.bkp'+str(num_bkp-i-1)+'.'+master_FF.split('/')[-1]
            old_FF      = master_path+'/.bkp'+str(num_bkp-i)+'.'+master_FF.split('/')[-1]
        shutil.copyfile(current_FF,old_FF)

    # copy .tmp file to master
    # Note: the overwrite of the original data is performed last to
    #       avoid losing information if there is an interuption.
    shutil.copyfile(master_FF+'.tmp',master_FF)

    # remove .tmp files
    os.remove(master_FF+'.tmp')
    return


#set the missing option to default
def miss2default(option): 
   
    options = ['master','new_params','only','no_ct','besides','num_bkp','replace']
    missing = [ i in option for i in options]
    
    # set default
    default = {}
    default['only'] = ''
    default['no_ct'] = False
    default['besides'] = ''
    default['num_bkp'] = 5
    default['replace'] = 0
    
    # set missing option to default
    for count_i,i in enumerate(missing):
      if i is False:
         option[options[count_i]] = default[options[count_i]]

    return option


   


if __name__ == "__main__":
   main(sys.argv[1:])
   
