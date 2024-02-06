import os,sys
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/Lib')
sys.path.append('/'.join(os.path.abspath(__file__).split('/')[:-2])+'/FF_functions')
import generate_connolly

def main(argv):
    global user,run_dir,c
    print(argv)
    fitfolder = argv[0] # fitcharge_{}
    q = argv[1] # 0.5, -0.5
    c_config = argv[2] # nonpolar/inchi/charges/CHELPG_calcs/configs
    run_dir = argv[3] # getcwd... in the very original driver.py command

    
    i = fitfolder.split('_')[1]
    xyz = '{}/{}/{}.xyz'.format(c_config,i,i)
    total_pc,pc_layer1,pc_layer2,pc_layer3 = get_pc(fitfolder)
    rewrite_pot_wrap(xyz,fitfolder,total_pc,'_{}'.format(q))

def get_pc(calc_folder):

    ori_dir = os.getcwd()
    print(run_dir+'/'+calc_folder)
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

def rewrite_pot_wrap(xyz,calc_folder,total_pc,tag):
    config_folder = 'configs_connolly{}'.format(tag)
    config_ori_folder = 'configs{}'.format(tag)
    os.chdir('{}/{}/{}'.format(run_dir,calc_folder,config_folder))
    folders =  list(range(0,total_pc)) 
    folders.append('unperturbed')
    for i in folders:
         if os.path.isdir(str(i)) is False: os.makedirs(str(i))
         else: 
            if os.path.isfile('{}/{}_charges.vpot'.format(i,i)) is True: continue
         os.chdir(str(i))

         oldpot = '{}/{}/{}/{}/{}_charges.vpot'.format(run_dir,calc_folder,config_ori_folder,i,i)
         newpot = '{}_charges.vpot'.format(i)
         connolly_grid = get_connolly_grid(calc_folder,xyz)
         rewrite_potential(oldpot,newpot,connolly_grid)
         os.chdir('{}/{}/{}'.format(run_dir,calc_folder,config_folder))
    return

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

def get_connolly_grid(calc_folder,xyz):

    ori_dir = os.getcwd()
    # Genererating connolly surface grid
    os.chdir(run_dir+'/'+calc_folder)
    if os.path.isfile('map.xyz') is False:
         sublist = ('{} -d 1.41'.format(xyz)).split()
         connolly_grid = generate_connolly.main(sublist)
    else:
         with open('map.xyz','r') as f:
            for lc,lines in enumerate(f):
               fields = lines.split()
               if lc == 0:
                  connolly_grid = int(fields[0]) 
                  break
    os.chdir(ori_dir)

    return connolly_grid

if __name__ == "__main__":
   main(sys.argv[1:])
