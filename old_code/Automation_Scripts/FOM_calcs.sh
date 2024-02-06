#!/bin/bash

# Consistency check
FF_FOM=$1
if [ $# -lt 1 ]; then
    echo "ERROR in FOM_calcs.sh: A force-field filename must be supplied to the first argument." 
    exit
else
    FOM_name="FOM"
    NM_name="Normal_Modes"
fi

# Consistency check
if [ $# -lt 2 ]; then
    echo "ERROR in FOM_calcs.sh: A configureation filename must be supplied to the second argument." 
    exit
else
    config_file=$2
fi

# IMPORT SYSTEM SPECIFIC PATHS AND BATCH ARGUMENTS
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"   # Saves the location of this script to a variable
TAFFI_PATH=${DIR}/../                                     # It is easier to define a lot of the commands relative to the taffi root
PYTHON_ABS=$(which python)                                # Saves the absolute path to the favored python dist to variable
ORCA_SUBMIT=${DIR}/orca_submit.py
SHELL_SUBMIT=${DIR}/shell_submit.sh
LAMMPS_SUBMIT=${DIR}/lammps_submit.py
LAMMPS_EXE=$( awk -v v="LAMMPS_EXE" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
ORCA_EXE=$( awk -v v="ORCA_EXE" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )

# FOM SPECIFIC ARGUMENTS
FOM_GEOOPT_PROCS=$( awk -v v="FOM_GEOOPT_PROCS" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
FOM_GEOOPT_WT=$( awk -v v="FOM_GEOOPT_WT" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
FOM_GEOOPT_Q=$( awk -v v="FOM_GEOOPT_Q" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
FOM_GEOOPT_SCHED=$( awk -v v="FOM_GEOOPT_SCHED" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
FOM_GEOOPT_PPN=$( awk -v v="FOM_GEOOPT_PPN" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )  
FOM_GEOOPT_SIZE=$( awk -v v="FOM_GEOOPT_SIZE" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} ) 

FOM_MD_GAS_PROCS=$( awk -v v="FOM_MD_GAS_PROCS" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )   
FOM_MD_GAS_WT=$( awk -v v="FOM_MD_GAS_WT" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )      
FOM_MD_GAS_Q=$( awk -v v="FOM_MD_GAS_Q" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )       
FOM_MD_GAS_NPP=$( awk -v v="FOM_MD_GAS_NPP" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )     
FOM_MD_GAS_SCHED=$( awk -v v="FOM_MD_GAS_SCHED" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )   
FOM_MD_GAS_PPN=$( awk -v v="FOM_MD_GAS_PPN" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )     
FOM_MD_GAS_SIZE=$( awk -v v="FOM_MD_GAS_SIZE" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )     

FOM_MD_COND_PROCS=$( awk -v v="FOM_MD_COND_PROCS" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )  
FOM_MD_COND_WT=$( awk -v v="FOM_MD_COND_WT" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )     
FOM_MD_COND_Q=$( awk -v v="FOM_MD_COND_Q" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )      
FOM_MD_COND_NPP=$( awk -v v="FOM_MD_COND_NPP" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )    
FOM_MD_COND_SCHED=$( awk -v v="FOM_MD_COND_SCHED" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )  
FOM_MD_COND_PPN=$( awk -v v="FOM_MD_COND_PPN" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )    
FOM_MD_COND_SIZE=$( awk -v v="FOM_MD_COND_SIZE" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )   
FOM_TEMPS=$( awk -v v="FOM_TEMPS" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
if [ ! -z ${FOM_TEMPS} -a ! -f ${FOM_TEMPS} ]; then
    echo "ERROR: the supplied temperature file ${FOM_TEMPS} does not exist. Exiting..."
    exit
fi

######################

# Waits for completion of jobs in the queue, expects an 
# array of job ids to be supplied as the first argument
function check_queue () {
    arr=("$@")
    result=('1')
    count=1
    while [ ${#result[@]} -gt 0 ]; do
	tmp=$( ${PYTHON_ABS} ~/qstat.py | grep ${USER} | awk '{ if ( NF > 4 && $5 != "C" ) { print $1 }}')
	readarray -t ids_curr <<< "${tmp}"
	l2=" ${ids_curr[*]} "                    # add framing blanks
	result=()
	for item in ${arr[@]}; do
	    if [[ $l2 =~ " $item " ]] ; then    # use $item as regexp
		result+=($item)
	    fi
	done	
	if [ ${#result[@]} -gt 0 ]; then
	    if [ ${count} == 1 ]; then
		echo "waiting for jobs (${result[*]}) to complete...."
		count=$(( ${count} +  1 ))
	    fi	    
	    sleep 100
	fi
    done
} 

# Inserts escape characters before parentheses so that the script can
# handle filenames/folders that involve these
function esc_str () {
    a=${1}
    a=$( echo ${a} | sed 's/(/\\(/' )
    a=$( echo ${a} | sed 's/)/\\)/' )
    echo ${a}
}

###############################
# Run normal mode comparisons #
###############################
if [ ! -d ${NM_name} ]; then
    mkdir ${NM_name}
    cp *.xyz ${NM_name}/.
fi

# Create geometry optimization runs 
cd ${NM_name}
for i in *.xyz; do
    name=$( echo ${i} | awk -F '[.]' '{ print $1 }' )	
    if [ ! -d ${name}_freq ]; then
	# transify the geometry then geoopt it
	${PYTHON_ABS} ${TAFFI_PATH}/Lib/transify.py ${i}
	obminimize -ff uff -sd -c 1e-20 -n 100000 straightened.xyz > ${name}.opt.xyz 2> /dev/null
	awk 'BEGIN{flag=0}{ if ( flag == 0 && $1 == "WARNING:" ) { next } else { flag = 1 }; if ( flag == 1 ) { print $0 } }' ${name}.opt.xyz >> tmp.${name}.opt.xyz # Remove warnings
	mv tmp.${name}.opt.xyz ${name}.opt.xyz
	${PYTHON_ABS} ${TAFFI_PATH}FF_functions/xyz_to_orca.py ${name}.opt.xyz -p ${FOM_GEOOPT_PROCS} -o ${name}_freq.in
	mkdir ${name}_freq
	mv ${name}_freq.in ${name}_freq/.
	mv ${name}.opt.xyz ${name}_freq/.
	rm straightened.xyz
    fi
done

# Submit jobs (script avoids running any jobs with existing data)
ids_sub=()
ids_sub+=$(${PYTHON_ABS} ${ORCA_SUBMIT} -f "*freq.in" -p ${FOM_GEOOPT_PROCS} -t ${FOM_GEOOPT_WT} -ppn ${FOM_GEOOPT_PPN} -q ${FOM_GEOOPT_Q} -sched ${FOM_GEOOPT_SCHED} -o frequencies -path_to_exe ${ORCA_EXE} -size ${FOM_GEOOPT_SIZE} --silent)$'\n'

# Wait until jobs complete
check_queue "${ids_sub[@]}"	

echo "COMPLETE: Reference frequency calculations."

# Remove job submission files
if [ -f frequencies.0.submit ]; then	
    rm frequencies.*
fi

# Calculate the classical modes
ids_sub=()
submits_list=()
for i in *.xyz; do

    name=$( echo ${i} | awk -F '[.]' '{ print $1 }' )	
    name_sub=$( echo ${name} | sed 's/(/-/' )
    name_sub=$( echo ${name_sub} | sed 's/)/-/' )
    name_safe=$( esc_str "${name}" )
    submits_list+=("parsenormals_${name_sub}")

    # Check for incompletion of the minimization step
    if [ -d ${name}_minimized -a ! -f ${name}_minimized/minimized.xyz ]; then
	rm -r ${name}_minimized
    fi

    # Run the minimization and comparison
    if [ ! -d ${name}_minimized ]; then

	# Check if the copy folder or file are missing
	if [ ! -d ${name}_freq -o ! -f ${name}_freq/${name}_freq.out ]; then
	    echo "An error occured when trying to process ${name}_freq. Exiting..."
	    exit
	fi

	str=""
	str="module load gcc &> /dev/null\n"
	str+="${PYTHON_ABS} ${TAFFI_PATH}FF_functions/lammps_minimize_molecule.py '${name}_freq/geo_opt.xyz' '${FF_FOM}' -o '${name}_minimized' -lammps_exe ${LAMMPS_EXE} > /dev/null \n"	
	str+="if [ ! -f ${name_safe}_minimized/minimized.xyz ]; then\n"
	str+="    echo 'An error occured during the minimization of ${name}.xyz'\n"
	str+="    exit\n"
	str+="fi\n\n"
	str+="# Calculate the FF normal modes, and create an aligned geometry\n"
	str+="${PYTHON_ABS} ${TAFFI_PATH}Parsers/parse_normal_modes.py -geo '${name}_minimized/minimized.xyz' -FF '${FF_FOM}' -o '${name}_FF' -out '${name}_freq/${name}_freq.out' > /dev/null \n"
	str+="${PYTHON_ABS} ${TAFFI_PATH}Parsers/align_points.py '${name}_freq/geo_opt.xyz' '${name}_minimized/minimized.xyz' -o '${name}_aligned' > /dev/null  \n"
	ids_sub+=$( ${SHELL_SUBMIT} "${str}" parsenormals_${name_sub} -p 1 -t 1 -q ${FOM_GEOOPT_Q} )
    fi
done

# Wait until jobs complete
check_queue "${ids_sub[@]}"
echo "COMPLETE: Normal mode comparisons."

# Clean up run scripts
for i in ${submits_list[*]}; do    
    if [ -f ${i}.submit ]; then
	rm ${i}*
    fi
done

# Plot the normal mode comparisons
if [ ! -f nm.txt ]; then
    ${PYTHON_ABS} ${TAFFI_PATH}Parsers/plot_normal_modes.py >> nm.log
fi
cd ..

############################################################
# Run density, dielectric, and enthalpy of vap assessments #
############################################################
if [ ! -d ${FOM_name} ]; then
    mkdir ${FOM_name}
    cp *.xyz ${FOM_name}/.
fi

# Generate the single molecule and condensed phase trajectories
cd ${FOM_name}
for i in *.xyz; do

    name=$(echo "${i}" | awk -F '[.]' '{ print $1 }' )
    N_mol=$( awk '{ if ( NR == 1 ) { print $1; exit } }' ${i} )
    N_mol=$( echo "5000/${N_mol}" | bc )

    # Insert appropriate escape characters
    i=$( echo $i | sed  's/(/\\\(/g' )
    i=$( echo $i | sed  's/)/\\\)/g' )
    oname=$( echo $name | sed  's/(/\\\(/g' )
    oname=$( echo $oname | sed  's/)/\\\)/g' )

    # Optional read-temperature from file functionality. 
    # If the FOM_TEMPS file is unspecified, the default behavior falls back on 298K
    if [[ -z "${FOM_TEMPS// }" ]]; then 
	temps="298"
    else 
	temps=$( ${PYTHON_ABS} ${DIR}/get_temps.py ${FOM_TEMPS} ${name} )	    

	# If the molecule is not found in the FOM_TEMPS file, then the script falls back on 298K
	if [[ -z "${temps// }" ]]; then
	    temps="298"
	fi
    fi

    # If the *_gas folder doesn't exist, run the job generation script
    if [ ! -d "${name}_gas" ]; then
	mkdir ${name}_gas
	# Generate the five trajectories at each temperature 
	for t in ${temps[*]}; do
	    if [ ! -d "${name}_gas/${t}_K" ]; then
		mkdir ${name}_gas/${t}_K
		str=""
		for j in $( seq 1 5 ); do
		    str+="${PYTHON_ABS} ${TAFFI_PATH}FF_functions/gen_mol_traj.py ${i} '${FF_FOM}' -T ${t} -T_A ${t} -o ${oname}_gas/${t}_K/run_${j}\n"
		done
		ids_sub+=( $( ${SHELL_SUBMIT} "${str}" ${name}_gas -p 1 -t 4 -q ${FOM_GEOOPT_Q} ) )
	    fi
	done
    fi

    # If the *_cond folder doesn't exist, run the job generation script
    if [ ! -d "${name}_cond" ]; then
	mkdir ${name}_cond

	# Generate the five trajectories at each temperature 
	for t in ${temps[*]}; do
	    if [ ! -d "${name}_cond/${t}_K" ]; then
		mkdir ${name}_cond/${t}_K
		str=""
		for j in $( seq 1 5 ); do
		    str+="${PYTHON_ABS} ${TAFFI_PATH}FF_functions/gen_md_for_sampling.comet.py ${i} -FF '${FF_FOM}' -N ${N_mol} -T ${t} -T_A 10 -t_A 1E5 -t 0 -t_ext 2E5 -o ${oname}_cond/${t}_K/run_${j} -mixing wh --tail '--safe_rebuild\n"
		done
		echo -e "string:\n${str}"
		ids_sub+=( $( ${SHELL_SUBMIT} "${str}" ${name}_cond -p 1 -t 4 -q ${FOM_GEOOPT_Q} ) )
	    fi
	done
    fi
done
exit
# Wait until all of the md generation jobs have completed
check_queue "${ids_sub[@]}"
echo "COMPLETE: MD job generation."

# Run single molecule and condensation jobs
ids_sub=()    
ids_sub+=( $(${PYTHON_ABS} ${LAMMPS_SUBMIT} -f "run_*.in.init" -d "gas" -p ${FOM_MD_GAS_PROCS} -t ${FOM_MD_GAS_WT} -q ${FOM_MD_GAS_Q} -npp ${FOM_MD_GAS_NPP} -o single_runs -sched ${FOM_MD_GAS_SCHED} -ppn ${FOM_MD_GAS_PPN} -path_to_exe ${LAMMPS_EXE} -size ${FOM_MD_GAS_SIZE} --no_gpu --silent) )

ids_sub+=( $(${PYTHON_ABS} ${LAMMPS_SUBMIT} -f "run_*.in.init" -d "cond" -p ${FOM_MD_COND_PROCS} -t ${FOM_MD_COND_WT} -q ${FOM_MD_COND_Q} -npp ${FOM_MD_COND_NPP} -o cond_runs -sched ${FOM_MD_COND_SCHED} -ppn ${FOM_MD_COND_PPN} -path_to_exe ${LAMMPS_EXE} -size ${FOM_MD_COND_SIZE} --no_gpu --silent) )

# Wait until jobs complete
check_queue "${ids_sub[@]}"

# Run the md extensions
for i in $( seq 55 ); do
    ids_sub=()    
    submits_list=()
    for j in *.xyz; do

	name=$(echo "${j}" | awk -F '[.]' '{ print $1 }' )
	name_sub=$( echo ${name} | sed 's/(/-/' )
	name_sub=$( echo ${name_sub} | sed 's/)/-/' )
	name_safe=$( esc_str "${name}" )
	submits_list+=("ext_runs_${name_sub}")
	num=$(( ${i} - 1 ))

	# Optional read-temperature from file functionality. 
	# If the FOM_TEMPS file is unspecified, the default behavior falls back on 298K
	if [[ -z "${FOM_TEMPS// }" ]]; then 
	    temps="298"
	else 
	    temps=$( ${PYTHON_ABS} ${DIR}/get_temps.py ${FOM_TEMPS} ${name} )	    

	    # If the molecule is not found in the FOM_TEMPS file, then the script falls back on 298K
	    if [[ -z "${temps// }" ]]; then
		temps="298"
	    fi
	fi

	for k in $( seq 1 5 ); do	    
	    for t in ${temps[*]}; do
		if [ ! -f ${name}_cond/${t}_K/run_${k}/${num}.thermo.avg -a ! -f ${name}_cond/${t}_K/run_${k}/equil.lammpstrj ]; then
		    ids_sub+=( $(${PYTHON_ABS} ${LAMMPS_SUBMIT} -f "*extend.in.init" -d "${name}_cond ${t} run_${k}" -p ${FOM_MD_COND_PROCS} -t ${FOM_MD_COND_WT} -q ${FOM_MD_COND_Q} -npp ${FOM_MD_COND_NPP} -o ext_runs_${name_sub} -sched ${FOM_MD_COND_SCHED} -ppn ${FOM_MD_COND_PPN} -size ${FOM_MD_COND_SIZE} -path_to_exe ${LAMMPS_EXE} --no_gpu --silent --overwrite) )    
		fi
	    done
	done
    done
    check_queue "${ids_sub[@]}"
done

# Clean up run scripts
for i in ${submits_list[*]}; do    
    if [ -f ${i}.submit ]; then
	rm ${i}*
    fi
done
echo "COMPLETE: Benchmark MD simulations with final parameters."
exit
# Combine the trajectories
for i in *.xyz; do

    name=$(echo "${i}" | awk -F '[.]' '{ print $1 }' )
    cd ${i}

    # Optional read-temperature from file functionality. 
    # If the FOM_TEMPS file is unspecified, the default behavior falls back on 298K
    if [[ -z "${FOM_TEMPS// }" ]]; then 
	temps="298"
    else 
	temps=$( ${PYTHON_ABS} ${DIR}/get_temps.py ${FOM_TEMPS} ${name} )	    

	# If the molecule is not found in the FOM_TEMPS file, then the script falls back on 298K
	if [[ -z "${temps// }" ]]; then
	    temps="298"
	fi
    fi

    for t in ${temps[*]}; do
	cd ${t}
	for j in $( seq 1 5 ); do
	    cd run_${j}
	    if [ ! -f equil.lammpstrj ]; then
		${TAFFI_PATH}/Automation_Scripts/combine_ext.sh '*.sys.lammpstrj' 0 54
		mv combined.lammpstrj equil.lammpstrj
	    fi

	    if [ ! -f thermo.avg ]; then
		${TAFFI_PATH}/Automation_Scripts/combine_ext.sh '*.thermo.avg' 0 54 --thermo
		mv combined.thermo thermo.avg
	    fi 
	    cd ..
	done
	cd ..
    done
    cd ..
done
echo "COMPLETE: Collating run data."
exit

# Parse the self-diffusion coefficients from the condensed phase runs
ids_sub=()    
for i in *cond/; do
    cd "${i}"
    name=$( echo ${i} | awk -F '[/]' '{ print $1 }' )
    name_sub=$( echo ${name} | sed 's/(/-/' )
    name_sub=$( echo ${name_sub} | sed 's/)/-/' )
    name_safe=$( esc_str "${name}" )

    cp ${name}.data equil.data
    cp ${name}.map equil.map

    # Check if data from a previous dffusivity parse is present and act appropriately
    if [ ! -d diff ]; then
	ids_sub+=( $( ${SHELL_SUBMIT} "${PYTHON_ABS} ${TAFFI_PATH}Parsers/msd_parse_multi.py -name equil.lammpstrj -start 1000 -every 5 -mols 0 -T 298 -folder diff" ${name}_diff -p 1 -t 4 -q ${FOM_GEOOPT_Q} ) )
    elif [ ! -f diff/molecule_0_msd.txt ]; then
	rm -r diff
	ids_sub+=( $( ${SHELL_SUBMIT} "${PYTHON_ABS} ${TAFFI_PATH}Parsers/msd_parse_multi.py -name equil.lammpstrj -start 1000 -every 5 -mols 0 -T 298 -folder diff" ${name}_diff -p 1 -t 4 -q ${FOM_GEOOPT_Q} ) )
    fi

    cd ..
done

# Wait until jobs complete
check_queue "${ids_sub[@]}"
echo "COMPLETE: Parsed self-diffusion, density, and enthalpy of vaporization or model compound(s)."

# Print FOM summary the densities and print summary of values
kb=0.001987204118
T=298
kcal2kJ=4.184
printf "\n%40s   %20s   %40s %20s" "Molecule" "density (g/cm^3)" "diff (1E5 * cm^2/s)" "H_vap (kcal/mol)" > summary.txt
for i in *cond/; do
    name=$( echo ${i} | awk -F '_cond' '{ print $1 }' )
    density=$(${TAFFI_PATH}Parsers/avg_thermo.sh ${name}_cond/thermo.avg 1000 60000 3)
    diff=$( awk '{ if ( $1 == "The" && $2 == "Self-Diffusion" && $3 == "Coefficient" ) { printf("%12.6f",$9*1.0E5) }}' ${name}_cond/diff/msd_parse.log)
    U_vap=$( ${TAFFI_PATH}Parsers/avg_thermo.sh ${name}_gas/thermo.equil.avg 1000 0 5 )
    U_bulk=$( ${TAFFI_PATH}Parsers/avg_thermo.sh ${name}_cond/thermo.avg 1000 0 5 )
    N_mol=$( awk '{ if ( NR == 1 ) { print $1; exit } }' ${name}.xyz )
    N_mol=$( echo "1000/${N_mol}" | bc )
    U_bulk=$( echo "${U_bulk} / ${N_mol}" | bc -l )
    H_vap=$( echo "(${U_vap} - ${U_bulk} + ${kb} * ${T} ) * ${kcal2kJ}" | bc -l )
    printf "\n%40s %12.4f %40s %23.4f" ${name} ${density} ${diff} ${H_vap} >> summary.txt
done
echo "" >> summary.txt

# Return to outer folder
#cd ../..

