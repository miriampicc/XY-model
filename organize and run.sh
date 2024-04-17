#!/bin/bash

############### SET MAIN FOLDERS ###############
BASEDIR=${HOME}/XY-model
SCRIPT_DIR=${BASEDIR}/launch_script

#############################################

time_limit="1-00:00:00"

LLIST="8 12 16 20 24 32 40 48"

################ Input Parameters for the Monte Carlo simulation #################

nsteps=5000
transient=1000
tau=32
T=0.558
restart=0
K=5.0
J1=1.0
J2=2.0
e=0.1
theta_box=0.78539816339
charged_fluid=1 #if we have charged_fluid=1 then we will have a superconductor, if charged_fluid=0 then we have a superfluid, it would be good if we have the option of both


EXECUTE_DIR="../build/release-conan"

################################################################################3

for L in $LLIST; do

############Creation of the output folder and of the two files of initialization####################

    cd ${BASEDIR}/Output_TBG


    if [ ! -d ./SK_${K} ]; then

    mkdir -p K_${K}
    fi

    cd K_${K}

    if [ ! -d ./Se_${e} ]; then
    mkdir -p e_${e}
    fi

    cd e_${e}

    if [ ! -d ./SL${L}_K${K}_e${e} ]; then
    mkdir -p L${L}_K${K}_e${e}
    fi

    cd L${L}_K${K}_e${e}

    if [ ! -d ./ST_${T} ]; then
    mkdir -p T_${T}
    fi

    DIR_OUT=${BASEDIR}/Output_TBG/K_${K}/e_${e}/L${L}_K${K}_e${e}/T_${T}
    DIR_IN="${DIR_OUT}"


    #################Creation of the submit_runs script#########################

    jobname="T_${T}_L${L}_K${K}_e${e}"
    nnodes=1
    ntasks=1 #parallel tempering over ntasks temperatures


    cd ${SCRIPT_DIR}
    DIR_PAR="${DIR_OUT}"

    #SEED= If I want to repeat exactly a simulation I could initialize the random number generator exactly at the same way

    EXECUTE_DIR="../build/release-conan"

    #SBATCH --nodes=${nnodes}               # Number of nodes

    echo "#!/bin/bash
    #SBATCH --job-name=${jobname}          # Name of the job
    #SBATCH --time=${time_limit}               # Allocation time
    #SBATCH --mem-per-cpu=2000              # Memory per allocated cpu
    #SBATCH --nodes=${nnodes}               # Number of nodes
    #SBATCH --ntasks=${ntasks}
    #SBATCH --output=${DIR_PAR}/logs/log_${jobname}.o
    #SBATCH --error=${DIR_PAR}/logs/log_${jobname}.e

    srun ${EXECUTE_DIR}/CMT ${L} ${nsteps} ${transient} ${tau} ${T} ${restart} ${K} ${J1} ${J2} ${theta_box} ${e} ${DIR_IN} ${DIR_OUT}

    " >  submit_run

    #Submission of the work --> sbatch submit_runs

    mkdir -p ${DIR_PAR}/logs

    sbatch submit_run

done
