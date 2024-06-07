#!/bin/bash

############### SET MAIN FOLDERS ######
BASEDIR=${HOME}/XY-model/relative_density_fluctuations
#BASEDIR=/Users/mirimi/Desktop/hihi/KTH/XY-model
SCRIPT_DIR=${BASEDIR}/launch_script

#############################################

time_limit="5-0:00:00"

LLIST="8 12 20 24 32 40 48 64 96"

################ Input Parameters for the Monte Carlo simulation #################

#Hamiltonian input parameters

#Monte Carlo parameters
nsteps=500000
transient=100000
tau=32
T=0.3
restart=0
K=5.0
J1=0.0
J2=0.0
e=0.1
beta_high=1.9    #1.754   #1.818  #T=0.57 #0.55
beta_low=1.6   #T=0.59
theta_box=0.78539816339
theta_box_A=0.1

################################################################################3

for L in $LLIST; do

############Creation of the output folder and of the two files of initialization####################

    cd ${BASEDIR}/Output_TBG_rdf || exit


    if [ ! -d ./SK_${K} ]; then

    mkdir -p K_${K}
    fi

    cd K_${K} || exit

    if [ ! -d ./Se_${e} ]; then
    mkdir -p e_${e}
    fi

    cd e_${e} || exit

    if [ ! -d ./SL${L}_K${K}_e${e}_bmin${beta_low}_bmax${beta_high} ]; then
    mkdir -p L${L}_K${K}_e${e}_bmin${beta_low}_bmax${beta_high}
    fi

    DIR_OUT=${BASEDIR}/Output_TBG_rdf/K_${K}/e_${e}/L${L}_K${K}_e${e}_bmin${beta_low}_bmax${beta_high}

    #################Creation of the submit_runs script#########################

    jobname="L${L}_K${K}_e${e}_bmin${beta_low}_bmax${beta_high}_rdf"
    nnodes=2
    ntasks=64 #parallel tempering over ntasks temperatures

    #I create ntasks folder: one for each rank.

    cd ${DIR_OUT} || exit

    for ((rank=0; rank<${ntasks}; rank++)); do

      if [ ! -d ./Sbeta_${rank} ]; then
        mkdir -p beta_${rank}
      fi

    done

    cd ${SCRIPT_DIR} || exit
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

srun ${EXECUTE_DIR}/CMT ${L} ${nsteps} ${transient} ${tau} ${T} ${restart} ${K} ${J1} ${J2} ${e} ${beta_high} ${beta_low} ${theta_box} ${theta_box_A} ${DIR_OUT} &> ${DIR_PAR}/logs/log_${jobname}.o

" >  submit_run

    #Submission of the work --> sbatch submit_runs

    mkdir -p ${DIR_PAR}/logs

    sbatch submit_run

done
