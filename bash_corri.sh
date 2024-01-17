#!/bin/bash

#to run it: ./bash_corri.sh
#to have the permission to run it: chmod +x bash_corri.sh

BASEDIR="${HOME}/Desktop/hihi/KTH/XY-model"

#input parameters
L=12
nsteps=100000
restart=0
K=1.0
J1=1.0
J2=1.0

BASEDIR_OUT="${HOME}/Desktop/hihi/KTH/XY-model/Output"
DIR_IN=${BASEDIR_OUT}

EXCUTE_DIR=${BASEDIR}/build/release-conan-apple-clang

T_LIST=" 2.5 2.0 1.7 1.5 1.45 1.4 1.3 1.27 1.25 1.22 1.2 1.17 1.15 1.13 1.1 1.09  1.07 1.06
   1.05 1.04 1.02 1.0 0.98 0.96 0.95 0.94 0.92 0.9 0.87 0.85 0.82 0.8 0.77 0.75 0.72 0.7 0.67 0.65 0.62 0.6 0.57
   0.55 0.52 0.5 0.45 0.4 0.37 0.35 0.32 0.3 0.27
      0.25 0.22 0.2 0.17 0.15 0.14 0.13 0.11 0.1 0.09 0.08 0.07 0.05 0.03 0.02 0.01 0.009"


#T_LIST="5.0 4.5 4.0 3.5 3.0 2.5 2.4 2.3 2.2 2.1 2.0 1.9 1.8 1.7 1.6 1.5 1.45
 #       1.4 1.37 1.35 1.33 1.3 1.29 1.27 1.25 1.23 1.2 1.18 1.17 1.15 1.14 1.12 1.1 1.09 1.07
  #       1.05 1.04 1.02 1.0 0.98 0.96 0.95 0.94 0.92 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.55 0.5 0.45 0.4 0.35 0.3
   #     0.25 0.2 0.15 0.14 0.13 0.11 0.1 0.09 0.08 0.07 0.05 0.03 0.02 0.01 0.009"


#T_LIST="0.95 0.94 0.92 0.9 0.85"
#T_LIST="0.9 0.35 0.01"

        #10.0 9.0 8.0 7.0 6.0

#T_LIST=" 0.65 0.6 0.5  0.4  0.3"
#T_LIST=" 0.25 0.2 0.15 0.13 0.1"
#T_LIST=" 2.2 1.4"

for T in ${T_LIST}; do

    if [ ! -d ${BASEDIR_OUT}/ST_${T} ]; then
         mkdir -p ${BASEDIR_OUT}/T_${T}
    fi

    DIR_OUT="${BASEDIR_OUT}/T_${T}"

    ${EXCUTE_DIR}/CMT ${L} ${nsteps} ${T} ${restart} ${K} ${J1} ${J2} ${DIR_IN} ${DIR_OUT}

    restart=1
    DIR_IN=${DIR_OUT}

done

