#!/bin/bash

#to run it: ./bash_corri.sh
#to have the permission to run it: chmod +x bash_corri.sh

BASEDIR="${HOME}/Desktop/hihi/KTH/XY-model"

#input parameters
L=16
nsteps=10000
#nsteps=100000 this is the right amount
#nsteps=1000000
restart=0

BASEDIR_OUT="${HOME}/Desktop/hihi/KTH/XY-model/Output"
DIR_IN=${BASEDIR_OUT}

EXCUTE_DIR=${BASEDIR}/build/release-conan-gcc-13

T_LIST="5.0 4.5 4.0 3.5 3.0 2.5 2.4 2.3 2.2 2.1 2.0 1.9 1.8 1.7 1.6 1.5 1.45
        1.4 1.35 1.3 1.25 1.2 1.15 1.1 1.05 1.0 0.95 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.55 0.5 0.45 0.4 0.35 0.3
        0.25 0.2 0.15 0.14 0.13 0.11 0.1 0.09 0.08 0.07 0.05 0.03 0.02 0.01 0.009"

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

    ${EXCUTE_DIR}/CMT ${L} ${nsteps} ${T} ${restart} ${DIR_IN} ${DIR_OUT}

    restart=1
    DIR_IN=${DIR_OUT}

done

