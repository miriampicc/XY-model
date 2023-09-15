#!/bin/bash

BASEDIR=${HOME}/Desktop/hihi/KTH/XY\-model

#input parameters
L=8
nsteps=10000
restart=0

BASEDIR_OUT=${HOME}/Desktop/hihi/KTH/XY\-model/Output

T_LIST="1.5 1.4 1.3 1.2 1.1 1.0 0.95 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.5 0.4 0.2 0.1"

EXCUTE_DIR=${BASEDIR}/build/release\-conan\-gcc\-13
DIR_IN=${BASEDIR_OUT}

for T in ${T_LIST}; do

    cd ${BASEDIR_OUT}

    if [ ! -d ./ST_${T} ]; then
         mkdir -p T_${T}
    fi

    DIR_OUT=${BASEDR_OUT}/T_${T}

    cd ${EXCUTE_DIR}

    ./CMT ${L} ${nsteps} ${T} ${restart} ${DIR_IN} ${DIR_OUT}

    restart=1
    DIR_IN=${DIR_OUT}

done

