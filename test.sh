#!/bin/bash

set -e

### SETTINGS

FDSIZE=31
EP=2.7
ND=6400
NODESET="md079.06400"
#ND=655362
#NODESET="icos655362.mat"

NODESET2="${ND}_${FDSIZE}.mat"
GALEW="par = struct('test', 'galew', 'n', ${ND}, 'fd', ${FDSIZE}, 'ep', ${EP}, 'order', 4, 'gamma_c', -0.05);"
TC5="par = struct('test', 'tc5',   'n', ${ND}, 'fd', ${FDSIZE}, 'ep', ${EP}, 'order', 4, 'gamma_c', -0.05);"
GALEWN="galew-${ND}-${FDSIZE}-ep${EP}-o4-gc-0.05"
TC5N="tc5-${ND}-${FDSIZE}-ep${EP}-o4-gc-0.05"
DT="600"

echo "${GALEW} ${GALEWN}"
echo "${TC5} ${TC5N}"

# setup
mkdir -p nodesets
mkdir -p tree
mkdir -p cartfd
mkdir -p data
mkdir -p resdata

TESTS=$1

if [ "x${TESTS}x" == "xx" ]; then
  TESTS="01234"
fi

### 0: COMPILE

if [[ ${TESTS} == *0* ]]; then
  make bin/sw2.prod
fi


### 1: PREPROCESS NODES

MLCODE=""
if [[ ${TESTS} == *1* ]]; then
  rm -f "nodesets/${NODESET2}"
  rm -f "tree/${NODESET2}"
  MLCODE="${MLCODE} mt_preprocess('../orgnodesets/${NODESET}', ${FDSIZE});"
fi


### 2: GENERATE DATA

if [[ ${TESTS} == *2* ]]; then
  rm -rf "data/${TC5N}"
  rm -rf "data/${GALEWN}"
  MLCODE="${MLCODE} clear; ${GALEW} mt_save(par);"
  MLCODE="${MLCODE} clear; ${TC5} mt_save(par);"
fi

if [ -n "${MLCODE}" ]; then
  MLCODE="${MLCODE} exit;"
  ( cd matlab ; matlab -nojvm -r "${MLCODE}" )
fi


### 3: RUN

if [[ ${TESTS} == *3* ]]; then

  rm -rf "resdata/${GALEWN}-dt${DT}"
  rm -rf "resdata/${TC5}-dt${DT}"

  echo ./run ${GALEWN} 400 ${DT} 518400
  ./run ${GALEWN} 400 ${DT} 518400
  echo "size = 16 chunk_size= 400 num_chunks= 16 cpus= 4 time= 5399834449 dbg=[ -27.9535 -87.1889 -302.851 6.27435e+08 -1609.46 -856.794 -399.075 88968.8 2139.55 5001.25 310.509 6.27435e+08 ] <-- EXPECTED"

  ./run ${TC5N} 400 ${DT} 1296000
  echo "size = 16 chunk_size= 400 num_chunks= 16 cpus= 4 time= 8641111484 dbg=[ 203.336 134.071 -124.541 -1.75073e+07 -1508.3 7.12457 -296.373 -1.75073e+07 869.626 2906.25 173.942 -7895.51 ] <-- EXPECTED"
fi


### 4: PLOT RESULTS

if [[ ${TESTS} == *4* ]]; then
  MLCODE="clear; ${GALEW} saveplot(par, ${DT}); "
  MLCODE="${MLCODE}; clear; ${TC5} saveplot(par, ${DT}); "
  MLCODE="${MLCODE}; exit;"
  ( cd matlab ; matlab -nosplash -nodesktop -r "${MLCODE}" )
fi
