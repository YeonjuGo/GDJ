#!/bin/bash
DATE=`date +%Y%m%d`
VER=$1
VAR=$1

if [ $# -lt 1 ]; then
  echo "Usage: ./bash/run_phoTagJetRaa_jetEnergy.sh <version> <systematic>"
  exit 1
fi

if [ $# -eq 2 ]; then
    VAR=$1"_"$2
fi
echo $VAR

#### jet energy for unfolding (response matrix) 
./bin/phoTaggedJetRaa_jetEnergy.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VAR.config >& ./log/$VER/phoTagJetRaa_jetEnergy_PbPbMC_${VAR}_${DATE}.log &
./bin/phoTaggedJetRaa_jetEnergy.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VAR.config >& ./log/$VER/phoTagJetRaa_jetEnergy_PPMC_${VAR}_${DATE}.log 

echo 'DONE ./bin/phoTaggedJetRaa_jetEnergy.exe'
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetEnergy/draw_jetEnergy_responseMatrix_v2.C("PP", "'${VAR}'")' & 
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetEnergy/draw_jetEnergy_responseMatrix_v2.C("PbPb", "'${VAR}'")' 
