#!/bin/bash
DATE=`date +%Y%m%d`

if [ $# -lt 1 ]; then
  echo "Usage: ./bash/run_phoTagJetRaa_jetEnergy_2DUnfolding.sh <version> <systematic>"
  exit 1
elif [ $# -eq 1 ]; then
    SYS="nominal"
elif [ $# -eq 2 ]; then
    SYS=$2
fi

VER=$1
VAR=$VER"_"$SYS
echo $VAR

echo "./bin/phoTaggedJetRaa_jetEnergy_2DUnfolding.exe for" $VAR

#### jet energy for unfolding (response matrix) 
./bin/phoTaggedJetRaa_jetEnergy_2DUnfolding.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VAR.config >& ./log/$VER/phoTagJetRaa_jetEnergy_2DUnfolding_PbPbMC_${VAR}_${DATE}.log &
./bin/phoTaggedJetRaa_jetEnergy_2DUnfolding.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VAR.config >& ./log/$VER/phoTagJetRaa_jetEnergy_2DUnfolding_PPMC_${VAR}_${DATE}.log 
wait $(jobs -p)

echo 'DONE ./bin/phoTaggedJetRaa_jetEnergy_2DUnfolding.exe'

#mv /direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetEnergy/figures/*.pdf /direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetEnergy/figures/backup
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetEnergy/draw_jetEnergy_responseMatrix_2DUnfolding_v2.C("PP", "'${VER}'", "'${SYS}'")' & 
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetEnergy/draw_jetEnergy_responseMatrix_2DUnfolding_v2.C("PbPb", "'${VER}'", "'${SYS}'")' 
