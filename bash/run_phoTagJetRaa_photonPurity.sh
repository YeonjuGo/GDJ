#!/bin/bash
DATE=`date +%Y%m%d`

if [ $# -lt 1 ]; then
  echo "Usage: ./bash/run_phoTagJetRaa_photonPurity.sh <version> <systematic>"
  exit 1
elif [ $# -eq 1 ]; then
    SYS="nominal"
elif [ $# -eq 2 ]; then
    SYS=$2
fi

VER=$1
VAR=$VER"_"$SYS
echo "./bin/phoTaggedJetRaa_photonPurity.exe for" $VAR

###### photon purity
./bin/phoTaggedJetRaa_photonPurity.exe input/phoTagJetRaa/phoTagJetRaa_PbPbData_$VAR.config >& ./log/$VER/phoTagJetRaa_photonPurity_PbPbData_${VAR}_${DATE}.log &
./bin/phoTaggedJetRaa_photonPurity.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VAR.config >& ./log/$VER/phoTagJetRaa_photonPurity_PbPbMC_${VAR}_${DATE}.log &
./bin/phoTaggedJetRaa_photonPurity.exe input/phoTagJetRaa/phoTagJetRaa_PPData_$VAR.config >& ./log/$VER/phoTagJetRaa_photonPurity_PPData_${VAR}_${DATE}.log &
./bin/phoTaggedJetRaa_photonPurity.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VAR.config >& ./log/$VER/phoTagJetRaa_photonPurity_PPMC_${VAR}_${DATE}.log 
wait $(jobs -p)

echo 'DONE ./bin/phoTaggedJetRaa_photonPurity.exe'
mv /direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonPurity/figures/*.pdf /direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonPurity/figures/backup

root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonPurity/draw_photon_purity_v3.C("PP", "'${VER}'", "'${SYS}'")' &
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonPurity/draw_photon_purity_v3.C("PbPb", "'${VER}'", "'${SYS}'")'
