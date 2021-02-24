#!/bin/bash
export DOGLOBALDEBUGROOT=0
DATE=`date +%Y%m%d`

if [ $# -lt 1 ]; then
  echo "Usage: ./bash/run_phoTagJetRaa_photonEnergy.sh <version> <systematic>"
  exit 1
elif [ $# -eq 1 ]; then
    SYS="nominal"
elif [ $# -eq 2 ]; then
    SYS=$2
fi

VER=$1
VAR=$VER"_"$SYS
echo "./bin/phoTaggedJetRaa_photonEnergy.exe for" $VAR

##### photon energy
./bin/phoTaggedJetRaa_photonEnergy.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VAR.config >& ./log/$VER/phoTagJetRaa_photonEnergy_PbPbMC_${VAR}_${DATE}.log &
./bin/phoTaggedJetRaa_photonEnergy.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VAR.config >& ./log/$VER/phoTagJetRaa_photonEnergy_PPMC_${VAR}_${DATE}.log 
wait $(jobs -p)
echo 'DONE ./bin/phoTaggedJetRaa_photonEnergy.exe'

mv /direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonEnergy/figures/*.pdf /direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonEnergy/figures/backup
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonEnergy/draw_photon_energy_v1.C("PP", "'${VER}'", "'${SYS}'")' &
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonEnergy/draw_photon_energy_v1.C("PbPb", "'${VER}'", "'${SYS}'")'
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonEnergy/draw_photon_energy_pp_pbpb_allTogether_v1.C("'${VER}'", "'${SYS}'")'
