#!/bin/bash
DATE=`date +%Y%m%d`
VER=$1
VAR=$1

if [ $# -eq 2 ]; then
    VAR=$1"_"$2
fi
echo $VAR

##### photon energy
./bin/phoTaggedJetRaa_photonEnergy.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VAR.config >& ./log/$VER/phoTagJetRaa_photonEnergy_PbPbMC_${VAR}_${DATE}.log &
./bin/phoTaggedJetRaa_photonEnergy.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VAR.config >& ./log/$VER/phoTagJetRaa_photonEnergy_PPMC_${VAR}_${DATE}.log 
echo 'DONE ./bin/phoTaggedJetRaa_photonEnergy.exe'
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonEnergy/draw_photon_energy_v1.C("PP", "'${VAR}'")' &
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonEnergy/draw_photon_energy_v1.C("PbPb", "'${VAR}'")'
