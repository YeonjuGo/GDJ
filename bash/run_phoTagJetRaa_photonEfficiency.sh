#!/bin/bash
DATE=`date +%Y%m%d`
VER=$1
VAR=$1

if [ $# -eq 2 ]; then
    VAR=$1"_"$2
fi
echo $VAR

##### photon efficiency
./bin/phoTaggedJetRaa_photonEfficiency.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VAR.config >& ./log/$VER/phoTagJetRaa_photonEfficiency_PbPbMC_${VAR}_${DATE}.log &
./bin/phoTaggedJetRaa_photonEfficiency.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VAR.config >& ./log/$VER/phoTagJetRaa_photonEfficiency_PPMC_${VAR}_${DATE}.log 
echo 'DONE ./bin/phoTaggedJetRaa_photonEfficiency.exe'
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonEfficiency/draw_photon_efficiency_v4.C("PP", "'${VAR}'")' &
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonEfficiency/draw_photon_efficiency_v4.C("PbPb", "'${VAR}'")'
