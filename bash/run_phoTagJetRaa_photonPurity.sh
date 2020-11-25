#!/bin/bash
DATE=`date +%Y%m%d`
VER=$1
VAR=$1

if [ $# -eq 2 ]; then
    VAR=$1"_"$2
fi
echo $VAR

##### photon purity
./bin/phoTaggedJetRaa_photonPurity.exe input/phoTagJetRaa/phoTagJetRaa_PbPbData_$VAR.config >& ./log/$VER/phoTagJetRaa_photonPurity_PbPbData_${VAR}_${DATE}.log &
./bin/phoTaggedJetRaa_photonPurity.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VAR.config >& ./log/$VER/phoTagJetRaa_photonPurity_PbPbMC_${VAR}_${DATE}.log &
./bin/phoTaggedJetRaa_photonPurity.exe input/phoTagJetRaa/phoTagJetRaa_PPData_$VAR.config >& ./log/$VER/phoTagJetRaa_photonPurity_PPData_${VAR}_${DATE}.log &
./bin/phoTaggedJetRaa_photonPurity.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VAR.config >& ./log/$VER/phoTagJetRaa_photonPurity_PPMC_${VAR}_${DATE}.log &
echo 'DONE ./bin/phoTaggedJetRaa_photonPurity.exe'
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonPurity/draw_photon_purity_v3.C("PP", "'${VAR}'")' &
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonPurity/draw_photon_purity_v3.C("PbPb", "'${VAR}'")'
