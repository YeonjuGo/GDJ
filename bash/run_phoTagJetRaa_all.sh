#!/bin/bash
export DOGLOBALDEBUGROOT=0
DATE=`date +%Y%m%d`
VER=v1
VAR=
echo "version = $VER"
echo "systematics = $VAR"
mkdir -p log/$VER/

##### photon energy
#bash bash/run_phoTagJetRaa_photonEnergy.sh $VER $VAR &
#bash bash/run_phoTagJetRaa_photonEfficiency.sh $VER $VAR &
#bash bash/run_phoTagJetRaa_photonPurity.sh $VER $VAR &
#bash bash/run_phoTagJetRaa_jetPt.sh $VER $VAR &
bash bash/run_phoTagJetRaa_jetEnergy.sh $VER $VAR 
wait $(jobs -p)
bash bash/run_unfolding.sh $VER $VAR 

