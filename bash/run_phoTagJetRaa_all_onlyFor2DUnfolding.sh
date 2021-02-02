#!/bin/bash
export DOGLOBALDEBUGROOT=0
DATE=`date +%Y%m%d`
VER=v1
SYS=nominal
echo "version = $VER"
echo "systematics = $SYS"
mkdir -p log/$VER/

#bash bash/run_phoTagJetRaa_photonEnergy.sh $VER $SYS &
#bash bash/run_phoTagJetRaa_photonEfficiency.sh $VER $SYS &
#bash bash/run_phoTagJetRaa_photonPurity.sh $VER $VAR &
bash bash/run_phoTagJetRaa_jetEnergy_2DUnfolding.sh $VER $VAR &
bash bash/run_phoTagJetRaa_jetPt_2DUnfolding.sh $VER $VAR 
wait $(jobs -p)
bash bash/run_phoTagJetRaa_unfolding_2DUnfolding.sh $VER $VAR &
wait $(jobs -p)
bash bash/run_phoTagJetRaa_finalPlot_2DUnfolding.sh $VER $VAR &

