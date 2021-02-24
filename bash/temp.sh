#!/bin/bash
export DOGLOBALDEBUGROOT=0
DATE=`date +%Y%m%d`

if [ $# -lt 1 ]; then
  echo "Usage: ./bash/run_phoTagJetRaa_all_onlyFor2DUnfolding.sh <version> <systematic>"
  exit 1
elif [ $# -eq 1 ]; then
    SYS="nominal"
elif [ $# -eq 2 ]; then
    SYS=$2
fi

VER=$1
VAR=$VER"_"$SYS
echo "./bash/run_phoTagJetRaa_all_onlyFor2DUnfolding.sh for" $VAR

mkdir -p output/$VER/
mkdir -p log/$VER/

#bash bash/run_phoTagJetRaa_photonEnergy.sh $VER $SYS &
#bash bash/run_phoTagJetRaa_photonEfficiency.sh $VER $SYS &
#bash bash/run_phoTagJetRaa_photonPurity.sh $VER $SYS &
#bash bash/run_phoTagJetRaa_jetPt_2DUnfolding.sh $VER $SYS 
#wait $(jobs -p)
bash bash/run_phoTagJetRaa_jetEnergy_2DUnfolding.sh $VER $SYS &
wait $(jobs -p)
bash bash/run_phoTagJetRaa_unfolding_2DUnfolding.sh $VER $SYS &
wait $(jobs -p)
bash bash/run_phoTagJetRaa_finalPlot_2DUnfolding.sh $VER $SYS 
