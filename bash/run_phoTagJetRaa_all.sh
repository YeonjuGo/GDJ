#!/bin/bash
export DOGLOBALDEBUGROOT=0
DATE=`date +%Y%m%d`
VER=v1
VAR=
echo "version = $VER"
echo "systematics = $VAR"
mkdir -p log/$VER/

##### photon energy
bash bash/run_phoTagJetRaa_photonEnergy.sh $VER $VAR
bash bash/run_phoTagJetRaa_photonEfficiency.sh $VER $VAR
bash bash/run_phoTagJetRaa_photonPurity.sh $VER $VAR
bash bash/run_phoTagJetRaa_jetPt.sh $VER $VAR

##### jet pT for signal photon
#time ./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbData_$VER.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbData_${VER}_${DATE}_signal.log &
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VER.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbMC_${VER}_${DATE}_signal.log &
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PPData_$VER.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PPData_${VER}_${DATE}_signal.log &
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VER.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PPMC_${VER}_${DATE}_signal.log &
#
###### jet pT for background photon
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbData_$VER.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbData_${VER}_${DATE}_bkg.log &
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VER.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbMC_${VER}_${DATE}_bkg.log &
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PPData_$VER.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PPData_${VER}_${DATE}_bkg.log &
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VER.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PPMC_${VER}_${DATE}_bkg.log &

#### jet energy for unfolding (response matrix) 
#./bin/phoTaggedJetRaa_jetEnergy.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VER.config >& ./log/$VER/phoTagJetRaa_jetEnergy_PbPbMC_${VER}_${DATE}.log &
#./bin/phoTaggedJetRaa_jetEnergy.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VER.config >& ./log/$VER/phoTagJetRaa_jetEnergy_PPMC_${VER}_${DATE}.log &

wait $(jobs -p)
