#!/bin/bash
DATE=`date +%Y%m%d`
VER=$1
VAR=$1

if [ $# -eq 2 ]; then
    VAR=$1"_"$2
fi
echo $VAR

##### jet pT for signal photon
time ./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbData_$VAR.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbData_${VAR}_${DATE}_signal.log &
./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VAR.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbMC_${VAR}_${DATE}_signal.log &
./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PPData_$VAR.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PPData_${VAR}_${DATE}_signal.log &
./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VAR.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PPMC_${VAR}_${DATE}_signal.log &

##### jet pT for background photon
./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbData_$VAR.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbData_${VAR}_${DATE}_bkg.log &
./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VAR.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbMC_${VAR}_${DATE}_bkg.log &
./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PPData_$VAR.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PPData_${VAR}_${DATE}_bkg.log &
./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VAR.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PPMC_${VAR}_${DATE}_bkg.log &

echo 'DONE ./bin/phoTaggedJetRaa_jetPt.exe'
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityEfficiencyCorrection_v2.C("PP", "'${VAR}'", 0)' &
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityEfficiencyCorrection_v2.C("PbPb", "'${VAR}'", 0)' &
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityEfficiencyCorrection_v2.C("PP", "'${VAR}'", 1)' &
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityEfficiencyCorrection_v2.C("PbPb", "'${VAR}'", 1)'
