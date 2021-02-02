#!/bin/bash
DATE=`date +%Y%m%d`

if [ $# -lt 1 ]; then
  echo "Usage: ./bash/run_phoTagJetRaa_jetPt_1DUnfolding.sh <version> <systematic>"
  exit 1
elif [ $# -eq 1 ]; then
    SYS="nominal"
elif [ $# -eq 2 ]; then
    SYS=$2
fi

VER=$1
VAR=$VER"_"$SYS
echo $VAR

##### jet pT for signal photon
#time ./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbData_$VAR.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbData_${VAR}_${DATE}_signal.log &
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VAR.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbMC_${VAR}_${DATE}_signal.log &
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PPData_$VAR.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PPData_${VAR}_${DATE}_signal.log &
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VAR.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PPMC_${VAR}_${DATE}_signal.log &
#
###### jet pT for background photon
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbData_$VAR.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbData_${VAR}_${DATE}_bkg.log &
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VAR.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbMC_${VAR}_${DATE}_bkg.log &
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PPData_$VAR.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PPData_${VAR}_${DATE}_bkg.log &
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VAR.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PPMC_${VAR}_${DATE}_bkg.log
#wait $(jobs -p)

echo 'DONE ./bin/phoTaggedJetRaa_jetPt_1DUnfolding.exe'
mv /direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/figures/*1DUnfolding\*.pdf /direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/figures/backup_1DUnfolding
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityEfficiencyCorrection_v2.C("PP", "'${VER}'", "'${SYS}'", 0)' &
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityEfficiencyCorrection_v2.C("PbPb", "'${VER}'", "'${SYS}'", 0)' &
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityEfficiencyCorrection_v2.C("PP", "'${VER}'", "'${SYS}'", 1)' &
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityEfficiencyCorrection_v2.C("PbPb", "'${VER}'", "'${SYS}'", 1)'

#root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityCorrection_v2.C("PP", "'${VER}'", "'${SYS}'", 0)' &
#root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityCorrection_v2.C("PbPb", "'${VER}'", "'${SYS}'", 0)' &
#root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityCorrection_v2.C("PP", "'${VER}'", "'${SYS}'", 1)' &
#root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityCorrection_v2.C("PbPb", "'${VER}'", "'${SYS}'", 1)'

#root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityEfficiencyCorrection_v2.C("PP", "'${VER}'", "'${SYS}'", 0)' &
#root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityEfficiencyCorrection_v2.C("PbPb", "'${VER}'", "'${SYS}'", 0)' &
#root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityEfficiencyCorrection_v2.C("PP", "'${VER}'", "'${SYS}'", 1)' &
#root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityEfficiencyCorrection_v2.C("PbPb", "'${VER}'", "'${SYS}'", 1)'
