#!/bin/bash
DATE=`date +%Y%m%d`
VER=$1
VAR=$1

if [ $# -lt 1 ]; then
  echo "Usage: ./bash/run_phoTagJetRaa_jetPt.sh <version> <systematic>"
  exit 1
fi

if [ $# -eq 2 ]; then
    VAR=$1"_"$2
fi
echo $VAR

##### jet pT for signal photon
./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VAR.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbMC_${VAR}_${DATE}_signal.log &
./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VAR.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PPMC_${VAR}_${DATE}_signal.log &

##### jet pT for background photon
./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_$VAR.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbMC_${VAR}_${DATE}_bkg.log &
./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_$VAR.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PPMC_${VAR}_${DATE}_bkg.log

echo 'DONE ./bin/phoTaggedJetRaa_jetPt.exe'
mv /direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/figures/*.pdf /direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/figures/backup
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityEfficiencyCorrection_v2.C("PP", "'${VAR}'", 1)' &
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityEfficiencyCorrection_v2.C("PbPb", "'${VAR}'", 1)'
