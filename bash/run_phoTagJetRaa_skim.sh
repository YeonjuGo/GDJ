#!/bin/bash
export DOGLOBALDEBUGROOT=0
DATE=`date +%Y%m%d`

#if [ $# -lt 1 ]; then
#  echo "Usage: ./bash/run_phoTagJetRaa_skim.sh"
#  exit 1
#elif [ $# -eq 1 ]; then
#    SYS="nominal"
#elif [ $# -eq 2 ]; then
#    SYS=$2
#fi
#VER=$1
#VAR=$VER"_"$SYS
echo "./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe"

##### MC small stats 
./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PbPbMC.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PbPbMC_${DATE}.log
./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PPMC.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PPMC_${DATE}.log

##### data 
./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PbPbData.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PbPbData_${DATE}.log
./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PPData.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PPData_${DATE}.log

##### PYTHIA+OVERLAY for TEST! 
#./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PbPbMC_test.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PbPbMC_${DATE}.log
#./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PbPbData_test.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PbPbData_${DATE}.log
