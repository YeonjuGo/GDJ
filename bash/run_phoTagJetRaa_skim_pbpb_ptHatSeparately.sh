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
./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PbPbMC_35to50.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PbPbMC_${DATE}_35to50.log
./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PbPbMC_50to70.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PbPbMC_${DATE}_50to70.log
./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PbPbMC_70to140.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PbPbMC_${DATE}_70to140.log
./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PbPbMC_140to280.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PbPbMC_${DATE}_140to280.log

#./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PbPbMC.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PbPbMC_${DATE}.log
#./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PPMC.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PPMC_${DATE}.log

##### data 
#./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PbPbData.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PbPbData_${DATE}.log
#./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PPData.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PPData_${DATE}.log

##### PYTHIA+OVERLAY for TEST! 
#./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PbPbMC_test.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PbPbMC_${DATE}.log
#./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe input/ntuplePreProc/ntuplePreProc_phoTaggedJetRaa_PbPbData_test.config >& ./log/skim/ntuplePreProc_phoTaggedJetRaa_PbPbData_${DATE}.log
