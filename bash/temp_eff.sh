#!/bin/bash
DATE=`date +%Y%m%d`

if [ $# -lt 1 ]; then
  echo "Usage: ./bash/run_phoTagJetRaa_photonEfficiency.sh <version> <systematic>"
  exit 1
elif [ $# -eq 1 ]; then
    SYS="nominal"
elif [ $# -eq 2 ]; then
    SYS=$2
fi

VER=$1
VAR=$VER"_"$SYS
echo "./bin/phoTaggedJetRaa_photonEfficiency.exe for" $VAR

SYS=('jetPtSys_nCentMix80_nPsiMix16_nVZ3_nMix50' 'jetPtSys_nCentMix80_nPsiMix8_nVZ3_nMix50')
for ((i=0; i< ${#SYS[@]};i++))
do

##### photon efficiency
./bin/phoTaggedJetRaa_photonEff.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_${VER}_${SYS}.config >& ./log/$VER/phoTagJetRaa_photonEfficiency_PbPbMC_${VER}_${SYS}_${DATE}.log &
./bin/phoTaggedJetRaa_photonEff.exe input/phoTagJetRaa/phoTagJetRaa_PPMC_${VER}_${SYS}.config >& ./log/$VER/phoTagJetRaa_photonEfficiency_PPMC_${VER}_${SYS}_${DATE}.log 

wait $(jobs -p)

echo 'DONE ./bin/phoTaggedJetRaa_photonEfficiency.exe'

mv /direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonEfficiency/figures/*.pdf /direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonEfficiency/figures/backup
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonEfficiency/draw_photon_efficiency_v4.C("PP", "'${VER}'", "'${SYS}'")' &
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/photonEfficiency/draw_photon_efficiency_v4.C("PbPb", "'${VER}'", "'${SYS}'")'

