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
    SYS=$2
fi
echo $VAR

#SYS=('jetPtSys_nCentMix10_nPsiMix16_nVZ5_nMix100' 'jetPtSys_nCentMix40_nPsiMix16_nVZ5_nMix200' 'jetPtSys_nCentMix40_nPsiMix16_nVZ5_nMix100')
#SYS=('jetPtSys_nCentMix10_nPsiMix16_nVZ5_nMix50' 'jetPtSys_nCentMix10_nPsiMix8_nVZ5_nMix50' 'jetPtSys_nCentMix10_nPsiMix16_nVZ5_nMix10' 'jetPtSys_nCentMix20_nPsiMix8_nVZ5_nMix50' 'jetPtSys_nCentMix20_nPsiMix4_nVZ5_nMix50' 'jetPtSys_nCentMix40_nPsiMix4_nVZ5_nMix50' 'jetPtSys_nCentMix40_nPsiMix4_nVZ3_nMix50' 'jetPtSys_nCentMix40_nPsiMix4_nVZ3_nMix10')
#SYS=('jetPtSys_nCentMix40_nPsiMix4_nVZ3_nMix5' 'jetPtSys_nCentMix40_nPsiMix4_nVZ3_nMix3' 'jetPtSys_nCentMix40_nPsiMix4_nVZ3_nMix1' 'jetPtSys_nCentMix80_nPsiMix4_nVZ3_nMix10' 'jetPtSys_nCentMix80_nPsiMix4_nVZ3_nMix5')
#SYS=('jetPtSys_nCentMix10_nPsiMix16_nVZ5_nMix10' 'jetPtSys_nCentMix40_nPsiMix4_nVZ3_nMix1' 'jetPtSys_nCentMix20_nPsiMix8_nVZ5_nMix50' 'jetPtSys_nCentMix80_nPsiMix4_nVZ3_nMix5')
#SYS=('jetPtSys_nCentMix80_nPsiMix4_nVZ1_nMix5' 'jetPtSys_nCentMix80_nPsiMix4_nVZ1_nMix10' 'jetPtSys_nCentMix80_nPsiMix4_nVZ1_nMix100' 'jetPtSys_nCentMix80_nPsiMix4_nVZ0_nMix100')
#SYS=('jetPtSys_nCentMix80_nPsiMix4_nVZ1_nMix10' 'jetPtSys_nCentMix80_nPsiMix4_nVZ1_nMix100' 'jetPtSys_nCentMix80_nPsiMix16_nVZ3_nMix100' 'jetPtSys_nCentMix80_nPsiMix4_nVZ0_nMix100')
#SYS=('jetPtSys_nCentMix80_nPsiMix16_nVZ3_nMix100_phoFinerBin' 'jetPtSys_nCentMix80_nPsiMix16_nVZ3_nMix100_phoCoarserBin')
#SYS=('jetPtSys_nCentMix80_nPsiMix16_nVZ3_nMix100_phoFinerBin' 'jetPtSys_nCentMix80_nPsiMix16_nVZ3_nMix100_phoFinerBin')
#SYS=('jetPtSys_nCentMix80_nPsiMix16_nVZ3_nMix100' 'jetPtSys_nCentMix80_nPsiMix16_nVZ3_nMix20' 'jetPtSys_nCentMix80_nPsiMix16_nVZ3_nMix10' 'jetPtSys_nCentMix80_nPsiMix16_nVZ3_nMix5' 'jetPtSys_nCentMix80_nPsiMix8_nVZ3_nMix100' 'jetPtSys_nCentMix80_nPsiMix8_nVZ3_nMix20' 'jetPtSys_nCentMix80_nPsiMix8_nVZ3_nMix10' 'jetPtSys_nCentMix80_nPsiMix8_nVZ3_nMix5')
#SYS=('jetPtSys_nCentMix80_nPsiMix16_nVZ3_nMix20' 'jetPtSys_nCentMix80_nPsiMix16_nVZ3_nMix10' 'jetPtSys_nCentMix80_nPsiMix16_nVZ3_nMix5' 'jetPtSys_nCentMix80_nPsiMix16_nVZ3_nMix1')
SYS=('jetPtSys_nCentMix80_nPsiMix16_nVZ3_nMix1')
for ((i=0; i< ${#SYS[@]};i++))
do
    #time ./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_${VER}_${SYS[$i]}.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbMC_${VER}_${SYS[$i]}_${DATE}_signal.log &
    ./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_${VER}_${SYS[$i]}.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbMC_${VER}_${SYS[$i]}_${DATE}_bkg.log &
done

#
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_${VER}_nominal.config 0 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbMC_${VER}_nominal_${DATE}_signal.log & 
#./bin/phoTaggedJetRaa_jetPt.exe input/phoTagJetRaa/phoTagJetRaa_PbPbMC_${VER}_nominal.config 1 >& ./log/$VER/phoTagJetRaa_jetPt_PbPbMC_${VER}_nominal_${DATE}_bkg.log
#
wait $(jobs -p)
echo 'DONE ./bin/phoTaggedJetRaa_jetPt.exe'

mv /direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/figures/*.pdf /direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/figures/backup
for ((i=0; i< ${#SYS[@]};i++))
do
    root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/jetPt/draw_jetPtDist_withPurityEfficiencyCorrection_v2.C("PbPb", "'${VAR}'", "'${SYS[$i]}'", 1)' & 
done
