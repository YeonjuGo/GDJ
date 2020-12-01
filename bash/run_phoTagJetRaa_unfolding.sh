#!/bin/bash
DATE=`date +%Y%m%d`
VER=$1
VAR=$1

if [ $# -lt 1 ]; then
  echo "Usage: ./bash/run_phoTagJetRaa_unfolding.sh <version> <systematic>"
  exit 1
fi

if [ $# -eq 2 ]; then
    VAR=$2
#VAR=$1"_"$2
fi
echo $VAR

for j in {1..30}
do
    python /usatlas/u/goyeonju/phoTaggedJetRaa/unfolding/unfold_jetPt_v2.py PbPb bayes $j $VER $VAR &
    python /usatlas/u/goyeonju/phoTaggedJetRaa/unfolding/unfold_jetPt_v2.py PP bayes $j $VER $VAR &
done

echo 'DONE unfolding'
echo 'DONE unfolding'
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/unfolding/drawUnfoldingFigures_v2.C("PP", "'${VAR}'")' &
root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/unfolding/drawUnfoldingFigures_v2.C("PbPb", "'${VAR}'")' &
