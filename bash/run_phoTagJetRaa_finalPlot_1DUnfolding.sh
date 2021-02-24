#!/bin/bash
DATE=`date +%Y%m%d`

if [ $# -lt 1 ]; then
  echo "Usage: ./bash/run_phoTagJetRaa_finalPlot_1DUnfolding.sh <version> <systematic>"
  exit 1
elif [ $# -eq 1 ]; then
    SYS="nominal"
elif [ $# -eq 2 ]; then
    SYS=$2
fi

VER=$1
VAR=$VER"_"$SYS
echo $VAR

#method=('bin')
method=('bayes' 'svd' 'inv' 'bin')
#method="bayes"
END=5

echo "Final plot! for" $VAR

for ((i=0; i< ${#method[@]};i++))
do
    if [ ${method[$i]} == 'bayes' ] || [ ${method[$i]} == 'svd' ]; then
        for j in $(seq 1 $END)
        do
            root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/finalPlot/draw_finalPlots_for1DUnfolding_v1.C("'${VER}'", "'${SYS}'", "'${method[$i]}'", '$j')' &
            root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/finalPlot/draw_finalPlots_for1DUnfolding_v1.C("'${VER}'", "'${SYS}'", "'${method[$i]}'", '$j')' &
        done
    else
        root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/finalPlot/draw_finalPlots_for1DUnfolding_v1.C("'${VER}'", "'${SYS}'", "'${method[$i]}'", 0)' &
        root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/finalPlot/draw_finalPlots_for1DUnfolding_v1.C("'${VER}'", "'${SYS}'", "'${method[$i]}'", 0)' &
    fi
done
