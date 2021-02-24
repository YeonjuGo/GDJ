#!/bin/bash
DATE=`date +%Y%m%d`

if [ $# -lt 1 ]; then
  echo "Usage: ./bash/run_phoTagJetRaa_unfolding_2DUnfolding.sh <version> <systematic>"
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
#method=('svd')
#method=('bayes')
method=('bayes' 'svd')
#method=('bayes' 'svd' 'inv' 'bin')
#method="bayes"
END=20

echo "2DUnfolding for" $VAR

for ((i=0; i< ${#method[@]};i++))
do
    if [ ${method[$i]} == 'bayes' ] || [ ${method[$i]} == 'svd' ]; then
        for j in $(seq 1 $END)
        do
            python /usatlas/u/goyeonju/phoTaggedJetRaa/unfolding/unfold_jetPt_for2DUnfolding_v2.py PbPb ${method[$i]} $j $VER $SYS &
            python /usatlas/u/goyeonju/phoTaggedJetRaa/unfolding/unfold_jetPt_for2DUnfolding_v2.py PP ${method[$i]} $j $VER $SYS & 
        done
    else
        python /usatlas/u/goyeonju/phoTaggedJetRaa/unfolding/unfold_jetPt_for2DUnfolding_v2.py PbPb ${method[$i]} 0 $VER $SYS &
        python /usatlas/u/goyeonju/phoTaggedJetRaa/unfolding/unfold_jetPt_for2DUnfolding_v2.py PP ${method[$i]} 0 $VER $SYS &
    fi
done
wait $(jobs -p)

echo 'DONE 2DUnfolding'
mv /direct/usatlas+u/goyeonju/phoTaggedJetRaa/unfolding/figures/*.pdf /direct/usatlas+u/goyeonju/phoTaggedJetRaa/unfolding/figures/backup

for ((i=0; i< ${#method[@]};i++))
do
    if [ ${method[$i]} == 'bayes' ] || [ ${method[$i]} == 'svd' ]; then
        for j in $(seq 1 $END)
        do
            root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/unfolding/drawUnfoldingFigures_2DUnfolding_v2.C("PP", "'${VER}'", "'${SYS}'", "'${method[$i]}'", '$j')' &
            root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/unfolding/drawUnfoldingFigures_2DUnfolding_v2.C("PbPb", "'${VER}'", "'${SYS}'", "'${method[$i]}'", '$j')' &
        done
    else
        root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/unfolding/drawUnfoldingFigures_2DUnfolding_v2.C("PP", "'${VER}'", "'${SYS}'", "'${method[$i]}'", 0)' &
        root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/unfolding/drawUnfoldingFigures_2DUnfolding_v2.C("PbPb", "'${VER}'", "'${SYS}'", "'${method[$i]}'", 0)' &
    fi
done

wait $(jobs -p)

for ((i=0; i< ${#method[@]};i++))
do
    if [ ${method[$i]} == 'bayes' ]; then
        for j in $(seq 1 $END)
        do
          root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/unfolding/unfoldingDiagnostics_2DUnfolding_v2.C("PP", "'${VER}'", "'${SYS}'", "'${method[$i]}'",'$END')' &
          root -l -b -q '/direct/usatlas+u/goyeonju/phoTaggedJetRaa/unfolding/unfoldingDiagnostics_2DUnfolding_v2.C("PbPb", "'${VER}'", "'${SYS}'", "'${method[$i]}'",'$END')' &
        done
    fi
done
