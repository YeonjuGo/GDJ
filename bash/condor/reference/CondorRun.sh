#!/bin/bash

echo "------------------setting up environment--------------------"
export Cur_dir=$(pwd)
echo "running area:" ${Cur_dir}
echo "setup ATLAS"
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
lsetup "root 6.08.00-x86_64-slc6-gcc49-opt"
echo "-------------------------running----------------------------"
cd ${Cur_dir}
ls
root -l -b -q runMaker.C
