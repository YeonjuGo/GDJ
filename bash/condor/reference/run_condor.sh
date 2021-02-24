#!/bin/bash

export TargetDir="$PWD"/condorRun

rm inputdata.txt

if [ -d ${TargetDir} ]; then
  rm -rf ${TargetDir}/OutDir*
else
  mkdir ${TargetDir}
fi


#for file in /usatlas/scratch/bseidlit/pp_decor_run0973/user.bseidlit.ppDecorrelation_07.00340973.r11215_p3764_ANALYSIS.root/*
for file in /usatlas/scratch/bseidlit/pp_5TeV_w_tracksTowers/user.bseidlit.ppDecorrelation_13.00340973.r11215_p3764_ANALYSIS.root/*.root
do
cat >>inputdata.txt<< EOF
$file
EOF
done


j=252
tot_files=$( cat inputdata.txt | wc -l )
echo "total files: $tot_files"
rem=$(( $tot_files%$j ))
files_per_job=$(( $tot_files/$j ))
njob=$j
if [ $rem -ne 0 ]; then
  files_per_job=$(( $files_per_job+1 ))
fi
rem2=$(( $tot_files%$files_per_job ))
njob=$(( $tot_files/$files_per_job ))
if [ $rem2 -ne 0 ]; then
  njob=$(( ($tot_files/$files_per_job)+1 ))
fi
echo "files per job: $files_per_job"
echo "njob: $njob"

for((i=0;i<$njob;i++));
do

  mkdir ${TargetDir}/OutDir$i
  export WorkDir="${TargetDir}/OutDir$i"
  echo "WorkDir:" ${WorkDir}
  start_file=$(( $i*$files_per_job+1 ))
  end_file=$(( $start_file+$files_per_job-1 ))
  echo "start file: $start_file   end file: $end_file"

  sed -n $start_file\,${end_file}p inputdata.txt > tmp.txt
  mv tmp.txt ${WorkDir}/inputdata.txt

  pushd ${WorkDir}

  cp "$PWD"/../../CondorRun.sh CondorRunTC$i.sh
  cp "$PWD"/../../runMaker.C .
  cp "$PWD"/../../skim_pp.* .
  cp "$PWD"/../../trkEff_nominal.root .
  cp "$PWD"/../../pp_towers_gt4.0_array.txt .
  cp "$PWD"/../../pp_towers_gt4.0_array_phi.txt .
  cp "$PWD"/../../pp_towers_calib_cuts.txt .
  cp "$PWD"/../../rootFiles/pp_tower_calib_cuts.root .

  cat >>ff.sub<< EOF
+JobFlavour                   = "workday"
transfer_input_files          = ${WorkDir}/runMaker.C,${WorkDir}/inputdata.txt,${WorkDir}/CondorRunTC$i.sh,${WorkDir}/skim_pp.C,${WorkDir}/trkEff_nominal.root,${WorkDir}/pp_towers_gt4.0_array.txt,${WorkDir}/pp_towers_gt4.0_array_phi.txt,${WorkDir}/pp_towers_calib_cuts.txt,${WorkDir}/pp_tower_calib_cuts.root
Executable                    = CondorRunTC$i.sh
Universe                      = vanilla
Notification                  = Never
GetEnv                        = True
Priority                      = +20
Output                        = test.out
Error                         = test.err
Log                           = test.log
Notify_user                   = blair.daniel.seidlitz@cern.ch
accounting_group              = group_atlas.boulder

Queue
EOF

  condor_submit ff.sub
  popd
done
