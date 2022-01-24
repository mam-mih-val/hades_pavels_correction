#!/bin/bash

format='+%Y/%m/%d-%H:%M:%S'

date $format

job_num=$(($SLURM_ARRAY_TASK_ID))

filelist=$lists_dir/$(ls $lists_dir | sed "${job_num}q;d")

cd $output_dir
mkdir -p $job_num
cd $job_num

qn_tools=/lustre/nyx/hades/user/mmamaev/QnAnalysis/build-centos7/src
root_620=/lustre/nyx/hades/user/mmamaev/install/root-6.20.04-centos7-cxx17/bin/thisroot.sh

module load /cvmfs/vae.gsi.de/centos7/modules/linux-centos7-x86_64/gcc-8.1.0-gcc-4.8.5-oyp4lmr

echo "loading " $root_620
source $root_620

$qn_tools/QnAnalysisCorrect/QnAnalysisCorrect -i $filelist \
                                              -t hades_analysis_tree \
                                              --yaml-config-file=/lustre/nyx/hades/user/mmamaev/QnAnalysis/setups/hades/correction-auau-123-plains.yml \
                                              --yaml-config-name=hades_analysis \
                                              -n 1000 \
                                              --cuts-macro Hades/AuAu1.23.C \
                                              --event-cuts hades/auau/1.23/event_cuts/standard/pt3

echo "executing $build_dir/yield -i $filelist \
                                -t hades_analysis_tree \
                                -n 1000 -o yield.root \
                                --cuts-macro Hades/AuAu1.23.C \
                                --pdg-code $pdg_code \
                                --q-vector-file=correction_out.root \
                                --q-vector-name=W1_PLAIN \
                                --cuts-macro Hades/AuAu1.23.C \
                                --event-cuts hades/auau/1.23/event_cuts/standard/pt3"

$build_dir/yield -i $filelist \
                -t hades_analysis_tree \
                -n 1000 -o yield.root \
                --cuts-macro Hades/AuAu1.23.C \
                --pdg-code $pdg_code \
                --q-vector-file=correction_out.root \
                --q-vector-name=W1_PLAIN \
                --cuts-macro Hades/AuAu1.23.C \
                --event-cuts hades/auau/1.23/event_cuts/standard/pt3

date $format
echo JOB FINISHED!