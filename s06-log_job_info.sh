#!/bin/bash

jobid=$1

log_file=$(ls *${jobid}.log)
if [ ".log" != "${log_file:$(( ${#log_fi} - 4 ))}" ]; then
    exit
fi
prefix=${log_file%_${jobid}.log}    # extract prefix
hash=${prefix:$(( ${#pref} - 10 ))} # then the hash is 10 chars

echo $hash
echo $prefix

state=$(seff 303262 | grep -i "state: ")

if [ ${val#"State: "} != "RUNNING" ]; then
    seff $jobid > ${log_file:0:-4}.info
    python mv_n_tar_after_log.py $hash
fi
# mv -v *${jobid}.* data/*$hash
