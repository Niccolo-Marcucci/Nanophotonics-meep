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

seff $jobid > ${log_file:0:-4}.info

mv -v *${jobid}.* data/*$hash
