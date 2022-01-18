#!/bin/bash

jobid=$1

log_file=$(ls *${jobid}.log)

prefix=${log_file:0:-11} 	    # jobid + extension should take 11 characters
hash=${prefix:$(( ${#pref} - 10 ))} # then the hash is 10 chars

seff $jobid > ${log_file:0:-4}.info

mv -v *${jobid}.* data/*$hash
