#!/bin/bash

t0=$(date +%s.%N)
python test_embarassing_parallelization.py bash 0
t1=$(date +%s.%N)
t_ex0=$(echo "$t1 - $t0" | bc)

t0=$(date +%s.%N)
python test_embarassing_parallelization.py serial
t1=$(date +%s.%N)
t_ex1=$(echo "$t1 - $t0" | bc)

python test_embarassing_parallelization.py
t2=$(date +%s.%N)
t_ex2=$(echo "$t2 - $t1" | bc)

mpirun -np 6 python test_embarassing_parallelization.py
t3=$(date +%s.%N)
t_ex3=$(echo "$t3 - $t2" | bc)

for i in 0 1 2 3 4 6
do 
python test_embarassing_parallelization.py bash $i & 
done
wait
t4=$(date +%s.%N)
t_ex4=$(echo "$t4 - $t3" | bc)

echo "done"

echo "single simulation took $t_ex0 s, serial took $t_ex1 s, Pool() took $t_ex2 s, mpirun took $t_ex3 s, bash took $t_ex4 s"
