#!/bin/sh
for i in `seq 0 3`;
do
    python calc_tanimoto.py tanimoto$i.csv &
done
wait
