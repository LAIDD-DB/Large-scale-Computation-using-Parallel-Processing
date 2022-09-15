#!/bin/sh
for i in `seq 0 3`;
do
    python calc_fingerprint.py ligands$i.csv fingerprint$i.csv &
done
wait
