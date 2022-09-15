#!/bin/sh
./split_data.py
./calc_fingerprint_parallel.sh
./merge_csv.py merged_fingerprint.csv fingerprint*.csv
./split_tanimoto.py
./calc_tanimoto_parallel.sh
./merge_csv.py merged_tanimoto.csv tanimoto*.csv
