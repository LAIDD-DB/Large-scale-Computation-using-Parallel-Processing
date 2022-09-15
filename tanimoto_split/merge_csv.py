#!/usr/bin/env python
import sys
import pandas as pd

outcsv = sys.argv[1]

dfs = []
for csv in sorted(sys.argv[2:]):
    df = pd.read_csv(csv)
    dfs.append(df)

df = pd.concat(dfs)
df.to_csv(outcsv, index=False)

