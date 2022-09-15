#!/usr/bin/env python
import numpy as np
import pandas as pd

df = pd.read_csv('./merged_fingerprint.csv')
mols = df['mol'].tolist()

out = []
for i in range(len(df)-1):
    mol_i = mols[i]
    for j in range(i+1,len(df)):
        mol_j = mols[j]
        out.append((mol_i,mol_j))

df = pd.DataFrame(out)
df.columns = ['i','j']
inds = np.split(np.arange(len(df)),4)

for i, ind in enumerate(inds):
    df.iloc[ind].to_csv(f'tanimoto{i}.csv', index=False)
