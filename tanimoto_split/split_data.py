#!/usr/bin/env python
import numpy as np
import pandas as pd

df = pd.read_csv('../tanimoto/ligands.csv')
inds = np.split(np.arange(len(df)),4)

for i, ind in enumerate(inds):
    df.iloc[ind].to_csv(f'ligands{i}.csv', index=False)
