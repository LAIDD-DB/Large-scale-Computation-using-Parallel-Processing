#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np

def calc_tanimoto(fp1, fp2):
    n1 = np.count_nonzero(fp1)
    n2 = np.count_nonzero(fp2)
    comm = np.count_nonzero(np.logical_and(fp1,fp2))

    tanimoto = comm / (n1 + n2 - comm)
    return tanimoto


df = pd.read_csv('./merged_fingerprint.csv')
fps = df['fingerprint'].tolist()     
fps = [ np.array(list(map(int,list(x)))) for x in fps ]
mols = df['mol'].tolist()

fps_dict = { mols[x]: fps[x] for x in range(len(mols)) }

csv = sys.argv[1]
df = pd.read_csv(csv)

out = []
for i, row in df.iterrows():
    i_mol = row['i']
    j_mol = row['j']
    t = calc_tanimoto(fps_dict[i_mol], fps_dict[j_mol])
    out.append((i_mol, j_mol, t))

df = pd.DataFrame(out)
df.columns = ['i','j','tanimoto']
df.to_csv(csv, index=False)
