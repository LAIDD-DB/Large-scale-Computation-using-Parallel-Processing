#!/usr/bin/env python
import pandas as pd
import numpy as np

def calc_tanimoto(fp1, fp2):
    n1 = np.count_nonzero(fp1)
    n2 = np.count_nonzero(fp2)
    comm = np.count_nonzero(np.logical_and(fp1,fp2))

    tanimoto = comm / (n1 + n2 - comm)
    return tanimoto


df = pd.read_csv('fingerprint.csv')
fps = df['fingerprint'].to_list()     
fps = [ np.array(list(map(int,list(x)))) for x in fps ]
mols = df['mol'].tolist()

out = []
for i in range(len(fps)-1):
    fp1 = fps[i]
    for j in range(i+1,len(fps)):
        fp2 = fps[j]
        t = calc_tanimoto(fp1,fp2)
        out.append((mols[i],mols[j],t))
df = pd.DataFrame(out)
df.columns = ['i','j','tanimoto']
df.to_csv('tanimoto.csv', index=False)
