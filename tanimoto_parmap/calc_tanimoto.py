#!/usr/bin/env python
import pandas as pd
import numpy as np
import parmap


def calc_tanimoto(pair, fps):
    i, j = pair
    fp1 = fps[i]
    fp2 = fps[j]
    n1 = np.count_nonzero(fp1)
    n2 = np.count_nonzero(fp2)
    comm = np.count_nonzero(np.logical_and(fp1,fp2))
    return comm / (n1 + n2 - comm)


if __name__ == '__main__':
    NPROC = 4

    df = pd.read_csv('fingerprint.csv')
    fps = df['fingerprint'].to_list()     
    fps = [ np.array(list(map(int,list(x)))) for x in fps ]
    mols = df['mol'].tolist()

    pairs = []
    tanimoto = []
    for i in range(len(fps)-1):
        for j in range(i+1,len(fps)):
            pairs.append((i,j))
            tanimoto.append((mols[i], mols[j]))

    ret = parmap.map(calc_tanimoto, pairs, fps, pm_pbar=True, pm_processes=NPROC) 

    df = pd.DataFrame(tanimoto)
    df.columns = ['i', 'j']
    df['tanimoto'] = ret
    df.to_csv('tanimoto.csv', index=False)
