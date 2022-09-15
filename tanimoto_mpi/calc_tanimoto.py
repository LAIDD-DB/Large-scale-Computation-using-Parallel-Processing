#!/usr/bin/env python
import pandas as pd
import numpy as np
from mpi4py import MPI

mpicomm = MPI.COMM_WORLD
mpirank = mpicomm.Get_rank()
mpisize = mpicomm.Get_size()


def calc_tanimoto(fp1, fp2):
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

    count = -1
    result = {}
    for i in range(len(fps)-1):
        fp1 = fps[i]
        for j in range(i+1,len(fps)):
            count += 1 
            if count % mpisize != mpirank:
                continue
            fp2 = fps[j]
            t = calc_tanimoto(fp1, fp2)
            result[count] = (mols[i], mols[j], t)

    all_result = mpicomm.gather(result, root=0)

    if mpirank == 0:
        merged_result = {}
        for x in (all_result):
            merged_result.update(x)
        result = []
        for i in sorted(merged_result):
            result.append(merged_result[i])

        df = pd.DataFrame(result)
        df.columns = ['i','j','tanimoto']
        df.to_csv('tanimoto.csv', index=False)
