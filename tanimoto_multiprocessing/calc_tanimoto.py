#!/usr/bin/env python
import pandas as pd
import numpy as np
import multiprocessing


def calc_tanimoto(pairs, mols, fps, proc_id, ret):
    tanimoto = []
    for i, j in pairs:
        fp1 = fps[i]
        fp2 = fps[j]

        i_mol = mols[i]
        j_mol = mols[j]

        n1 = np.count_nonzero(fp1)
        n2 = np.count_nonzero(fp2)
        comm = np.count_nonzero(np.logical_and(fp1,fp2))
        t = comm / (n1 + n2 - comm)
        tanimoto.append((i, j, i_mol, j_mol, t))
    ret[proc_id] = tanimoto


if __name__ == '__main__':
    NPROC = 4

    df = pd.read_csv('fingerprint.csv')
    fps = df['fingerprint'].to_list()     
    fps = [ np.array(list(map(int,list(x)))) for x in fps ]
    mols = df['mol'].tolist()

    manager = multiprocessing.Manager()
    ret = manager.dict()

    # pair 나누기
    split_pairs = [ [] for x in range(NPROC) ]
    count = 0
    for i in range(len(fps)-1):
        for j in range(i+1,len(fps)):
            ind = count % NPROC
            split_pairs[ind].append((i,j))
            count += 1

    procs = []
    for proc_id in range(NPROC):
        pairs = split_pairs[proc_id]
        proc = multiprocessing.Process(target=calc_tanimoto,
                args=(pairs, mols, fps, proc_id, ret))
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()

    tanimoto = []
    for key in sorted(ret.keys()):
        tanimoto_split = ret[key]
        for t in tanimoto_split:
            tanimoto.append(t)
             
    df = pd.DataFrame(tanimoto)
    df.columns = ['ind1', 'ind2', 'i','j','tanimoto']
    df = df.sort_values(by=['ind1','ind2'])[['i','j','tanimoto']]
    df.to_csv('tanimoto.csv', index=False)

