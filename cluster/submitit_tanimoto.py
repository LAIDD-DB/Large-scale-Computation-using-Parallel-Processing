#!/usr/bin/env python

import submitit
import pandas as df
import pandas as pd
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np


def calc_fingerprints(smiles):
    # fingerprint score
    fps = []
    for smi in smiles:
        try:
            mol = Chem.MolFromSmiles(smi)
        except:
            return None
        # Morgan fingerprint with radius 2
        fp = AllChem.GetMorganFingerprintAsBitVect(mol,2,1024)
        arr = np.zeros((1,),dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp,arr)
        fp = ''.join(list(map(str,arr.tolist())))
        fps.append(fp)
    return fps


def calc_tanimoto(mols, fps, n_work, work_ind):
    ret = []
    count = -1
    for i in range(len(fps)-1):
        fp1 = fps[i]
        for j in range(i+1,len(fps)):
            fp2 = fps[j]
            count += 1
            if count % n_work != work_ind: 
                continue
            n1 = np.count_nonzero(fp1)
            n2 = np.count_nonzero(fp2)
            comm = np.count_nonzero(np.logical_and(fp1,fp2))

            tanimoto = comm / (n1 + n2 - comm)
            ret.append({'i': mols[i], 'j': mols[j], 'tanimoto': tanimoto})
    return ret


n_split = 4

# fingerprint calculation
df = pd.read_csv('ligands.csv')
smiles = df['smiles'].tolist()
split_smiles = np.array_split(smiles,n_split)

executor = submitit.AutoExecutor(folder='temp_dir')
executor.update_parameters(slurm_array_parallelism=2)
jobs = []
with executor.batch():
    for s_smiles in split_smiles:
        job = executor.submit(calc_fingerprints, s_smiles)
        jobs.append(job)
result = []
for job in jobs:
    result += job.result()
df['fingerprint'] = result
df.to_csv('fingerprint.csv', index=False)


# tanimoto calculation
fps = df['fingerprint'].to_list()     
fps = [ np.array(list(map(int,list(x)))) for x in fps ]
mols = df['mol'].tolist()

executor = submitit.AutoExecutor(folder='temp_dir')
executor.update_parameters(slurm_array_parallelism=2)
jobs = []
with executor.batch():
    for work_ind in range(n_split):
        job = executor.submit(calc_tanimoto, mols, fps, n_split, work_ind)
        jobs.append(job)
result = []
for job in jobs:
    result += job.result()
    print(result)

df = pd.DataFrame(result)
df.to_csv('tanimoto.csv', index=False)
