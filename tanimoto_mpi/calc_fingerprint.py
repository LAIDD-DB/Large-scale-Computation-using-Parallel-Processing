#!/usr/bin/env python

import pandas as pd
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from mpi4py import MPI

mpicomm = MPI.COMM_WORLD
mpirank = mpicomm.Get_rank()
mpisize = mpicomm.Get_size()


def calc_fingerprint(smi):
    # fingerprint score
    try:
        mol = Chem.MolFromSmiles(smi)
    except:
        return None
    # Morgan fingerprint with radius 2
    fp = AllChem.GetMorganFingerprintAsBitVect(mol,2,1024)
    arr = np.zeros((1,),dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp,arr)
    fp = ''.join(list(map(str,arr.tolist())))

    return fp


if __name__ == '__main__':
    df = pd.read_csv('../tanimoto/ligands.csv')
    smiles = df['smiles'].tolist()

    fps_dict = {}
    for i, smi in enumerate(smiles):
        if i % mpisize != mpirank:
            continue
        fp = calc_fingerprint(smi)
        fps_dict[i] = fp

    all_fps = mpicomm.gather(fps_dict, root=0)

    if mpirank == 0:
        merged_fps = {}
        for x in (all_fps):
            merged_fps.update(x)
        fps = []
        for i in sorted(merged_fps):
            fps.append(merged_fps[i])

        df['fingerprint'] = fps
        df.to_csv('fingerprint.csv', index=False)
