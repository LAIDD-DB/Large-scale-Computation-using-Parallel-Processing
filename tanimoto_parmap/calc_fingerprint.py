#!/usr/bin/env python

import pandas as pd
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import parmap


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
    NPROC=4

    df = pd.read_csv('../tanimoto/ligands.csv')
    smiles = df['smiles'].tolist()

    fps = parmap.map(calc_fingerprint, smiles, pm_processes=NPROC)

    df['fingerprint'] = fps
    df.to_csv('fingerprint.csv', index=False)
