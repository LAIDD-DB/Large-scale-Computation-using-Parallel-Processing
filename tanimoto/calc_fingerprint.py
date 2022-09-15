#!/usr/bin/env python

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


df = pd.read_csv('ligands.csv')
smiles = df['smiles'].tolist()

fps = calc_fingerprints(smiles)

df['fingerprint'] = fps
df.to_csv('fingerprint.csv', index=False)
