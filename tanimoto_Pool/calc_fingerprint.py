#!/usr/bin/env python

import pandas as pd
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import multiprocessing


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


if __name__ == '__main__':
    NPROC=1

    df = pd.read_csv('../tanimoto/ligands.csv')
    smiles = df['smiles'].tolist()
    split_smiles = np.array_split(smiles, NPROC)

    # 병렬 프로세스 생성하고 서브루틴 실행
    with multiprocessing.Pool(NPROC) as p:
        ret = p.map(calc_fingerprints, split_smiles)

    # 결과 값 가지고 오기
    fps = []
    for fps_split in ret:
        for fp in fps_split:
            fps.append(fp)

    df['fingerprint'] = fps
    df.to_csv('fingerprint.csv', index=False)
