#!/usr/bin/env python

import pandas as pd
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import multiprocessing


def calc_fingerprints(smiles, proc_id, return_dict):
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

    return_dict[proc_id] = fps


# 아래의 __name__ == '__main__' 확인하는 부분이 필요하다.
# 멀티 프로세스 모듈이 스크립트를 중복 실행하기 때문에 처음 실행되는
# 경우에만 아래의 코드가 작동되도록 해야 하기 때문이다.
if __name__ == '__main__':
    NPROC=4

    df = pd.read_csv('../tanimoto/ligands.csv')
    smiles = df['smiles'].tolist()
    split_smiles = np.array_split(smiles, NPROC)

    # 결과 값을 저장하기 위한 manager 생성
    manager = multiprocessing.Manager()
    return_dict = manager.dict()

    # 병렬 프로세스 생성하고 서브루틴 실행
    procs = []
    for proc_id, smiles in enumerate(split_smiles):
        proc = multiprocessing.Process(target=calc_fingerprints,
                args=(smiles, proc_id, return_dict))
        procs.append(proc)
        proc.start()

    # 생성된 프로세스 종료하기
    for proc in procs:
        proc.join()

    # 결과 값 가지고 오기
    fps = []
    for key in sorted(return_dict.keys()):
        fps_split = return_dict[key]
        for fp in fps_split:
            fps.append(fp)

    df['fingerprint'] = fps
    df.to_csv('fingerprint.csv', index=False)
