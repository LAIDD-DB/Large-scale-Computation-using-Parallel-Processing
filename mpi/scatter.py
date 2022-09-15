#!/usr/bin/env python 
from mpi4py import MPI
import numpy as np

mpicomm = MPI.COMM_WORLD
mpisize = mpicomm.Get_size()
mpirank = mpicomm.Get_rank()

if mpirank == 0:
    data = np.array_split(np.arange(10), mpisize)
else:
    data = None
data = mpicomm.scatter(data,root=0)

print('rank', mpirank, data)
