#!/usr/bin/env python 
from mpi4py import MPI
import numpy as np

mpicomm = MPI.COMM_WORLD
mpisize = mpicomm.Get_size()
mpirank = mpicomm.Get_rank()

data = mpirank

data = mpicomm.gather(data,root=0)
#data = mpicomm.allgather(data)

print('rank', mpirank, data)
