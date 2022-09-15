#!/usr/bin/env python
from mpi4py import MPI
import numpy as np

mpicomm = MPI.COMM_WORLD
mpirank = mpicomm.Get_rank()
mpisize = mpicomm.Get_size()

if mpirank == 0:
    data = np.arange(10,dtype=int)
else:
    data = np.empty(10,dtype=int)
#mpicomm.Bcast([data,MPI.INT],root=0)
mpicomm.Bcast(data,root=0)

print( 'I am', mpirank, 'and I have', data )
