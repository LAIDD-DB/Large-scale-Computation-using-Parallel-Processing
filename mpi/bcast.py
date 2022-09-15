#!/usr/bin/env python
from mpi4py import MPI

mpicomm = MPI.COMM_WORLD
mpirank = mpicomm.Get_rank()
mpisize = mpicomm.Get_size()

if mpirank == 0:
    data = 'A GIFT'
    mpicomm.bcast(data,root=0)
else:
    data = mpicomm.bcast(None,root=0)

print('I am', mpirank, 'and I have', data)
