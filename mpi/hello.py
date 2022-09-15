#!/usr/bin/env python
from mpi4py import MPI
mpicomm = MPI.COMM_WORLD
mpirank = mpicomm.Get_rank()
mpisize = mpicomm.Get_size()

print('Hello I am', mpirank)
