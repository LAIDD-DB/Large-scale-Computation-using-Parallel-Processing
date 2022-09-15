#!/usr/bin/env python
import numpy as np
from mpi4py import MPI
import sys
mpicomm = MPI.COMM_WORLD
mpirank = mpicomm.Get_rank()
mpisize = mpicomm.Get_size()

icolor = mpirank%2 + 1

mpicomm_sub = mpicomm.Split(icolor,mpirank)
mpirank_sub = mpicomm_sub.Get_rank()
mpisize_sub = mpicomm_sub.Get_size()

print( 'Hello I am', mpirank, 'in mpicomm and', mpirank_sub, 'in mpicomm_sub.' )

