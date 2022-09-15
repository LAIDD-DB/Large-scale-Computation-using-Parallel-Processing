#!/usr/bin/env python
import numpy as np
from mpi4py import MPI
import sys

mpicomm = MPI.COMM_WORLD
mpirank = mpicomm.Get_rank()
mpisize = mpicomm.Get_size()

def main():
    if mpirank == 0: 
        # master prepare the jobs
        jobs = np.array( prepare_jobs(), dtype=object )
        print( 'Master prepares jobs' )
        print(jobs)
        split_job = np.array_split(jobs, mpisize)
        my_jobs = mpicomm.scatter(split_job,root=0)
    else:
        my_jobs = mpicomm.scatter(None,root=0)

    # perform jobs
    my_results = []
    for job in my_jobs:
        my_results.append( perform_job(job) )

    # gather results
    results = mpicomm.gather(my_results,root=0)
    if mpirank == 0:
        output = []
        for result in results:
            output += result
        for i in range(len(jobs)):
            print('job', jobs[i], 'result', output[i])

def perform_job(job):
    op, n1, n2 = job
    if op == '+':
        out = n1+n2
    elif op == '-':
        out = n1-n2
    elif op == '*':
        out = n1*n2
    elif op == '/':
        out = n1/n2
    print('Rank', mpirank, 'calculates', n1, op, n2, 'and return', out)
    return out

def prepare_jobs():
    operations = '+-*/'
    jobs = []
    for i in range(10):
        iop = np.random.randint(4)
        op = operations[iop]
        n1 = np.random.randint(1,11)
        n2 = np.random.randint(1,11)
        jobs.append([op,n1,n2])
    return jobs

if __name__ == '__main__':
    main()
