#!/usr/bin/env python
import numpy as np
from mpi4py import MPI
import sys
mpicomm = MPI.COMM_WORLD
mpirank = mpicomm.Get_rank()
mpisize = mpicomm.Get_size()

def main():
    if mpisize < 2:
        sys.exit()

    if mpirank == 0:
        master_work()
    else:
        slave_work()

def master_work():
    # master prepare the jobs
    jobs = prepare_jobs()
    print('Master prepares jobs')
    print(jobs)

    work_assigned = np.empty(mpisize,dtype=int)
    work_assigned[:] = -1

    output = [None]*len(jobs)
    ijob = 0
    # send jobs and recv results
    while ijob < len(jobs):
        for irank in range(1,mpisize):
            if work_assigned[irank] < 0:
                mpicomm.send(jobs[ijob],dest=irank)
                work_assigned[irank] = ijob
                ijob += 1
                if ijob >= len(jobs):
                    break
        status = MPI.Status()
        result = mpicomm.recv(source=MPI.ANY_SOURCE,status=status)
        src_rank = status.Get_source()
        output[work_assigned[src_rank]] = result
        work_assigned[src_rank] = -1

    # recv results
    for irank in range(1,mpisize):
        if work_assigned[irank] < 0:
            continue
        result = mpicomm.recv(source=irank)
        output[work_assigned[irank]] = result

    # send termination signal
    for irank in range(1,mpisize):
        mpicomm.send(None,dest=irank)

    for i in range(len(jobs)):
        print('job', jobs[i], 'result', output[i])

def slave_work():
    # slave
    while True:
        job = mpicomm.recv(source=0)
        if job is None:
            break
        out = perform_job(job)
        mpicomm.send(out,dest=0)

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
        jobs.append((op,n1,n2))
    return jobs

if __name__ == '__main__':
    main()
