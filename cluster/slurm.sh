#!/bin/sh
#SBATCH -J name_of_the_job
# 노드 갯수
#SBATCH -N 1
# cpu 갯수
#SBATCH -n 12
# stdout 출력 및 stderr 출력
#SBATCH -o output.log
#SBATCH -e output.log

# mpirun 대신 srun 적용 가능
srun your_mpi_job
