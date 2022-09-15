#!/bin/sh

# job 이름
#PBS -N name_of_the_job

# stdout 파일
#PBS -o stdout.log

# stdout stderr 출력을 하나의 파일로 만드는 설정
#PBS -j oe

#PBS -l walltime=24:00:00

# cpu 할당 갯수 설정 방법
# 클러스터에 따라서 변수 이름이 다를 수 있으므로 클러스터에서
# 제공하는 지침 확인 필요
#PBS -l ncpus=12

cat $PBS_NODEFILE

# core 갯수를 카운트 하여 mpi 를 실행하는 방법
# 멀티 노드 (노드 2개 이상 할당) 시에 적용 가능한지 확인 필요
NCPU=`wc -l < $PBS_NODEFILE`
mpirun -np $NCPU your_mpi_job

# PBS_NODEFILE 변수를 직접 활용하는 방법 mpi 프로그램에 따라서 옵션
# 설정 방법이 다를 수 있다
mpirun -f $PBS_NODEFILE your_mpi_job
