#!/usr/bin/env python

import submitit

def calc(a,b):
    # do some complicated job
    return a+b


# 소량의 job 제출시

# job 생성
executor = submitit.AutoExecutor(folder='temp_dir')
executor.update_parameters(timeout_min=1)
jobs = []
for pair in ((1,2), (3,4), (5,6), (7,8)):
    i, j = pair
    job = executor.submit(calc, i, j)
    jobs.append(job)

# 결과값 수집
results = []
for job in jobs:
    results.append(job.result())


# 대량의 job 제출시 arrayjob 활용
# argument 각각이 리스트 형식으로 준비됨
a = [ 1, 3, 5, 7 ]
b = [ 2, 4, 6, 8 ]
# 동시에 2개씩 작업 (slurm_array_parallelism)
executor = submitit.AutoExecutor(folder='temp_dir')
executor.update_parameters(slurm_array_parallelism=2)
jobs = executor.map_array(calc, a, b)
for job in jobs:
    print(job.result())


# 아래의 예제는 위와 동일한 결과를 준다
# with executor.batch() 형식으로 되어 있는 것이 차이점
# 또 다른 차이점은 argument를 개별적으로 줄 수 있다는 점
executor = submitit.AutoExecutor(folder='temp_dir')
executor.update_parameters(slurm_array_parallelism=2)
jobs = []
with executor.batch():
    for pair in ((1,2), (3,4), (5,6), (7,8)):
        i, j = pair
        job = executor.submit(calc, i, j)
        jobs.append(job)
for job in jobs:
    print(job.result())
