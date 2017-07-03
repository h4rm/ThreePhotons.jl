#llrun -f reconstruction/jobs/login_gpu_hydra.sh -w UNLIMITED
# @ shell=/bin/bash
# @ error = $(jobid).err
# @ output = $(jobid).out
# @ node = 1
# @ tasks_per_node = 1
# @ node_resources = ConsumableMemory(1gb) ConsumableCpus(1)
# @ requirements = (Feature=="gpu")
# @ task_affinity = core(1)
# @ node_usage = shared
# @ wall_clock_limit = 24:00:00
# @ initialdir = /u/bardenn
# @ queue

julia -e "Pkg.test(\"CUDArt\")"
