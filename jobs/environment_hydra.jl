# ENV_name = "hydra"
# ENV_root = "output_hydra"
# println("Environment: Hydra")
#
# function jobengine_head(dir::String, Ncores::Integer, gpu::Bool, num_steps::Integer)
#
#   function single_step(n::Integer)
#     """
#     # @ step_name = iteration$n
#     $(n > 1 ? "# @ dependency = (iteration$(n-1) == 0)" : "" )
#     # @ node = 1
#     # @ tasks_per_node = 1
#     # @ resources = ConsumableCpus($Ncores)
#     $(gpu ? "# @ requirements = (Feature==\"gpu\")" : "")
#     # @ wall_clock_limit = 24:00:00
#     # @ initialdir = /u/bardenn/reconstruction/$dir
#     # @ notification = never
#     # @ queue
#
#     """
#   end
#
#   steps = [single_step(n) for n = 1:num_steps]
#
#   head = """
#   # @ shell=/bin/bash
#   # @ job_name = $(jobname(dir))
#   # @ error = /u/bardenn/reconstruction/$dir/terminal-err
#   # @ output = /u/bardenn/reconstruction/$dir/terminal-out
#   # @ job_type = parallel
#
#   $(join(steps))
#   """
# end
#
# """Launches the script via qsub, in case dir doesnt exist, create it"""
# function launch_job(dir::String, Ncores::Integer, gpu::Bool, julia_script::String, successive_jobs::Integer=1)
#   githead = check_git_status()
#   head = jobengine_head(dir, Ncores, gpu, successive_jobs)
#
#   run(`mkdir -p $dir`)
#   open("$dir/job.sh", "w+") do file
#     write(file,"""#git-SHA: $githead
#     $head
#
#     poe julia << EOF
#       addprocs($(Ncores > 1 ? Ncores : 0)) #Let's go into parallel mode if necessary
#       $julia_script
#       flush(STDOUT) #make sure output is flushed to terminal files
#     EOF
#     """)
#   end
#
#   run(`llsubmit $dir/job.sh`)
# end
