ENV_name = "gwdg"
ENV_root = "parallel_gwdg"
println("Environment: GWDG")

function jobengine_head(dir::String, Ncores::Integer, gpu::Bool, number::Int64, hours::Int64=48)
  return """#!/bin/sh
  #BSUB -L /bin/bash
  #BSUB -q $(gpu ? "gpu" : "mpi")
  #BSUB -W $(hours):00
  #BSUB -e /usr/users/bardenn/reconstruction/$dir/terminal-err
  #BSUB -o /usr/users/bardenn/reconstruction/$dir/terminal-out
  #BSUB -n $Ncores
  #BSUB -R span[hosts=1]
  #BSUB -R "rusage[ngpus_shared=$Ncores]"
  #BSUB -a openmp
  #BSUB -cwd /usr/users/bardenn/reconstruction/$dir
  #BSUB -J $(jobname(dir,number))
  #$(number > 1 ? "BSUB -w \"done($(jobname(dir,number-1)))\"" : "" )
  """
end

"""Launches the script via qsub, in case dir doesnt exist, create it"""
function launch_job(dir::String, Ncores::Integer, gpu::Bool, julia_script::String, successive_jobs::Integer=1, hours::Int64=48)
  run(`mkdir -p $dir`)

  #Let's queue some jobs, hold_jid takes care of chaining
  for n = 1:successive_jobs
    jobfile = "job$n.sh"
    head = jobengine_head(dir, Ncores, gpu, n, hours)
    open("$dir/$jobfile", "w+") do file
      write(file,"""$head

      julia << EOF
        addprocs($(Ncores > 1 ? Ncores : 0)) #Let's go into parallel mode if necessary
        $julia_script
        flush(STDOUT) #make sure output is flushed to terminal files
      EOF
      """)
      flush(file)
    end
    submit_job("$dir/$jobfile")
  end
end

"""Submits a job to the queue"""
function submit_job(path)
  run(pipeline("$path", `bsub`))
end
