ENV_name = "owl"
ENV_root = "output_owl"
println("Environment: OWL")

function jobengine_head(name::String, dir::String, Ncores::Integer, gpu::Bool; hours::Int64=48, architecture::String="ivy-bridge|sandy-bridge|haswell|broadwell|skylake", memory::String="")
  return """
  #!/bin/bash
  #\$ -S /bin/bash
  $(Ncores > 1 ? "#\$ -pe openmp_fast $Ncores" : "")
  #\$ -q *
  #\$ -e $dir/terminal-err
  #\$ -o $dir/terminal-out
  #\$ -l cm=$(architecture)
  #\$ -N $(jobname(name))
  #\$ -M bardenn@gwdg.de
  #\$ -m n
  #\$ -l h_rt=$(hours):00:00
  #\$ -wd $dir
  #\$ -V
  #\$ -hold_jid $(jobname(name))
  $(memory != "" ? "#\$ -l h_data=$memory" : "")
  $(gpu ? "#\$ -l gpu=1" : "")
  """
end

"""Launches the script via qsub, in case dir doesnt exist, create it"""
function launch_job(dir::String, Ncores::Integer, gpu::Bool, julia_script::String, successive_jobs::Integer=1; hours::Int64=48, architecture::String="ivy-bridge|sandy-bridge|haswell|broadwell|skylake", memory::String="")
  githead = check_git_status()
  full_dir = "$(ENV["DETERMINATION_DATA"])/$(ENV_root)/$dir"
  head = jobengine_head(dir, full_dir, Ncores, gpu; hours=hours, architecture=architecture, memory=memory)

  run(`mkdir -p $(full_dir)`)
  open("$(full_dir)/job.sh", "w+") do file
    write(file,"""#git-SHA: $githead
    $head

    julia << EOF
      addprocs($(Ncores > 1 ? 8 : 0)) #Let's go into parallel mode if necessary
      $julia_script
      flush(STDOUT) #make sure output is flushed to terminal files
    EOF
    """)
  end

  #Let's queue some jobs, hold_jid takes care of chaining
  for n = 1:successive_jobs
    submit_job("$(full_dir)/job.sh")
  end
end

"""Submits a job to the queue"""
function submit_job(path)
  run(`qsub $path`)
end
