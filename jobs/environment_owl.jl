ENV_name = "owl"
ENV_root = "output_owl"
println("Environment: OWL")

function jobengine_head(dir::String, Ncores::Integer, gpu::Bool; hours::Int64=48, architecture::String="ivy-bridge|sandy-bridge|haswell|broadwell|skylake")
  return """
  #!/bin/bash
  #\$ -S /bin/bash
  $(Ncores > 1 ? "#\$ -pe openmp_fast $Ncores" : "")
  #\$ -q *
  #\$ -e $(DETERMINATION_DATA)/$(ENV_root)/$dir/terminal-err
  #\$ -o $(DETERMINATION_DATA)/$(ENV_root)/$dir/terminal-out
  #\$ -l cm=$(architecture)
  #\$ -N $(jobname(dir))
  #\$ -M bardenn@gwdg.de
  #\$ -m n
  #\$ -l h_rt=$(hours):00:00
  #\$ -wd $(DETERMINATION_DATA)/$(ENV_root)/$dir
  #\$ -V
  #\$ -hold_jid $(jobname(dir))
  $(gpu ? "#\$ -l gpu=1" : "")
  """
end

"""Launches the script via qsub, in case dir doesnt exist, create it"""
function launch_job(dir::String, Ncores::Integer, gpu::Bool, julia_script::String, successive_jobs::Integer=1; hours::Int64=48, architecture::String="ivy-bridge|sandy-bridge|haswell|broadwell|skylake")
  githead = check_git_status()
  head = jobengine_head(dir, Ncores, gpu; hours=hours, architecture=architecture)
  full_dir = "$(DETERMINATION_DATA)/$(ENV_root)/$dir"

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
