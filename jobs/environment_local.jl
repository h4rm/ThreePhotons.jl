ENV_name = "local"
ENV_root = "parallel_local"
println("Environment: Local")

"""Launches the script via qsub, in case dir doesnt exist, create it"""
function launch_job(dir::String, Ncores::Integer, gpu::Bool, julia_script::String, successive_jobs::Integer=1,  hours::Int64=48)

  cd("$(ENV["THREEPHOTONS_PATH"])")
  run(`mkdir -p $(ENV["THREEPHOTONS_PATH"])/$dir`)
  open("$(ENV["THREEPHOTONS_PATH"])/$dir/job.jl", "w+") do file
    write(file,"""
    addprocs($(Ncores > 1 ? Ncores : 0)) #Let's go into parallel mode if necessary
    $julia_script
    flush(STDOUT) #make sure output is flushed to terminal files
    """)
  end
  cd("$(ENV["THREEPHOTONS_PATH"])/$dir")
  include("$(pwd())/job.jl")
end
