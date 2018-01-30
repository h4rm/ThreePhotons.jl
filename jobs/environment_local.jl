ENV_name = "local"
ENV_root = "output"
println("Environment: Local")

"""Launches the script via qsub, in case dir doesnt exist, create it"""
function launch_job(dir::String, Ncores::Integer, gpu::Bool, julia_script::String, successive_jobs::Integer=1; architecture::String="", hours::Int64=48, fresh::Bool=false, memory::String="", gpu_mem::String="")
  full_dir = "$(ENV["DETERMINATION_DATA"])/$(ENV_root)/$dir"
  run(`mkdir -p $(full_dir)`)
  open("$(full_dir)/job.jl", "w+") do file
    write(file,"""
    $julia_script
    flush(STDOUT) #make sure output is flushed to terminal files
    """)
  end
  cd("$(full_dir)")
  include("$(pwd())/job.jl")
end
