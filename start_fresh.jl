Pkg.init()
ENV["JUPYTER"] = ""
ENV["PYTHON"] = ""
Pkg.add("IJulia")
Pkg.add("PyPlot")
Pkg.add("ProgressMeter")
Pkg.add("Optim")
Pkg.add("CUDArt")
Pkg.add("CUBLAS")
Pkg.add("Distributions")
Pkg.add("HDF5")
Pkg.add("Images")
