Pkg.init()
Pkg.add("Optim")
Pkg.add("CUDArt")
Pkg.add("CUBLAS")
Pkg.add("Distributions")

ENV["JUPYTER"] = ""
ENV["PYTHON"] = ""
Pkg.add("IJulia")
Pkg.add("PyPlot")
import Conda
Conda.add("seaborn")
Pkg.add("LaTeXStrings")

Pkg.add("HDF5")
Pkg.add("Images")
