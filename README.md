A Julia package for the structure determination from single molecule X-ray scattering experiments at very low photon counts. Detailed information about the context of this package can be found at <link to paper>.

To start from a fresh Julia installation, you may run:

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

Package Functionality
======================



Description of 3D Structures with Spherical Harmonics Basis
-----------------------------------------------------------
This package provides functions to describe 3D structures on a cubic grid or on a spherical grid, each shell expanded in spherical harmonics.

Synthetic Scattering Images And Histogramming
----------------------------------------------
This package covers the generation of synthetic scattering images (with and without noise) and subsequent two-photon and three-photon histogramming. For computational reasons, the scattering images are not cached.

Structure Determination
-----------------------

Dependencies
==============

Libraries
---------

* FFTW 3.3 (for SH calculations)
* Gnu Science Library (GSL) (for spherical harmonics basis functions and wigner symbols)

Julia (5.1)
------
* CUDArt.jl (CUDA 8.0)
* CUBLAS.jl
* Distributions.jl
* PyPlot.jl
* LaTeXStrings.jl
* ProgressMeter.jl

#### Author: Benjamin von Ardenne
#### Copyright 2017
