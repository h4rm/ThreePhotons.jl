A Julia package for the structure determination from single molecule X-ray scattering experiments with down to only three photons.

Detailed information about the context of this package can be found at <link to paper>.

General usage is described in <>.

To start from a fresh Julia installation, you may run:

```julia
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
```

Package Functionality
======================

Description of 3D Structures with Spherical Harmonics Basis
-----------------------------------------------------------
This package provides functions to describe 3D structures on a cubic grid or on a spherical grid, each shell expanded in spherical harmonics. PDB structures can be loaded directly.

```julia
  LMAX = 25 #Maximum expansion order of spherical harmonics expansion
  KMAX = 30 #Maximum number of shells used in the expansion

  #Description of the Crambin electron density, Fourier density, and Fourier intensity expanded in Spherical Harmonics
  density,fourier,intensity = createSphericalHarmonicsStructure("data/structures/crambin.pdb", LMAX, KMAX, float(KMAX))
  #Same Crambin structure on a cubic grid
  densityCube,fourierCube,intensityCube = createCubicStructure("data/structures/crambin.pdb", 2*KMAX+1, float(KMAX))
```

Synthetic Scattering Images And Histogramming
----------------------------------------------
This package covers the generation of synthetic scattering images (with and without noise) and subsequent two-photon and three-photon histogramming. For computational reasons, the scattering images are not cached.

```julia
  include("jobs/runs.jl")

  generate_histograms(;
  num_images = Integer(1e6), #Number of images to generate,
  max_triplets=Integer(0), #Alternatively, the maximum number of triplets can be limited
  Ncores=8, #Number of CPU cores used for the data generation
  N=32, #Alpha/beta discretization
  photons_per_image=10, #Photons per image
  batchsize=round(Int64,1e6/8),
  successive_jobs=1,
  use_cube=false, #Use cubic or spherical harmonics description for data generation
  qcut_ratio=1.0, #Fraction of maximum wave number
  kcut=38, #total number of shells
  rmax=float(38), #maximum radius in real space
  histogram_method="histogramCorrelationsInPicture_alltoall")
```

Structure Determination
-----------------------
Given a histogrammed two- and three-photon correlation, the structure can be retrieved without additional knowledge.

```julia
  include("jobs/runs.jl")
  num_images::Int64 = Integer(1e6) #number of images
  KMAX::Int64=38 #Maximum shell number used for two-photon inversion
  N::Int64=32 #Alpha/beta discretization
  L::Int64=18 #Maximum expansion order
  K::Int64=26 #Number of shells used for structure determination
  ppi::Int64=10 #Photons per image used for the histogram
  RMAX = float(KMAX)#Maximum radius of the reference structures
  name = histogram_name("", ppi, N, KMAX, float(KMAX), img, "") for img in image_list) #histogram file name

  run_determination(
  "runname", #Name of the run
  histograms=name, #Path to the histogram file

  #Expansion parameters (see above)
  kcut=K,
  lcut=L,
  kmax=KMAX,
  rmax=RMAX,
  N=N,

  #Monte Carlo simulated annealing parameters
  initial_stepsize=pi/4.0,
  optimizer="rotate_all_at_once",
  initial_temperature_factor=0.1,
  temperature_decay=0.99998,
  stepsizefactor=1.01
  measure="Bayes",

  #Misc parameters
  range=1000:1019,
  postprocess=true,
  gpu=true,
  Ncores=20,
)
```

Dependencies
==============

Julia (5.1)

Libraries
---------

* FFTW 3.3 (for SH calculations in s2kit)
* Gnu Science Library (GSL) (for spherical harmonics basis functions (Ylm) and Wigner-3j symbols)

Julia Packages
------
* CUDArt.jl (CUDA 8.0)
* CUBLAS.jl
* Distributions.jl
* PyPlot.jl
* LaTeXStrings.jl
* ProgressMeter.jl

#### Author: Benjamin von Ardenne
#### Copyright 2017
