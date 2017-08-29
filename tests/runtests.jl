using Base.Test
addprocs() #Can be turned on for faster calculations
using ThreePhotons

function test_file(filename)
  print_with_color(:green, "Testing: $filename\n")
  include(filename)
end

#General test parameters
K = 20
L = 16
N = 16
KMAX = 38
LMAX = 25

test_file("test_structure.jl")
# test_file("test_cubic.jl")
# test_file("test_spherical_harmonics.jl")
# test_file("test_datagen.jl")
# test_file("test_correlations.jl")
test_file("test_cuda.jl")
# test_file("test_determination.jl")
# test_file("test_phases.jl")
# test_file("test_data_processing.jl")
# test_file("test_utilities.jl")

print_with_color(:green, "Testing complete!\n")
