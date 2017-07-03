#General test parameters
K = 10
L = 10
N = 10
KMAX = 12
LMAX = 11

#Structures

density,fourier,intensity = createSphericalHarmonicsStructure("../data/structures/crambin.pdb", LMAX, KMAX, float(KMAX))
densityCube,fourierCube,intensityCube = createCubicStructure("../data/structures/crambin.pdb", 2*KMAX+1, float(KMAX))
# @test_approx_eq_eps similarity(getCube(density), densityCube) 1.0 1e-1
# @test_approx_eq_eps similarity(getCube(intensity), intensityCube) 1.0 1e-2
