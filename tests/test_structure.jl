#Structures

density,fourier,intensity = createSphericalHarmonicsStructure("$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb", LMAX, KMAX, float(KMAX))
densityCube,fourierCube,intensityCube = createCubicStructure("$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb", 2*KMAX+1, float(KMAX))
# @test_approx_eq_eps similarity(getCube(density), densityCube) 1.0 1e-1
# @test_approx_eq_eps similarity(getCube(intensity), intensityCube) 1.0 1e-2
