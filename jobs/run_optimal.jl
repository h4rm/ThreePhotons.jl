"""Calculates the resolution of corresponding (dq,K,L) combination"""
function run_optimal(K2::Int64, K3::Int64, L3::Int64)
    name = "optimal/optimal_K2_$(K2)_K3_$(K3)_L3_$(L3)"

    julia_script = """
    using ThreePhotons
    deltar = 1.0 #2 Angstrom resolution maximum
    qmax2 = pi/deltar
    K2 = $K2
    rmax = K2*deltar
    deltaq = qmax2 / K2

    K3 = $K3
    L3 = $L3
    _,_,intensity = createSphericalHarmonicsStructure("$(ENV["THREEPHOTONS_PATH"])/structures/crambin.pdb", 25, K2, rmax)
    densityCube,fourierCube,intensityCube = createCubicStructure("$(ENV["THREEPHOTONS_PATH"])/structures/crambin.pdb", 2*K2+1, rmax)
    calculate_optimal(intensity,densityCube,fourierCube,intensityCube, K3, L3, ".", 16)

    #Now check the same with retrieving the structure
    checkRotation_basis = complexBasis(L3,32,25)
    c2hist_full = twoPhotons(intensity, checkRotation_basis, intensity.KMAX, true, false)
    start = retrieveSolution(c2hist_full, checkRotation_basis.L, intensity.LMAX, intensity.KMAX, intensity.rmax)
    energyfunc = (volume) -> 0.0
    reference_intensity = deleteTerms(intensity, K3, checkRotation_basis.L)
    rightCoeff,finalFSC = checkRotationSearch(start, reference_intensity, K3, checkRotation_basis.L, checkRotation_basis, save_structures=false, energy = energyfunc, iterations=4.0e4, reduce_stepsize=3000)
    saveCube(rightCoeff, "./retrievedIntensity.mrc")

    #Phase retrieved structure
    res = calculateSC(deleteTerms(rightCoeff, K3,checkRotation_basis.L), densityCube, fourierCube, intensityCube, 16)
    saveCube(res[1], "./retrievedDensity.mrc")
    serializeToFile("./retrieved.dat", res)
    """
    launch_job(name, 8, false, julia_script, 1)
end

function run_all_optimal()
    #K2 must be at least such that rmax > 35.0 Angstrom (Phasing)
    deltar = 1.0 #2 Angstrom resolution maximum
    qmax2 = pi/deltar
    K2_min = ceil(Int64, 35.0 / deltar)
    for K2 = K2_min:3:K2_min+12
      deltaq = qmax2 / K2
      for K3 = 20:2:min(K2,40)
        for L3 = 12:2:floor(Int64,(K2-1)/2)
          run_optimal(K2, K3, L3)
        end
      end
    end
end
