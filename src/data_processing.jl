#data processing
export fitStructures_full, fitStructures_random, getAllEvenRotations,
    sRAAR_bestfit,
    calculateSC,
    postprocess_run,
    get_qrange,
    calculate_cutoff,
    calculate_maximum_resolution,
    load_runs,
    analyse_ensemble,
    load_fitted_intensity,
    load_density,
    average_intensities,
    calculate_optimal,
    load_parameter_list,
    get_parameter_entry,
    get_optimal,
    load_sc_vs_triplets,
    substract_gamma,
    radial_difference_with_intensity,
    denoise_structure

"Fits structure `volume` to structure `reference` by exploring all possible SO(3) rotations"
function fitStructures_full(volume::SphericalHarmonicsVolume, reference::SphericalHarmonicsVolume, num::Int64, K_range::UnitRange{Int64}, L::Int64, ta::Float64 = 0.0, tb::Float64 = float(pi)/2, pa::Float64 =0.0, pb::Float64 = float(pi), ga::Float64 = 0.0, gb::Float64 = float(pi))
    println("Start with fitting")
    #check a single angle combination
    check_angles = function(theta::Float64,phi::Float64, gamma::Float64)
        rvolume = rotateStructure(volume, theta, phi, gamma, maximum(K_range), 2:2:L)
        sc = similarity(reference,rvolume, K_range)
        return (sc, theta, phi, gamma)
    end

    #generate all angle combinations
    angles = [ (theta, phi, gamma) for theta = ta:(tb-ta)/num:tb for phi = pa:(pb-pa)/num:pb for gamma = ga:(gb-ga)/num:gb]

    #check shell correlation for each angle combination
    results = pmap((angles)->check_angles(angles[1], angles[2], angles[3]), angles)

    #Retrieve the best result
    sort!(results, lt=(a,b)->a[1]>b[1])
    bestsc, bestt, bestp, bestg = results[1]
    println("Best intensity correlation: ", bestsc)
    return (rotateStructure(volume, bestt, bestp, bestg, maximum(K_range), 2:2:L), bestsc, bestt, bestp, bestg)
end

function fitStructures_full(volume::SphericalHarmonicsVolume, reference::SphericalHarmonicsVolume, num::Int64, K::Int64, L::Int64, ta::Float64 = 0.0, tb::Float64 = float(pi)/2, pa::Float64 =0.0, pb::Float64 = float(pi), ga::Float64 = 0.0, gb::Float64 = float(pi))
    fitStructures_full(volume, reference, num, 1:K, L, ta, tb, pa, pb, ga, gb)
end

"Fits structure `volume` to structure `reference` by exploring all possible SO(3) rotations"
function fitStructures_random(volume::Volume, reference::Volume, K_range::UnitRange{Int64}, L::Int64, stepsizefactor::Float64, repititions::Int64=4)
    stepsize = pi
    # results = []
    results = @parallel vcat for i = 1:repititions
        phi,theta,gamma = 2*pi*rand(), acos(rand()), 2*pi*rand()
        bestsc = 0.0

        while stepsize > pi/180.0
            new_theta = theta + randn()*stepsize/2
            new_phi = phi + randn()*stepsize
            new_gamma = gamma + randn()*stepsize
            rvolume = rotateStructure(volume, new_theta, new_phi, new_gamma, maximum(K_range), 2:2:L)
            sc = similarity(reference,rvolume, K_range)
            if sc > bestsc
                bestsc,theta,phi,gamma = sc, new_theta, new_phi, new_gamma
            end
            stepsize *= stepsizefactor
        end
        bestsc, theta, phi, gamma
    end
    sort!(results, lt=(a,b)->a[1]>b[1])
    bestsc,theta,phi,gamma = results[1]
    println("Best intensity correlation: ", bestsc)
    return (rotateStructure(volume, theta, phi, gamma, maximum(K_range), 2:2:L), bestsc, theta, phi, gamma)
end

function fitStructures_random(volume::Volume, reference::Volume, K::Int64, L::Int64, stepsizefactor::Float64, repititions::Int64=4)
    fitStructures_random(volume, reference, 1:K, L, stepsizefactor, repititions)
end

# "Calculates the power spectrum of a set of coefficients"
# function powerSpectrum(coeff)
#   sum = 0
#   for k = 1:KMAX
#     for l = 0:LMAX-1
#       # sum =  * sum
#       for m = -l:l
#         sum = sum + 1/(2*l+1)*abs(getc(coeff[k], l, m))^2
#       end
#     end
#   end
#   return sum
# end

"""Takes a structure and returns all 90 deg. rotations"""
function getAllEvenRotations(volume::CubeVolume)
    list = []
    da = pi
    for phi = 0:da:2.0*pi-da
        for theta = 0:da:2.0*pi-da
            let gamma = 0.0
                # for gamma = pi/2:pi/2:2.0*pi
                push!(list, rotateStructure(volume, phi, theta, gamma))
            end
        end
    end
    return list
end

"""Performs a spherical RAAR on the intensity and returns the best fit with the given densityCube (including mirror/rotated images)"""
function sRAAR_bestfit(referenceDensityCube::CubeVolume, intensity::SphericalHarmonicsVolume, iterations::Int64=1001, beta0::Float64 = 0.75, beta_max::Float64 = 0.90, tau::Float64 = 350.0; outputfile="", fourierCube::CubeVolume=CubeVolume(1, 1.0), cutoff_factor::Float64=0.5)

    #Do one phasing
    phased_density = sRAAR(intensity, iterations, beta0, beta_max, tau, cutoff_factor=cutoff_factor)

    #Get result as cube and all mirrored versions
    phased_densityCube = center_cube(getCube(phased_density))
    cubelist = [getAllEvenRotations(phased_densityCube); getAllEvenRotations(mirrorCube(phased_densityCube))]

    #Sort results by best similarity
    cube_cor_list = map((x)->(x,similarity(x,referenceDensityCube)),cubelist)
    sort!(cube_cor_list, lt=(a,b)->a[2]>b[2])

    #Choose best result and return
    best_cube,bestc = cube_cor_list[1]

    if outputfile != "" saveCube(best_cube, outputfile) end
    println("Correlation in realspace of fit: $bestc")#($([cube_cor_list[i][2] for i=1:length(cube_cor_list)]))

    if fourierCube.cubesize > 1
        best_resolution = calculate_maximum_resolution(FSC(best_cube, fourierCube), dr(fourierCube))
        println("Resolution of fit: ", best_resolution)
    end

    return best_cube
end

"""Calculates the FSC from a set of coefficients
Returning: the fitted density cube and the triple of ISC, FSC, ISC_nofitting"""
function calculateSC(volume::SphericalHarmonicsVolume, density::CubeVolume, fourier::CubeVolume, intensity::CubeVolume, num_phases::Int64=1, iterations::Int64=1001, beta0::Float64=0.75, beta1::Float64=0.90, tau::Float64=350.0)

    repetitions = nworkers() > 1 ? (nworkers() < num_phases ? num_phases : nworkers()) : num_phases
    #phasing tries on all available cores
    densities = pmap((vol)->sRAAR_bestfit(density, vol, iterations, beta0, beta1, tau, outputfile="", fourierCube=fourier), collect(Base.Iterators.repeated(volume,repetitions)))

    #Lets calculate the average density
    average_density = reduce(+, densities)

    #Calculate correlations for this result
    isc,fsc,isc_nofitting = ISC(average_density, intensity), FSC(average_density, fourier), shell_correlation_ISC(getCube(volume), intensity)

    println("Best resolution: ", calculate_maximum_resolution(fsc, dr(fourier)))
    println("Best density correlation: ", similarity(average_density,density))

    return average_density,isc,fsc,isc_nofitting
end

"""Performs a postprocessing on a run, including rotational fit in Fourier space and shell correlation calculations"""
function postprocess_run(params, state, reference_pdb_path::String, saving=false, num_points_sphere::Int64=35, sigma::Float64=0.0)

    #Load reference structures as spherical harmonics
    density,fourier,intensity = createSphericalHarmonicsStructure(reference_pdb_path, params["LMAX"], maximum(params["K2_range"]), qmax(maximum(params["K2_range"]), params["qmax"]))
    #Load reference structures as cube
    densCube,fourierCube,intensityCube = createCubicStructure(reference_pdb_path, 2*maximum(params["K2_range"])+1, qmax(maximum(params["K2_range"]), params["qmax"]))

    #intensity from runs
    input_intensity = deleteTerms(state["intensity"],maximum(params["K2_range"]),params["L"])

    #First check if we denoise?
    if sigma > 0.0
        input_intensity = denoise_structure(input_intensity, intensity, maximum(params["K2_range"]), sigma)[1]
    end

    if saving saveCube(input_intensity, "unfitted_intensity.mrc") end

    #Rotational fit in Fourier space
    # state["fittedIntensity"], bestsc, bestt, bestp, bestg = fitStructures(input_intensity, intensity, num_points_sphere, params["K"],params["L"], 0.0, float(pi), 0.0, 2.0*pi, 0.0, 2.0*pi)

    state["fittedIntensity"], bestsc, bestt, bestp, bestg = fitStructures_random(input_intensity, intensity, maximum(params["K2_range"]),params["L"] , 0.99)

    @everywhere gc()

    #Save the fittet intensity
    if saving saveCube(state["fittedIntensity"], "intensity.mrc") end

    #Calculate SC
    fittedDensityCube,isc,fsc,isc_nofitting = calculateSC(state["fittedIntensity"], densCube, fourierCube, intensityCube, 2*nworkers())

    #We are done
    state["sc"] = (isc,fsc,isc_nofitting)
    state["resolution"] = calculate_maximum_resolution(fsc, dr(intensityCube))
    state["state"] = "finished"

    #Save the shell correlations and the phased density
    if saving
        saveState(params, state)
        saveCube(fittedDensityCube,"density.mrc")
    end

    println("Done postprocessing.")
end

##########################################################################################

"""Calculates the corresponding q range for a dq and number of shells K"""
function get_qrange(K::Int64, dq::Float64)
    return collect(1:K) * dq
end

"Given a single fsc curve, calculates k_max for which 0.5=FSC(k_max)"
function calculate_cutoff(sc_curve::Array{Float64}, cutoff::Float64=0.5)
    for k = 1:length(sc_curve)-1
        x1 = k
        x2 = k+1
        y1 = sc_curve[x1]
        y2 = sc_curve[x2]
        if y1 > cutoff && y2 <= cutoff
            m = (y1-y2)/(x1-x2)
            b = y2-m*x2
            kcrossing = (cutoff-b)/m
            return kcrossing
        end
    end
    return length(sc_curve)
end

"Given a single fsc curve, calculates the resolution for which 0.5=FSC(k_max)"
function calculate_maximum_resolution(fsc::Array{Float64}, dq::Float64, cutoff::Float64=0.5)
    return 2.0*pi / ( calculate_cutoff(fsc, cutoff) * dq)
end

"""Alternative approach for resolution calculation"""
function calculate_maximum_resolution(fsc_list::Array{Array{Float64,1},1}, dq::Float64)
    if length(fsc_list) == 0 return 0.0, 0.0 end
    #Calculate the crossing of the mean
    resolutions = Float64[calculate_maximum_resolution(fsc, dq) for fsc in fsc_list]
    return mean(resolutions), std(resolutions)/sqrt(length(fsc_list))
end

"""Given a run in `dir`, returns the isc, fsc, isc_nofitting"""
function load_runs(dir::String, densityCube::CubeVolume, fourierCube::CubeVolume, intensityCube::CubeVolume; range=1000:1019)
    list = Any[]
    energylist = Float64[]
    density_list = CubeVolume[]
    intensity_list = CubeVolume[]
    for i in range
        try
            state = deserializeFromFile("$dir/$i/state.dat")
            push!(list, state["sc"])
            push!(energylist, state["E"])

            density = loadCube("$dir/$i/density.mrc")
            push!(density_list, density)

            intensity = loadCube("$dir/$i/intensity.mrc")
            push!(intensity_list, intensity)
        catch
            println("Couldn't load $dir/$i")
        end
    end

    return analyse_ensemble(density_list, intensity_list, densityCube, fourierCube, intensityCube)
end

# """Takes an ensemble of spherical correlations and calculates max, min, mean and std deviation"""
# function analyse_ensemble(sc_list::Array{Array{Float64}}, K::Int64)
#     if length(sc_list) > 0
#       maximum_sc = Float64[ maximum([length(sc_list[i]) > 1 ? sc_list[i][k] : 0.0 for i = 1:length(sc_list)]) for k=1:K]
#       minimum_sc = Float64[ minimum([length(sc_list[i]) > 1 ? sc_list[i][k] : 0.0 for i = 1:length(sc_list)]) for k=1:K]
#       mean_sc = Float64[ mean([length(sc_list[i]) > 1 ? sc_list[i][k] : 0.0 for i = 1:length(sc_list)]) for k=1:K]
#       stderr_sc = Float64[ std(Float64[length(sc_list[i]) > 1 ? sc_list[i][k] : 0.0 for i = 1:length(sc_list)]) for k=1:K]/sqrt(length(sc_list))
#       return Dict("max"=>maximum_sc,"min"=>minimum_sc,"mean"=>mean_sc,"std"=>stderr_sc)
#     else
#       return Dict("max"=>Float64[], "min"=>Float64[], "mean"=>Float64[], "std"=>Float64[])
#     end
#
# end

"""Takes an ensemble of densities/intensities and calculates the correlations"""
function analyse_ensemble(density_list::Array{CubeVolume}, intensity_list::Array{CubeVolume}, densityCube::CubeVolume, fourierCube::CubeVolume, intensityCube::CubeVolume)

    average_density = reduce(+, density_list)
    average_intensity = reduce(+, intensity_list)

    #Calculate averaged lines
    isc,fsc,isc_nofitting = ISC(average_density, intensityCube), FSC(average_density, fourierCube), shell_correlation_ISC(average_intensity, intensityCube)

    isc_list = map((density)-> ISC(density, intensityCube), density_list)
    fsc_list = map((density)-> FSC(density, fourierCube), density_list)
    isc_nofitting_list = map((intensity)->shell_correlation_ISC(intensity, intensityCube), intensity_list)

    res = calculate_maximum_resolution(fsc, dr(intensityCube))
    _,res_err = calculate_maximum_resolution(fsc_list, dr(intensityCube))

    get_err = (list) -> Float64[ std(Float64[list[i][k] for i = 1:length(list)]) for k=1:length(list[1])]/sqrt(length(list))

    isc_err = get_err(isc_list)
    fsc_err = get_err(fsc_list)
    isc_nofitting_err = get_err(isc_nofitting_list)
    return Dict("isc"=>isc, "isc_err"=>isc_err, "fsc"=>fsc, "fsc_err"=>fsc_err, "isc_nofitting"=>isc_nofitting, "isc_nofitting_err"=>isc_nofitting_err, "res"=>res, "res_err"=>res_err)

end

function load_fitted_intensity(path::String, energy_limit::Float64=1.0e10)
    statepath = "$path/state.dat"
    if isfile(statepath)
        res = deserializeFromFile(statepath)
        if res["E"] < energy_limit && haskey(res, "fittedIntensity")
            return res["fittedIntensity"]
        end
    end
end

function load_density(path::String)
    densitypath = "$path/density.mrc"
    if isfile(densitypath)
        return loadCube(densitypath)
    end
end

#e.g ../parallel/paper_res_vs_L_rotate_all_at_once_Bayes_N32_K16_I0.05_TD0.9999_SF1.01/P163840000_L14
function average_intensities(root::String, output_path::String, densityCube::CubeVolume, fourierCube::CubeVolume, intensityCube::CubeVolume; iterations::Int64=1001, beta0::Float64=0.75, beta1::Float64=0.90, tau::Float64=350.0, energy_limit::Float64=1.0e10, do_phasing::Bool=true)
    try mkpath(output_path) end
    #Average intensity
    intensity_list = [load_fitted_intensity("$root/$num", energy_limit) for num = 1000:1019]
    intensity_list = intensity_list[intensity_list.!=nothing]

    density_list = [load_density("$root/$num") for num = 1000:1019]
    density_list = density_list[density_list.!=nothing]

    average_intensity = getSphericalHarmonicsVolume(reduce(+, map(getSurfaceVolume, intensity_list)))
    average_density = reduce(+, density_list)
    average_density_res = calculate_maximum_resolution(FSC(average_density, fourierCube), dr(fourierCube))

    println("Averaged density has resoution of: $(average_density_res)")
    saveCube(average_intensity, "$(output_path)/intensity.mrc")
    saveCube(average_density, "$(output_path)/averaged_density.mrc")

    average_res = ()
    #Phase averaged intensity
    if do_phasing
        average_res = calculateSC(average_intensity, densityCube, fourierCube, intensityCube, 6, iterations, beta0, beta1, tau)
        saveCube(average_res[1], "$(output_path)/density.mrc")
    end
    serializeToFile("$(output_path)/state.dat", Dict("sc"=>average_res, "fittedIntensity"=>average_intensity))
end

# """Calculates the phasing of the optimal spherical harmonics structure with the given cutoff `K` and `L`"""
function calculate_optimal(intensity::SphericalHarmonicsVolume, densityCube::CubeVolume, fourierCube::CubeVolume, intensityCube::CubeVolume, K::Int64, L::Int64, dir::String, num_phases::Int64=8)
    res = calculateSC(deleteTerms(intensity, K, L), densityCube, fourierCube, intensityCube, num_phases, 1001, 0.75, 0.90, 350.0)
    saveCube(res[1], "$dir/density.mrc")
    saveCube(deleteTerms(intensity, K, L), "$dir/intensity.mrc")
    serializeToFile("$dir/optimal.dat", res)
end

function load_parameter_list(reference_intensity::SphericalHarmonicsVolume)
    parameter_list = Dict[]
    #K2 must be at least such that rmax > 35.0 Angstrom (Phasing)
    deltar = 1.0 #2 Angstrom resolution maximum
    qmax2 = pi/deltar
    K2_min = ceil(Int64, 35.0 / deltar)
    for K2 = K2_min:3:K2_min+12
        deltaq = qmax2 / K2
        for K3::Int64 = 20:2:min(K2,40)
            for L3::Int64 = 12:2:floor(Int64,(K2-1)/2)
                dir = "../parallel/optimal/optimal_K2_$(K2)_K3_$(K3)_L3_$(L3)"
                try
                    optimal = deserializeFromFile("$dir/optimal.dat")
                    retrieved = deserializeFromFile("$dir/retrieved.dat")
                    density = loadCube("$dir/density.mrc")
                    intensity = getCube(deleteTerms(reference_intensity,K3,L3))#loadCube("$dir/intensity.mrc")
                    res_optimal = calculate_maximum_resolution(optimal[3], deltaq)
                    res_retrieved = calculate_maximum_resolution(retrieved[3], deltaq)
                    println("K2=$(K2), K3=$(K3), L3=$(L3), res=$(res_optimal), retrieved=$(res_retrieved)")
                    #Note: K3*(K3+1)*(K3+2)/6*L3^4 is a rough estimation for the total number of operations needed to calculate the model, num_parameters is falsely named at this point
                    push!(parameter_list, Dict("density"=>density, "intensity"=>intensity, "K2"=>K2, "K3"=>K3, "L3"=>L3, "res_optimal"=>res_optimal, "res_retrieved"=>res_retrieved, "num_parameters"=>K3*(K3+1)*(K3+2)/6*L3^4, "optimal"=>optimal, "res_max"=>2*pi/(K3*deltaq)))
                catch
                    println("Failed loading $dir")
                end
            end
        end
    end
    return parameter_list
end

function get_parameter_entry(parameter_list::Array{Dict}, K::Int64, L::Int64)
    filter((x)-> x["K3"] == K && x["L3"] == L,parameter_list)[1]
end

function get_optimal(parameter_list::Array{Dict}, K::Int64, L::Int64)
    return get_parameter_entry(parameter_list, K,L)
end

function load_sc_vs_triplets(dict::Dict, densityCube::CubeVolume, fourierCube::CubeVolume, intensityCube::CubeVolume, K::Int64, L::Int64, dir::String, parameter_list::Array, correlation_list::Dict)
    if !haskey(dict,K) dict[K] = Dict() end
    dict[K][L] = Dict( name => merge(info, load_runs("$dir/$(name)", densityCube, fourierCube, intensityCube; range=1000:1019)) for (name, info) in correlation_list)

    opt = get_optimal(parameter_list, K, L)
    dict[K][L]["optimal"] = merge(
    Dict("shortlabel"=>"optimal", "label"=>"optimal", "linestyle"=>"--", "color"=>"grey"),
    analyse_ensemble([opt["density"]], [opt["intensity"]], densityCube, fourierCube, intensityCube)
    )
end


##################################################################
## noise
#################################################################

function substract_gamma(volume::SphericalHarmonicsVolume, intensity::SphericalHarmonicsVolume, factor::Float64, sigma::Float64, K::Int64)
    new_volume = deepcopy(volume)
    factor *= 1/sigma^3/(2.0*pi)^(3/2)
    for k = 1:K
        setc(new_volume, k, 0, 0, getc(volume,k,0,0) - factor*gaussian_distribution(k*dr(intensity), 0.0, sigma))
    end
    return new_volume
end

function radial_difference_with_intensity(volume::SphericalHarmonicsVolume, intensity::SphericalHarmonicsVolume, K::Int64)
    a = [real(getc(volume,k,0,0)) for k=1:K]
    a /= norm(a)
    b = [real(getc(intensity,k,0,0)) for k=1:K]
    b /= norm(b)
    return norm(a-b)
end

"""Denoises a single structure by making radial part of l=0,m=0 components fit the known radial distribution"""
function denoise_structure(volume::SphericalHarmonicsVolume, intensity::SphericalHarmonicsVolume, K::Int64, sigma::Float64)
    mydiff(x) = radial_difference_with_intensity(substract_gamma(volume, intensity, x, sigma, K), intensity, K)
    result = optimize(mydiff, 0.0, 1.0)
    gammamin = Optim.minimizer(result)
    return substract_gamma(volume, intensity, gammamin, sigma, K), gammamin
end

##################################################################
## Exp. data
#################################################################

export complete_core, phase_completed_intensity

"""Averages intensities and completes core via fitting"""
function complete_core(name::String, c1::C1, center_range::UnitRange{Int64}, range::UnitRange{Int64}; plotting::Bool=true)
    intensities = [deserializeFromFile("$name/$i/intensity.dat") for i = 1000:1019]
    average_intensity_surf = reduce(+, map(getSurfaceVolume, intensities))

    #Filter negativity
    average_intensity_surf.surf = map((x)-> max.(real(x), 0.0), average_intensity_surf.surf)
    average_intensity = getSphericalHarmonicsVolume(average_intensity_surf)
    saveCube(average_intensity, "intensity_averaged.mrc")

    reference_intensity = deepcopy(average_intensity)
    reference_surf = getSurfaceVolume(reference_intensity)
    c1_coliphage = [sum(abs, reference_surf.surf[k]) for k in 1:maximum(range)]

    # curve = deepcopy(c1_coliphage)
    fitting_range=4:10
    curve = log(deepcopy(c1))
    curve[1:3] = zeros(3)
    a = poly_fit(collect(fitting_range), curve[fitting_range], 2)
    func(x,a) = a[1] + x*a[2]+x^2*a[3]
    center_fit = [func(k,a) for k in 1:10]

    corrected_intensity = deepcopy(average_intensity)
    corrected_surf = getSurfaceVolume(corrected_intensity)
    shift = log.(sum(abs, corrected_surf.surf[minimum(range)]))-center_fit[minimum(range)]
    println("shfit = $shift")

    for k in center_range
        corrected_surf.surf[k] = exp(center_fit[k]+shift)*ones(length(corrected_surf.surf[k])) / length(corrected_surf.surf[k])
        # corrected_surf.surf[k] = shift*c1[k]*ones(length(corrected_surf.surf[k])) / length(corrected_surf.surf[k])
    end

    corrected_intensity = getSphericalHarmonicsVolume(corrected_surf)
    c1_coliphage_corrected = [sum(abs, corrected_surf.surf[k]) for k in 1:maximum(range)]
    saveCube(getSphericalHarmonicsVolume(corrected_surf), "intensity_averaged_corrected.mrc")

    if plotting == true
        plot(collect(1:10), curve[1:10], lw=5, label="orig")
        plot(collect(1:10), center_fit, label="fit 1")
        legend()
        savefig("fig1.pdf")

        figure()
        title("Radial sum of coliphage")
        plot(collect(1:maximum(range)),log(c1_coliphage), label="Surf coliphage", lw=3)
        # plot(collect(1:10), center_fit + shift, label="center fit")
        plot(collect(1:maximum(range)), log(c1_coliphage_corrected), label="surf coliphage corr.")
        ylabel("Sum of Radial Part")
        xlim(1,maximum(range))
        # ylim(1.0e-2, 5.0e3)
        xlabel("Shell")
        legend()
        savefig("fig2.pdf")
    end

    # Extend and add zeros at the end
    extended_corrected_intensity = deepcopy(corrected_intensity)
    extend_range = maximum(range)+1:40
    for i in extend_range
        push!(extended_corrected_intensity.coeff, zeros(Complex{Float64}, Base.size(extended_corrected_intensity.coeff[1])))
    end
    extended_corrected_intensity.KMAX += length(extend_range)
    return extended_corrected_intensity
end

function phase_completed_intensity(extended_corrected_intensity::SphericalHarmonicsVolume, num_tries::Int64=8, beta_end::Float64=0.90)

    reference_density_SH = sRAAR(extended_corrected_intensity, 1001, 0.75, beta_end, 350.0, cutoff_factor=0.3)
    reference_density = center_cube(getCube(reference_density_SH))

#     densities = pmap((vol)->sRAAR_bestfit(getCube(extended_corrected_intensity), vol,  1001, 0.75, 0.95, 350.0, cutoff_factor=0.3), collect(Base.Iterators.repeated(extended_corrected_intensity,num_tries-1)))
    densities = pmap((vol)->sRAAR_bestfit(reference_density, vol,  1001, 0.75, beta_end, 350.0, cutoff_factor=0.3), collect(Base.Iterators.repeated(extended_corrected_intensity,num_tries-1)))
    push!(densities, reference_density)
    map!((vol)->center_cube(vol), densities)
    for i = 1:length(densities) saveCube(densities[i], "density_$(i).mrc") end

    averaged_density = reduce(+, densities)

    serializeToFile("density_averaged.dat", averaged_density)

    saveCube(averaged_density, "density_averaged.mrc")
    return averaged_density
end
