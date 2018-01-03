const minstepsize = 1.5*(pi/180.0) #1 degree minimum stepsize for largest L

"""Generates a random start structure for the rotational search"""
function randomStartStructure(volume::SphericalHarmonicsVolume, K::Int64, L::Int64)
    new = deepcopy(volume)
    #Initialize with random unitary matrices
    for l = 2:2:L
        rot = complex(random_rotation(2*l+1))
        ctr = Umat(l) #complex to real transformation
        rtc = ctranspose(ctr) #real to complex transformation
        rot = rtc*rot*ctr
        for k = 1:K
            cvec_set(new,k,l, rot*cvec_get(volume,k,l))
        end
    end
    return new
end

""""Randomly pertubates a given structure by applying a single random rotations to different l's."""
function randomPertubation(volume::SphericalHarmonicsVolume, K::Int64, L::Int64, stepsize=0.1)
    new = deepcopy(volume)
    #Initialize with random unitary matrices
    for l = 2:2:L
        rot = complex(random_rotation_step(2*l+1, stepsize))
        ctr = Umat(l) #complex to real transformation
        rtc = ctranspose(ctr) #real to complex transformation
        rot = rtc*rot*ctr
        for k = 1:K
            cvec_set(new, k,l, rot*cvec_get(volume, k,l))
        end
    end
    return new
end

#--------------------------------------------------------------------------------
# Iterative projections approach
#--------------------------------------------------------------------------------
# function project_C(coeff)
#     newCoeff = deepcopy(coeff)
#     L = Int64(floor((KMAX-1)/2))
#     for l = 2:2:L
#         u,s,v = svd((c[l]*transpose(b[l])) * diagm([1:KMAX;])^2 * coefficients_to_matix(coeff,l)')
#         newCoeffMatrix = b[l]*c[l]*u*v'
#
#         for k = 1:KMAX
#           start = KMAX - (2*l)
#             cvec_set(newCoeff[k], l, newCoeffMatrix[k,start:KMAX])
#         end
#
#     end
#     return real_to_comp(newCoeff)
# end

# function three_photon_walk(coeff, klow = 2, K=10, stepsize=pi/180.0*33.0)
#     coeff = project_C(coeff)
#     bestprob = energy(coeff,ablists, basisindices,basis,K,L)
#     step = deepcopy(coeff)
#     for i = 1:4
# #         surf = map(backtransform, coeff)
# #         surf = map(abs, surf)
# #         coeff = map(transform, surf)
#         for l = 2:2:L
#             rot = complex(random_rotation_step(2*l+1, gaussian(l*stepsize, l*stepsize/2)))
#             rot =  rtc[l]*rot*ctr[l]#rtc*rot*ctr
#             for k = 1:KMAX cvec_set(step[k],l, rot*cvec_get(coeff[k],l)) end
#         end
#         prob = energy(step, K, L)
#         if prob > bestprob
#             bestprob = prob
#             coeff = deepcopy(step)
#         end
#     end
#     # return deleteTerms(coeff,K,L)
#     return coeff
# end

#Projects real electron density shells onto the contraints given by the 2 and 3 photon correlation
# function B2(density,klow=2, K=10)
#     Alm = sphBesselTransform(map(transform,density), true)
#     Almsurf = map(backtransform, Alm)
#
#     #Get amplitudes from projection
#     Ilm = getIntensityCoefficientsIndirect(Alm)
#     pIlm = three_photon_walk(Ilm, klow, K)
#     # pIlm = project_C(Ilm)
#     # pIlm = deleteTerms(pIlm, KMAX, L)
#     amplitudes = map((x) -> sqrt(max(real(backtransform(x)),0)),pIlm)
#
#     #Apply amplitudes
#     for k = 1:KMAX
#         for i = 1:length(Almsurf[k])
#             Almsurf[k][i] = Almsurf[k][i] / abs(Almsurf[k][i]) * amplitudes[k][i]
#         end
#     end
#
#     #go back to real space
#     newrlm = sphBesselTransform(map(transform, Almsurf), false)
#     return map(backtransform, newrlm)
# end

# function projection_refinement(;ablistsfilename="expdata/ablists32.dat", L = 8,klow=2, K = 10, number = 166)
#
#     #loading triplet data
#     loadHistograms(; klow=klow, K=K )
#
#     #precalc rotations
#     global ctr = Dict([l => Umat(l) for l=2:2:L])
#     global rtc = Dict([l => ctranspose(ctr[l]) for l=2:2:L])
#
#     #start with random spherical harmonics coefficients representing a function
#     rlm = getDensityCoefficients( (x) -> x <= rmax/2.0 ? ((p,t) -> 1 + rand()) : ((p,t) -> 0.0) )
#     d = map(backtransform,rlm)
#
#     iterations = 1000
#     beta0 = 0.75
#     beta_max = 0.75
#     tau= 350
#     global Koff = KMAX / 2
#
#   	tmp1 = 2* B2(d, klow, K) - d #--reflection
#     tmp2 = 0.0
#     # diffa = deepcopy(tmp1)
#
#     progress = Progress(iterations, 1, "Structure determination...", 50)
#   	for n = 1:iterations
#       	beta = exp((-n/tau)^3.0)*beta0 + (1-exp((-n/tau)^3.0))*beta_max
#     		tmp3 = A(tmp1)
#     		tmp_u = 0.5*(beta*(2.0*tmp3 -tmp1) + (1-beta)*tmp1 + d)
#     		tmp2 = B2(tmp_u, klow, K)
#
#       	# for k = 1:length(d)
#       	# 	diffa[k] = d[k] - tmp_u[k]
#       	# end
#       	# diff = sumabs(map(sumabs,diffa))/sumabs(map(sumabs,d))
#         diff = 0.0
#
#     		d = deepcopy(tmp_u)
#     		tmp1 = 2*tmp2 - tmp_u
#
#       	progress.desc = "$diff: "
#       	next!(progress)
#   	end
#   	d = deepcopy(tmp2)
#   	print("Phasing complete...")
#     coeff = map(transform, d)
#     saveCube(getCube(coeff, cubesize), cubesize, "parallel/result$(number).mrc", qmax)
# end

#------------------------------------------------------------------------------

"""Main rotation search method
`params` - rotation search parameters
`state` - starting state"""
function rotation_search(params = Dict("reference_pdb_path"=>"crambin.pdb","stepsizefactor"=>1.02, "initial_stepsize" => pi/180.0 * 180.0, "L"=>8, "K3_range" => 1:8, "N"=>32, "histograms"=>"expdata/correlations_N32_K25.dat", "optimizer"=>rotate_all_at_once, "initial_temperature_factor"=>1.0, "measure"=>"Bayes", "temperature_decay"=>0.99, "LMAX"=>25, "K2_range"=>1:35, "qmax"=>1.0, "lambda"=>0.0, "include_negativity"=>false), state = Dict{Any,Any}("newRun"=>true) )

    #Don't start a finished run
    if haskey(state, "state") && state["state"] == "finished_structure"
        return
    end

    #load the histogrammed dataset into the global scope
    c2ref_full, c2ref, c3ref_full, c3ref = loadHistograms(maximum(params["K3_range"]), maximum(params["K3_range"]), params["histograms"], false)

    #Open output file for logging
    out = open("det.out", state["newRun"] ? "w+" : "a")
    println("Starting on machine $(gethostname())")
    flush(STDOUT)

    if state["newRun"]

        #Mark it as "non"-new run for future continuations
        state["newRun"] = false

        #Lets start with the whole range by default (not applied for hierarchical approach)
        state["K"] = maximum(params["K3_range"])

        #Save reference intensity for later use
        if haskey(params, "reference_pdb_path") && params["reference_pdb_path"] != ""
            density,fourier,intensity = createSphericalHarmonicsStructure(params["reference_pdb_path"], params["LMAX"], maximum(params["K2_range"]), qmax(maximum(params["K2_range"]), params["qmax"]))
            state["reference_intensity"] = intensity
        end

        #Retrieve initial structure from 2p-correlation if not provided
        if !haskey(state,"intensity")
            # state["intensity"] = retrieveSolution(c2ref_full/sumabs(c2ref_full), params["L"], params["LMAX"], params["K2_range"], params["qmax"], params["lambda"])


            #TODO: added this to check determination run, becaus retrieveSolution was unstable
            c2_theo = twoPhotons(state["reference_intensity"], BasisType(params["N"], params["L"], params["LMAX"], maximum(params["K3_range"]), params["lambda"], dq(state["reference_intensity"])), maximum(params["K2_range"]), true, false)

            state["intensity"] = retrieveSolution(c2_theo, params["L"], params["LMAX"], params["K2_range"],  params["qmax"], params["lambda"])

            state["intensity"] = randomStartStructure(state["intensity"], state["intensity"].KMAX, state["intensity"].LMAX)
        end

        state["stepsizes"] = Dict()
        state["L"] = 2
        state["i"] = 0
        state["state"] = "running"

        #Make some logging things
        write(out,"""
        #Starting structure determination with the following parameters:
        #$params
        """)

        saveState(params, state)
    end

    try
        CUDA_init()
    catch x
        println(x)
        println("!!! Init of CUDA failed")
        println("!!! Fallback to CPU-backend.")
    end
    flush(STDOUT)

    #begin optimization with given optimizer
    params["optimizer"](out, params, state, c3ref)

    flush(out)
    close(out)

    #Save the final state
    state["state"] = "finished_structure"
    saveState(params, state)
    return params,state
end

"""Creates a Monte Carlo step with given stepsize"""
function get_MC_step(intensity::SphericalHarmonicsVolume, basis::AbstractBasisType, stepsizes::Dict)
    step = deepcopy(intensity)
    for l = 2:2:basis.L
        rot =  basis.rtc[l]*complex(random_rotation_step(2*l+1, gaussian(0.0, stepsizes[l])))*basis.ctr[l]
        for k = 1:intensity.KMAX
            cvec_set(step,k,l, rot*cvec_get(intensity,k,l))
        end
    end
    return step
end

"""Given a particular configuration, calculates the average energy of changes to the structure"""
function calculate_temperature(state::Dict, params::Dict, basis::AbstractBasisType, c3ref::C3, iterations::Int64=50)
    if params["initial_temperature_factor"] > 0.0
        println("Explore energy landscape for temperature... ")
        #Measure average energy variation for 100 steps with the initial_stepsize
        energy_sample = Float64[]

        for i = 1:iterations
            step = get_MC_step(state["intensity"], basis, state["stepsizes"])
            e = energy(step, basis, c3ref, params["K3_range"], params["measure"])
            println("Energy calculation step $i: $e")
            push!(energy_sample,e)
        end
        return std(energy_sample)*params["initial_temperature_factor"]
    else
        return 0.0
    end

end

"""Saves current state and paramters to a file."""
function saveState(params, state)
    serializeToFile("params.dat", params)
    serializeToFile("state.dat", state)
end

"""Regularly do something"""
function regular_action(out, state, params)
    #save out information every 100 iterations
    if mod(state["i"], 500) == 0
        @everywhere gc()
        saveState(params, state)
        flush(out)
        flush(STDOUT)
    end
end

"""Optimizes structure by changing all rotation matrices at once and comparing the structure's three photon correlation with the experimental one"""
function rotate_all_at_once(out, params::Dict, state::Dict, c3ref::C3)

    state["L"] = params["L"]
    d_basis = complexBasis_choice(state["L"], params["LMAX"], params["N"], maximum(params["K3_range"]), params["lambda"], dq(state["intensity"]))

    #Check if this is a continuation or a new run by checking stepsizes
    if !haskey(state["stepsizes"], state["L"])

        #gradually increase K
        state["K"] = maximum(params["K3_range"])

        #set stepsize for this resolution
        #This makes sure that the initial stepsize for L=2 stays as intended and the other L stepsizes are hierarchical higher
        # state["stepsizes"] = Dict(l=>params["initial_stepsize"]*2^(l/2.0-1.0) for l=2:2:state["L"])
        state["stepsizes"] = Dict(l=>params["initial_stepsize"]*(l-1) for l=2:2:state["L"])

        #calculate the initial energy with that resolution
        state["E"] = energy(state["intensity"], d_basis, c3ref, params["K3_range"], params["measure"])

        #Calculate reference energy
        state["reference_energy"] = haskey(state, "reference_intensity") ? energy(state["reference_intensity"], d_basis, c3ref, params["K3_range"], params["measure"]) : 0.0

        println("Reference energy: ", state["reference_energy"])

        #reset the temperature for that resolution
        state["T"] = calculate_temperature(state, params, d_basis, c3ref)

        if params["include_negativity"]
            state["negativity_factor"] = log((state["E"]+state["T"])/state["E"])/0.3 #0.3 is expected initial negativity

            #With new negativity_factor, recalculate energy
            state["E"] = energy(state["intensity"], d_basis, c3ref, params["K3_range"], params["measure"], state["negativity_factor"])
        else
            state["negativity_factor"] = 0.0
        end

        #Let's acknowledge that higher expansion limits L require slower temperature decay
        state["temperature_decay"] = params["temperature_decay"]
    end

    #run until we have converged with the stepsize (stepsize for largest L must below threshold)
    while state["stepsizes"][state["L"]] > minstepsize

        #increase overall step counter
        state["i"] += 1

        #adjust temperature each step with T_new = decay * T_old
        state["T"] *= state["temperature_decay"]

        #Monte Carlo step
        state["step"] = get_MC_step(state["intensity"], d_basis, state["stepsizes"])

        #evaluate step
        @time evaluate_step(out, params, state, d_basis, c3ref)

        regular_action(out, state, params)
    end
end

# "Iterativly optimizing a structure by increasing the angular resolution and adding more shells that are taken into consideration"
# function rotate_hierarchical(out, params::Dict, state::Dict, c3ref::C3)
#
#     for state["L"] = state["L"]:2:params["L"]
#
#         d_basis = complexBasis_choice(state["L"], params["LMAX"], params["N"], maximum(params["K3_range"]), params["lambda"], dq(state["intensity"]))
#
#         #Check if this is a continuation or a new run by checking stepsizes
#         if !haskey(state["stepsizes"], state["L"])
#
#             #gradually increase K
#             state["K"] = maximum(params["K3_range"]) - params["L"] + state["L"]
#
#             #set stepsize for this resolution
#             state["stepsizes"] = Dict(l=>params["initial_stepsize"]/(state["L"]-l+1) for l=2:2:state["L"])
#
#             #calculate the initial energy with that resolution
#             state["E"] = energy(state["intensity"], d_basis, c3ref, params["measure"])
#
#             #Calculate reference energy
#             state["reference_energy"] = haskey(state, "reference_intensity") ? energy(state["reference_intensity"], d_basis, c3ref, params["measure"]) : 0.0
#
#             #reset the temperature for that resolution
#             state["T"] = calculate_temperature(state, params, d_basis, c3ref)
#
#             #Let's acknowledge that higher expansion limits L require slower temperature decay
#             state["temperature_decay"] = params["temperature_decay"]
#         end
#
#         #run until we have converged with the stepsize (stepsize for largest L must below threshold)
#         while state["stepsizes"][state["L"]] > minstepsize
#
#             #increase overall step counter
#             state["i"] += 1
#
#             #adjust temperature each step with T_new = decay * T_old
#             state["T"] *= state["temperature_decay"]
#
#             #Monte Carlo step
#             state["step"] = get_MC_step(state["intensity"], d_basis, state["stepsizes"])
#
#             #evaluate step
#             @time evaluate_step(out, params, state, d_basis, c3ref)
#
#             #perform regular action such as logging
#             regular_action(out, state, params)
#         end
#     end
# end

# "MC sampling based on integration of the Hamilton equations stemming from negative log-probability"
# function rotate_HMC(new, step, out, initial_stepsize, iterations, stepsizecheck, Tstart, L, klow, K, N, E, number)
#
#   #variables for MC annealing
#   # num_accepted_steps = 0
#   # num_MC_steps = 0
#   # stepsize = initial_stepsize
#   dim = 10
#   q = rand(dim)
#   p = rand(dim)
#
#   grad_p = function(p) return p end
#   grad_q = function(q)
#     return q
#   end
#
#   # global coeff = deepcopy(new)
#   # global L_optim = L
#   # numangles = sum([(2*l+1)*l for l=2:2:L])
#   # q = rand(numangles)
#   # p = rand(numangles)
#   # eps = initial_stepsize
#   # E = opt_energy(q)
#   #
#   # grad_p = function(p) return p end
#   # grad_q = function(q)
#   #     derivative = deepcopy(q)
#   #     opt_energy_derivative!(q, derivative)
#   #     return derivative
#   # end
#
#   propose_step = function(q,p, iterations)
#     p = p - 0.5 * eps * grad_q(q)
#
#     for j = 1:iterations
#         q = q + eps * grad_p(p)
#         p = p - eps * grad_q(q)
#     end
#
#     q = q + eps * grad_p(p)
#     p = p - 0.5 * eps * grad_q(q)
#     return q,p
#   end
#
#   for i = 1:iterations
#
#       # # Monitoring acceptance rate and modify stepsize accordingly
#       # if num_MC_steps == stepsizecheck
#       #     rate = num_accepted_steps/num_MC_steps
#       #     num_MC_steps = 0
#       #     num_accepted_steps = 0
#       #     if rate < 0.25
#       #         stepsize /= 2
#       #     # end
#       #     elseif rate > 0.75
#       #         stepsize *= 2
#       #     end
#       # end
#
#       new_q,new_p = propose_step(q,p,100)
#
#       E_new = opt_energy(new_q)
#       delta_E = E_new - E
#       # T =  Tstart*exp(-i/(iterations/10))
#       ap  = exp(-(delta_E)) #acceptance energy
#       accepted = false
#       sign = 0.0
#       #Check if the step improved the triple correlation energy
#       if E_new < E || ( E_new > E && rand() <  ap)
#
#           sign = E_new < E ? 1.0 : -1.0
#           accepted = true
#
#           #save the rotational information
#           E = E_new
#           q = deepcopy(new_q)
#       end
#
#       println("$i\t$sign\t$E_new\t$E\t$delta_E\t$T\t$ap\t$eps")
#
#       if accepted num_accepted_steps += 1 end
#
#       num_MC_steps += 1
#   end
#   return new
# end

"""Increases the stepsize for all L"""
function increase_stepsize(state::Dict, params::Dict)
    for l in keys(state["stepsizes"])
        state["stepsizes"][l] *= params["stepsizefactor"]
    end
end

"""Decreases the stepsize for all L"""
function decrease_stepsize(state::Dict, params::Dict)
    for l in keys(state["stepsizes"])
        state["stepsizes"][l] /= params["stepsizefactor"]
    end
end

function evaluate_step(out, params::Dict, state::Dict, basis::AbstractBasisType, c3ref::C3)
    accepted = false
    sign = 0.0

    #Calculate the triple correlation and return difference to reference triple correlation
    E_new = energy(state["step"], basis, c3ref, params["K3_range"], params["measure"], state["negativity_factor"])
    delta_E = E_new - state["E"]
    ap  = exp(-(delta_E)/state["T"]) #acceptance energy

    #Check if the step improved the triple correlation energy
    if delta_E < 0.0 || ( delta_E > 0.0 && rand() <  ap)
        sign = delta_E < 0.0 ? 1.0 : -1.0
        increase_stepsize(state, params)
        #save the rotational information
        state["E"] = E_new
        state["intensity"] = deepcopy(state["step"])
    else
        decrease_stepsize(state, params)
    end

    #Logging current step
    line = "$(state["i"])\t$sign\t$E_new\t$(state["E"])\t$delta_E\t$(state["reference_energy"])\t$(state["E"]-state["reference_energy"])\t$(state["T"])\t$ap\t$(state["stepsizes"])\t$(state["L"])\t$(state["K"])\n"
    write(out, line)
    println(line)
end

#----------------------------------------------------------------------------------------

function checkRotationSearch(reference::SphericalHarmonicsVolume, K::Int64, L::Int64, basis::AbstractBasisType; histogram="../expdata/correlations_N32_K25_P2048000.dat", save_structures::Bool=false)
    c2full,c2,c3full,c3 = loadHistograms(K, K, histogram)
    start = retrieveSolution(c2full, L, reference.LMAX, reference.KMAX, reference.rmax)
    start = randomStartStructure(deleteTerms(start, K, L), K, L)
    checkRotationSearch(start, deleteTerms(reference, K, L), 1:K, L, basis, save_structures=save_structures)
end

"""Tries out the rotational search scheme by comparing given structures with the known original structure
This yields convergence within the bounds of the scheme"""
function checkRotationSearch(start::SphericalHarmonicsVolume, reference::SphericalHarmonicsVolume, K_range::UnitRange{Int64}, L::Int64, basis::AbstractBasisType; iterations=1.0e4, reduce_stepsize=1000, plotting=false, include_negativity=false, energy=(x)->0.0, save_structures::Bool=false)

    new = deepcopy(start)
    step = deepcopy(new)

    FSC_list = Float64[]
    E_list = Float64[]

    FSC = similarity(reference, start, K_range)
    stepsizes = Dict(l=>pi*(l-1) for l=2:2:L)

    E = energy(new)

    for i = 0:round(Int64,iterations)

        if mod(i, reduce_stepsize) == 0
            E = energy(new)
            push!(E_list, E)
            push!(FSC_list, FSC)

            for l in keys(stepsizes)
                stepsizes[l] /= 2.0
            end
            println("$i\tFSC=$FSC\tE=$E")
            if save_structures
                saveCube(new, "output/rotation_check_$i.mrc")
            end
        end

        #Modify the structure a little bit
        step = get_MC_step(new, basis, stepsizes)

        #Calculate FSC of this new structure
        FSC_new = similarity(reference, step, K_range)

        #Check if we accept this step
        if FSC_new > FSC
            FSC = FSC_new
            new = deepcopy(step)
        end
    end

    return new, FSC
end

# """Takes a set of coefficients and applies rotations to it"""
# function rotate_coefficients(coeff::SphericalHarmonicsVolume, L::Int64, angles::Array{Float64})
#   rotated_coeff = deepcopy(coeff)
#
#   #Rotate the coefficients by the given angles
#   for l = 2:2:L
#
#       num_l_angles = (2*l+1)*l
#       #Unpacking the angles from the vector and calculating the rotation matri
#       R = SON_parametrized(2*l+1,angles[ai:ai+num_l_angles-1])
#
#       #Rotating the individual coefficients for all shells
#       for k = 1:KMAX cvec_set(rotated_coeff,k,l, R*cvec_get(coeff, k,l)) end
#   end
#
#   return rotated_coeff
# end

# end
