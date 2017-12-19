export
    C1, C2, C3,
    C2Shared, C3Shared,
    AbstractBasisType,
    BasisType,
    twoPhotons,
    retrieveSolution,
    alpharange,
    calculate_basis,
    FullCorrelation_parallized,
    energy,
    flab,
    calculate_triple_products,
    calculate_triple_products_fast,
    FullC3,
    c3_slice,
    integrate_c3_shell,
    renormalize_correlation,
    symmetrize_correlation,
    add_Gaussian_filter,
    postprocess_correlations

typealias C1 Vector{Float64}
typealias C1Shared SharedArray{Float64, 1}
typealias C2 Array{Float64, 3}
typealias C2Shared SharedArray{Float64, 3}
typealias C3 Array{Float64, 5}
typealias C3Shared SharedArray{Float64, 5}

abstract AbstractBasisType


"""Datatype for precalculated three photon correlation basis function"""
type BasisType <: AbstractBasisType
    wignerlist::Array{Float64,1}
    indices::Array{Int64,2}
    PAcombos::Array{Int64,2}
    B::Array{Float64,2}
    h_P::Array{Float32,2}
    basislen::Int64
    N::Int64
    L::Int64
    LMAX::Int64
    lrange::StepRange
    ctr::Dict
    rtc::Dict
    K::Int64
    lambda::Float64
    dq::Float64

    function BasisType( wignerlist::Array{Float64,1}, indices::Array{Int64,2}, PAcombos::Array{Int64,2}, B::Array{Float64,2}, h_P::Array{Float32,2}, basislen::Int64, N::Int64, L::Int64, LMAX::Int64, lrange::StepRange, ctr::Dict, rtc::Dict, K::Int64, lambda::Float64, dq::Float64)
        new(wignerlist, indices, PAcombos, B, h_P, basislen, N, L, LMAX, lrange, ctr, rtc, K, lambda, dq)
    end

    function BasisType(N::Int64,L::Int64,LMAX::Int64,K::Int64,lambda::Float64,dq::Float64)
        lrange = (0:2:L)
        ctr = Dict(l => Umat(l) for l=lrange)
        rtc = Dict(l => ctranspose(ctr[l]) for l=lrange)

        new(zeros(Float64,1), zeros(Int64,2,2), zeros(Int64,2,2), zeros(Float64,2,2), Array(Float32,2,2), 0, N, L, LMAX, lrange, ctr, rtc, K, lambda, dq)
    end
end



"""Helper function to calculate basis"""
function calculate_combolist(K::Int64, mcombolist)
    PAcombos = Int64[]
    ki::Int64 = 1
    num::Int64 = 0
    for k1 = 1:K
        for k2 = 1:k1
            for k3 = 1:k2
                for (l1,l2,l3,mcombos, jstart) in mcombolist
                    num+=1
                    append!(PAcombos, [k1,k2,k3,ki,l1,l2,l3,jstart,mcombos])
                end
                ki+=1
            end
        end
    end
    return reshape(PAcombos, 9, num)
end

"""Precalculated values for calculating the three-photon and two-photon correlation"""
function calculate_basis(L::Int64, LMAX::Int64, N::Int64, K::Int64, lambda::Float64=0.0, dq::Float64=0.0, forIntensity=true)

    println("Calculating complex basis with N=$N L=$L K=$K (LMAX=$LMAX, lambda=$lambda, dq=$dq).")

    lrange = forIntensity ? (0:2:L) : (0:L)
    basislen = sum(abs(wigner(l1,m1, l2,m2, l3,-m3)) > 100*eps() ? 1 : 0 for l1=lrange for l2=lrange for l3=lrange for m1=-l1:l1 for m2=-l2:l2 for m3=-l3:l3)
    klength = Integer(K*(K+1)*(K+2)/6)
    qlist = Float64[acos(k*lambda*dq/(4*pi)) for k=1:K]

    wignerlist = Array(Float64, basislen)
    indiceslist = Array(Int64, 9, basislen)
    B = SharedArray(Float64, 2*N^2, basislen)
    P = SharedArray(Float64, klength, basislen)
    # B = zeros(Float64, 2*N^2, basislen)
    # P = zeros(Float64, klength, basislen)
    mcombolist = Vector{Int64}[]

    ctr = Dict(l => Umat(l) for l=lrange)
    rtc = Dict(l => ctranspose(ctr[l]) for l=lrange)

    i = 1
    for l1 = lrange
        for l2 = lrange
            for l3 = lrange
                mcombos = 0
                mcombo_istart = i
                for m1 = -l1:l1
                    for m2 = -l2:l2
                        for m3 = -l3:l3
                            w = wigner(l1,m1, l2,m2, l3,-m3)
                            if abs(w) > 100*eps()
                                wignerlist[i] = w
                                indiceslist[:,i] = Int64[seanindex(l1,m1,LMAX), seanindex(l2,m2,LMAX),seanindex(l3,-m3,LMAX), l1,m1, l2,m2, l3,-m3]
                                i+=1
                                mcombos += 1
                            end
                        end
                    end
                end
                if mcombos>0 push!(mcombolist, Int64[l1,l2,l3,mcombos,mcombo_istart]) end
            end
        end
    end

    PAcombos = calculate_combolist(K,mcombolist)

    # @time begin
    tic()
        @sync @parallel for i=1:basislen
            l1,m1,l2,m2,l3,m3 = indiceslist[4:9,i]
            w = wignerlist[i]
            #Here, only the real part of the complex exponential plays a role
            B[:,i] = reshape(Float64[cos(m2*a + m3*b) for a in alpharange(N), b in alpharange_2pi(2*N)], 2*N^2)
            P[:,i] = reshape(Float64[w*sphPlm(l1,m1,qlist[k1]) * sphPlm(l2,m2,qlist[k2]) * sphPlm(l3,m3,qlist[k3]) for k1=1:K for k2=1:k1 for k3=1:k2], klength)
        end
    # end
    toc()
    h_P = convert(Array{Float32}, transpose(sdata(P)))

    println("Calculation complete ($basislen basislen).")

    return BasisType(wignerlist, indiceslist, PAcombos, sdata(B), h_P, basislen, N, L, LMAX, lrange, ctr, rtc, K, lambda, dq)
end

complexBasis_choice = calculate_basis

function calculate_triple_products_fast(intensity::SphericalHarmonicsVolume, basis::BasisType)
    coeff = intensity.coeff
    PAcombos = basis.PAcombos
    indices = basis.indices
    wignerlist = basis.wignerlist
    P = basis.h_P
    PA = SharedArray(Float64, Base.size(P))

    @sync @parallel for i = 1:Base.size(PAcombos)[2]
        k1,k2,k3,ki,l1,l2,l3,jstart,mcombos = PAcombos[:,i]
        ck1,ck2,ck3 = coeff[k1],coeff[k2],coeff[k3]
        As::Float64=0.0
        for n = 1:mcombos
            j = jstart + n - 1
            @inbounds @fastmath As += wignerlist[j]*real(ck1[indices[1, j]]*ck2[indices[2, j]]*ck3[indices[3, j]])
        end
        @inbounds @fastmath PA[jstart:jstart+mcombos-1,ki] = As * P[jstart:jstart+mcombos-1,ki]
    end
    return sdata(PA)
end

"Calculates the three photon correlation from spherical harmonics coefficients in a paralllized way."
function FullCorrelation_parallized(intensity::SphericalHarmonicsVolume, basis::BasisType, minimal::Bool=true, normalize::Bool=false, return_raw::Bool=false)

    PA = calculate_triple_products_fast(intensity, basis)
    res = basis.B*PA

    if return_raw return res end

    t = zeros(Float64, basis.N, 2*basis.N, basis.K, basis.K, basis.K)
    i = 1
    for k1=1:basis.K for k2=1:k1 for k3=1:k2
        t[:,:, k3, k2, k1] = reshape(res[:, i], basis.N, 2*basis.N)
        i += 1
    end end end
    return normalize ? t / sumabs(t) : t
end

"""Central energy calculation function"""
function energy(intensity::SphericalHarmonicsVolume, basis::AbstractBasisType, c3ref::C3, measure::String="Bayes", negativity_factor::Float64=0.0)
    return energy(intensity, basis, c3ref, 1:basis.K, measure, negativity_factor)
end

"""Central energy calculation function"""
function energy(intensity::SphericalHarmonicsVolume, basis::AbstractBasisType, c3ref::C3, K3_range::UnitRange{Int64}, measure::String="Bayes", negativity_factor::Float64=0.0)
    c3 = FullCorrelation_parallized(intensity, basis, true, true, true)
    c3 = max(c3, 1e-30) #Filter out results where a negative c3 value is expected. This happens in particular when starting structures are derived from sparse (histogrammed) photon correlations.

    # average_negativity = negativityCheck(deleteTerms(intensity,basis.K, basis.L))/basis.K
    # neg = exp(negativity_factor * average_negativity)
    neg = 1.0

    res = 0.0
    if measure == "Bayes"
        i = 1
        for k1 in K3_range
            for k2 = minimum(K3_range):k1
                for k3 = minimum(K3_range):k2
                    p = reshape(c3[:, i], basis.N, 2*basis.N)
                    i += 1
                    q = c3ref[:,:,k3,k2,k1]
                    p = p / sumabs(p)
                    q = q / sumabs(q)
                    @fastmath res += -1.0*sum(q.*log(p))*tripletFactor(k1,k2,k3)
                end
            end
        end
    end
    return neg*res
end

"""Calculates the two photon correlation from spherical harmonics coefficients"""
function twoPhotons(volume::SphericalHarmonicsVolume, basis::BasisType, K2::Int64, minimal::Bool=false, normalize::Bool=false)
    c = zeros(Float64,basis.N,K2,K2)
    mdq = dq(volume)
    for k1=1:K2
        for k2=1:(minimal ? k1 : K2)
            slice = zeros(Float64, basis.N)
            for l = basis.lrange
                #NOTE: No need for conj for complex dot products (it's included in the definition)
                fac = dot(cvec_get(volume,k1,l), cvec_get(volume,k2,l))
                slice += fac * Float64[ Plm(l,0,alpha_star(alpha, k1, k2, mdq, basis.lambda)) for alpha = alpharange(basis.N)]
            end
            c[:,k2,k1] = real(slice)
        end
    end
    return normalize ? c/sumabs(c) : 1/(4*pi)*c
end

"""Alpharange of theory"""
function alpharange(N)
    da = pi/N
    linspace(da/2,pi-da/2,N)
end

"""Alpharange of histograms"""
function alpharange_2pi(N)
    da = 2*pi/N
    linspace(da/2,2*pi-da/2,N)
end

"""Retrieves a set of spherical harmonics coefficeints
Modified version to allow for ranges
"""
function retrieveSolution(c2::C2, L::Int64, LMAX::Int64, K2_range::UnitRange{Int64}, qmax::Float64, lambda::Float64)
    N,_,_ = size(c2)
    K2 = length(K2_range)
    K2_high = maximum(K2_range)
    K2_low = minimum(K2_range)
    @assert L <= K2/2

    #Create empty Spherical Harmonics volume
    intensity = SphericalHarmonicsVolume(LMAX, K2_high, qmax)
    println("Extracting solution over K2_range=$(K2_range) with K2=$(K2) and L=$L.")
    mdq = dq(intensity)
    eigenvecs = Dict()
    eigenvals = Dict()
    Gmatrices = Dict()

    for l = 0:2:L
        G = zeros(Float32, K2, K2)
        for k1 in K2_range
            for k2 = K2_low:k1
                slice = c2[:,k2,k1]

                #symmetrize 2 photon correlation if lambda = 0.0
                if lambda == 0.0
                    slice = 0.5*(slice + reverse(slice))
                end
                A = Float64[ (1/(4*pi))*Plm(l,0,alpha_star(alpha, k1, k2, mdq, lambda)) for alpha = alpharange(N), l = 0:2:L]
                fac = A \ slice
                val = fac[round(Int64,l/2)+1]

                G[k2-K2_low+1,k1-K2_low+1] = val
                G[k1-K2_low+1,k2-K2_low+1] = val
            end
        end
        #Diagonalize the matrix
        F = eigfact(Symmetric(G),K2-(2*l+1)+1:K2)
        eigenval, eigenvectors = F[:values], F[:vectors]

        #Calculate the vectors
        eigenvalmatrix = diagm(sqrt(max(0.0, eigenval)))
        eigenvecs[l] = eigenvectors'
        eigenvals[l] = eigenvalmatrix
        Gmatrices[l] = G
        m = eigenvalmatrix*eigenvectors'

        for k in K2_range
            cvec_set(intensity,k,l,(m[:,k-K2_low+1]))
        end
    end

    intensity = real_to_comp(intensity)

    if negativityCheck(intensity) > 0.33
        for k = 1:intensity.KMAX intensity.coeff[k] *= -1.0 end
    end

    return intensity#,eigenvecs,eigenvals,Gmatrices
end

"""Retrieves a set of spherical harmonics coefficeints"""
function retrieveSolution(c2::C2, L::Int64, LMAX::Int64, K2::Int64, qmax::Float64, lambda::Float64)
    return retrieveSolution(c2, L, LMAX, 1:K2, qmax, lambda)
end

#------------------------------------------------------------------------

"Part of the three photon correlation basis functions"
function flab(k1::Int64, k2::Int64, k3::Int64, dq::Float64, lambda::Float64, l1::Int64, l2::Int64, l3::Int64, a::Float64, b::Float64)
    fl = zero(Complex{Float64})
    theta1 = acos(k1*lambda*dq/(4*pi))
    theta2 = acos(k2*lambda*dq/(4*pi))
    theta3 = acos(k3*lambda*dq/(4*pi))
    if !(abs(l1-l2) <= l3 <= (l1+l2)) return 0.0 end

    for m1p=-l1:l1
        for m2p=-l2:l2
            for m3p=-l3:l3

                fac = wigner(l1,m1p, l2,m2p, l3,-m3p)

                if abs(fac) > 100*eps()
                    cont = sphPlm(l1,m1p,theta1) * sphPlm(l2,m2p,theta2) * sphPlm(l3,-m3p,theta3) * exp(1im*(m2p*a - m3p*b))

                    #The following is wrong, the curvature is already in the thetas, no need to adjust the alpha/betas (phew)
                    # cont = sphPlm(l1,m1p,theta1) * sphPlm(l2,m2p,theta2) * sphPlm(l3,-m3p,theta3) * exp(1im*(m2p*alpha_star(a, k1, k2, dq, lambda) - m3p*alpha_star(b, k1, k3, dq, lambda)))

                    fl += fac* cont
                end
            end
        end
    end
    return Float64(real(fl))
end

"Part of the three photon correlation basis functions"
function flab(l1::Int64, l2::Int64, l3::Int64, a::Float64, b::Float64)
    return flab(0, 0, 0, 0.0, 0.0, l1, l2, l3, a, b)
end

function c3_slice(ck1::Vector{Complex{Float64}}, ck2::Vector{Complex{Float64}}, ck3::Vector{Complex{Float64}}, N::Int64, L::Int64, LMAX::Int64, k1::Int64, k2::Int64, k3::Int64, dq::Float64, lambda::Float64)
    slice = zeros(Complex{Float64}, N, 2*N)
    for l1 = 0:2:L
        for l2 = 0:2:L
            for l3 = 0:2:L
                m = Float64[ flab(k1,k2,k3,dq,lambda, l1,l2,l3, a,b) for a=alpharange(N), b=alpharange_2pi(2*N)]
                if norm(m) > 3*eps()
                    for m1 = -l1:l1
                        for m2 = -l2:l2
                            for m3 = -l3:l3
                                w = wigner(l1,m1, l2,m2, l3,-m3)
                                if abs(w) > 3*eps()
                                    slice += m*w * ck1[seanindex(l1,m1,LMAX)] * ck2[seanindex(l2,m2,LMAX)] * ck3[seanindex(l3,-m3,LMAX)]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return real(slice)
end

"""Calculate complete three-photon correlation
This is just for testing because it's a slow calculation."""
function FullC3(volume::SphericalHarmonicsVolume, L::Int64, K::Int64, N::Int64, LMAX::Int64)

    coeff = volume.coeff
    c3 = zeros(Float64, N, 2*N, K, K, K)

    for k1 = 1:K
        ck1 = coeff[k1]
        for k2 = 1:k1
            ck2 = coeff[k2]
            for k3 = 1:k2
                ck3 = coeff[k3]
                c3[:,:, k3,k2,k1] = c3_slice(ck1, ck2, ck3, N, L)
            end
        end
    end
    return c3 #/ sumabs(c3)
end

# """From a two-photon correlation with k1>k2 restriction, calculates the full version."""
# function complete_two_photon_correlation(c2::C2)
#     N,K,_= Base.size(c2)
#     Float64[k1 >= k2 ? c2[a,k2,k1] : c2[a,k1,k2] for a = 1:N, k1=1:K, k2=1:K]
# end
#
# """From a three-photon correlation with k1>k2>k3 restriction, calculates the full version.
# This is used to complete the histogrammed three-photon correlation before integration."""
# function complete_three_photon_correlation(c3::C3)
#     N,_,K,_,_ = Base.size(c3)
#     c3new = deepcopy(c3)
#     for k1 = 1:K
#         for k2 = 1:K
#             for k3 = 1:K
#                 if     k1 >= k3 && k3 >= k2 c3new[:,:,k3, k2, k1] = c3[:,:,k2,k3,k1]'
#
#                 elseif k3 >= k1 && k1 >= k2
#                     b = c3[:,:,k2,k1,k3]
#                     c3new[:,:,k3, k2, k1] = [b[j,mod(j-i,N)+1] for i = 1:N, j=1:N]
#                 elseif k3 >= k2 && k2 >= k1
#                     b = c3[:,:,k1,k2,k3]
#                     c3new[:,:,k3, k2, k1] = [b[mod(j-i,N)+1,j] for i = 1:N, j=1:N]
#
#                 elseif k2 >= k3 && k3 >= k1
#                     b = c3[:,:,k1,k3,k2]
#                     c3new[:,:,k3, k2, k1] = [b[mod(i-j,N)+1,i] for i = 1:N, j=1:N]
#                 elseif k2 >= k1 && k1 >= k3
#                     b = c3[:,:,k3,k1,k2]
#                     c3new[:,:,k3, k2, k1] = [b[i,mod(i-j,N)+1] for i = 1:N, j=1:N]
#                 end
#             end
#         end
#     end
#     return c3new
# end

#----------------------------------------------------------------

function integrateShell_2pc_alt(intensity::SphericalHarmonicsVolume, N::Int64, K::Int64, iterations::Int64=Integer(1e7))
    da = pi/N
    c2 = zeros(N, K, K)
    c2counts = copy(c2)

    rotations = Matrix{Float64}[random_rotation(3) for k = 1: 10000]
    dq = dr(intensity)
    surf = getSurfaceVolume(intensity)

    for i = 1:iterations
        k1 = Int64(rand(1:K))
        k2 = Int64(rand(1:K))

        p1 = k1*dq*rotations[rand(1:length(rotations))]*[0,1,0]
        p2 = k2*dq*rotations[rand(1:length(rotations))]*[0,1,0]

        a = angle_between(p1,p2)
        ai = Int64(mod(floor(Int64, a/da),N)+1)

        @inbounds c2[ai,k2,k1] += real(getVolumeInterpolated(surf, p1))* real(getVolumeInterpolated(surf, p2))
        @inbounds c2counts[ai,k2,k1] += 1.0

    end
    return c2./c2counts
end

function integrateShell_3pc(intensity::SphericalHarmonicsVolume, N::Int64, K::Int64, iterations::Int64=Integer(1e7))
    da = pi/N
    c3 = zeros(N,N, K, K,K)
    c3counts = copy(c3)

    rotations = Matrix{Float64}[random_rotation(3) for k = 1: 10000]
    dq = dr(intensity)
    surf = getSurfaceVolume(intensity)

    for i = 1:iterations
        k1 = Int64(rand(1:K))
        k2 = Int64(rand(1:K))
        k3 = Int64(rand(1:K))

        a = rand()*2.0*pi
        b = rand()*2.0*pi

        rot = rotations[rand(1:length(rotations))]
        p1 = k1*dq*rot*[0,1,0]
        p2 = k2*dq*rot*[cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 0]*[0,1,0]
        p3 = k3*dq*rot*[cos(b) sin(b) 0; -sin(b) cos(b) 0; 0 0 0]*[0,1,0]

        alpha,beta = mod(angle_between(p1,p2), pi),mod(angle_between(p1,p3), pi)
        ai,bi = Int64(mod(floor(Int64, a/da),N)+1),Int64(mod(floor(Int64, b/da),N)+1)

        @inbounds c3[ai,bi,k3,k2,k1] += real(getVolumeInterpolated(surf, p1))* real(getVolumeInterpolated(surf, p2))* real(getVolumeInterpolated(surf, p3))
        @inbounds c3counts[ai,bi,k3,k2,k1] += 1.0
    end
    #Removing the cout normalization doesn't change anything
    return Float64[c3counts[a,b,k3,k2,k1] > 0 ? c3[a,b,k3,k2,k1] /c3counts[a,b,k3,k2,k1] : 0.0  for a=1:N,b=1:N,k3=1:K,k2=1:K,k1=1:K]
end

function integrateShell_2pc(intensity::SphericalHarmonicsVolume, N::Int64, K::Int64, iterations::Int64=Integer(1e7))
    da = pi/N
    c2 = zeros(N, K, K)
    c2counts = copy(c2)

    rotations = Matrix{Float64}[random_rotation(3) for k = 1: 10000]
    dq = dr(intensity)
    surf = getSurfaceVolume(intensity)

    for i = 1:iterations
        k1 = Int64(rand(1:K))
        k2 = Int64(rand(1:K))

        a = rand()*2.0*pi
        rot = rotations[rand(1:length(rotations))]
        p1 = k1*dq*rot*[0,1,0]
        p2 = k2*dq*rot*[cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 0]*[0,1,0]

        ai = Int64(mod(floor(Int64, a/da),N)+1)

        @inbounds c2[ai,k2,k1] += real(getVolumeInterpolated(surf, p1))* real(getVolumeInterpolated(surf, p2))
        @inbounds c2counts[ai,k2,k1] += 1.0
    end
    return c2./c2counts
end

function integrateShell_3pc_alt(intensity::SphericalHarmonicsVolume, N::Int64, K::Int64, iterations::Int64=Integer(1e7))
    rotations = Matrix{Float64}[random_rotation(3) for k = 1: 300]
    dq = dr(intensity)
    surf = getSurfaceVolume(intensity)
    da = pi/N
    c3 = zeros(N,N, K, K, K)
    #     c3counts = copy(c3)

    up = Float64[0.0, 1.0, 0.0]

    for k1::Int64=1:1K
        p1 = k1*dq*up
        for k2::Int64=1:K
            for k3::Int64=1:K

                for ai::Int64 = 1:N
                    #                     ai = Int64(mod(floor(Int64, a/da),N)+1)
                    a = Float64(da*ai-da/2)
                    p2 = k2*dq*[cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 0]*up

                    for bi::Int64 = 1:N
                        #                         bi = Int64(mod(floor(Int64, b/da),N)+1)
                        b = Float64(da*bi-da/2)
                        p3 = k3*dq*[cos(b) sin(b) 0; -sin(b) cos(b) 0; 0 0 0]*up

                        sum = 0.0
                        for k = 1:length(rotations)
                            rot = rotations[k]
                            sum += real(getVolumeInterpolated(surf, rot*p1))* real(getVolumeInterpolated(surf, rot*p2))* real(getVolumeInterpolated(surf, rot*p3))
                        end
                        c3[ai, bi, k3, k2, k1] = sum
                        #                         c3counts[ai, bi, k3, k2, k1]
                    end
                end
            end
        end
    end
    return c3 / sumabs(c3)
end

"""Integrate out a single shell for comparison with theory"""
function integrate_c3_shell(intensity::SphericalHarmonicsVolume, k1::Int64, k2::Int64, k3::Int64, N::Int64, lambda::Float64, iterations::Int64=Integer(1e6))
    da = pi/N
    mdq = dq(intensity)
    surf = getSurfaceVolume(intensity)

    c3,c3counts = @sync @parallel ( (a,b) -> (a[1]+b[1], a[2]+b[2])) for i = 1:nworkers()
        c3_local = zeros(Float64,N,2*N)
        c3counts_local = zeros(Float64,N,2*N)
        rotations = Matrix{Float64}[random_rotation(3) for k = 1: 10000]

        for j = 1:iterations

            a = rand()*2.0*pi
            b = rand()*2.0*pi

            rot = rotations[rand(1:length(rotations))]
            p1 = k1*mdq*[0,1,0]
            p2 = k2*mdq*[cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 0]*[0,1,0]
            p3 = k3*mdq*[cos(b) sin(b) 0; -sin(b) cos(b) 0; 0 0 0]*[0,1,0]

            alpha,beta = angle_between(p1,p2),angle_between(p1,p3)
            if alpha > pi
                alpha = 2*pi - alpha
                beta = 2*pi - beta
            end

            ai,bi = Int64(mod(floor(Int64, alpha/da),N)+1),Int64(mod(floor(Int64, beta/da),2*N)+1)

            p1e = rot*detector_to_Ewald_sphere(p1, lambda)
            p2e = rot*detector_to_Ewald_sphere(p2, lambda)
            p3e = rot*detector_to_Ewald_sphere(p3, lambda)

            @inbounds c3_local[ai,bi] += real(getVolumeInterpolated(surf, p1e))* real(getVolumeInterpolated(surf, p2e))* real(getVolumeInterpolated(surf, p3e))
            @inbounds c3counts_local[ai,bi] += 1.0

        end
        (c3_local,c3counts_local)
    end
    #Removing the cout normalization doesn't change anything
    return sdata(c3)./sdata(c3counts)
end

#----------------------------

function renormalize_correlation(c2::C2)
    N,K2,_ = Base.size(c2)
    c2_sym = deepcopy(c2)
    for k1=1:K2
        for k2=1:K2
            c2_sym[:,k2,k1] = sumabs(c2[:,k2,k1]) > eps() ? c2[:,k2,k1]/sumabs(c2[:,k2,k1]) : ones(N)
        end
    end
    return c2_sym
end
function renormalize_correlation(c3::C3)
    _,_,K3,_,_ = Base.size(c3)
   c3_sym = deepcopy(c3)
    for k1=1:K3
        for k2=1:K3
            for k3=1:K3
                c3_sym[:,:,k3,k2,k1] = c3[:,:,k3,k2,k1]/sumabs(c3[:,:,k3,k2,k1])# + 1.0
            end
        end
    end
    return c3_sym
end

function symmetrize_correlation(c2::C2)
    _,K2,_ = Base.size(c2)
    c2_sym = deepcopy(c2)
    for k1=1:K2
        for k2=1:K2
            c2_sym[:,k2,k1] = 0.5*(c2[:,k2,k1] + reverse(c2[:,k2,k1]))
        end
    end
    return c2_sym
end

function symmetrize_correlation(c3::C3)
    _,_,K3,_,_ = Base.size(c3)
    c3_sym = deepcopy(c3)
    for k1=1:K3
        for k2=1:K3
            for k3=1:K3
                c3_sym[:,:,k3,k2,k1] = 0.5*(c3[:,:,k3,k2,k1] + rot180(c3[:,:,k3,k2,k1]))
            end
        end
    end
    return c3_sym
end

function add_Gaussian_filter(c2::C2, sigma::Float64=1.0)
    filtered_c2 = deepcopy(c2)
    _,K2,_ = Base.size(c2)
    for k1=1:K2
        for k2=1:K2
            filtered_c2[:,k2,k1] = gaussian_filter(c2[:,k2,k1], sigma)
        end
    end
    return filtered_c2
end

function add_Gaussian_filter(c3::C3, sigma::Float64=1.0)
    c3_filtered = deepcopy(c3)
    _,_,K3,_,_ = Base.size(c3)
    for k1=1:K3
        for k2=1:k1
            for k3=1:k2
                c3_filtered[:,:,k3,k2,k1] = gaussian_filter(c3_filtered[:,:,k3,k2,k1], 1.0)
            end
        end
    end
    return c3_filtered
end

"""Corrects the two-photon and three-photon correlation with the beamstop and adds a Gauss filter to compensate for pixel noise"""
function postprocess_correlations(c2::C2, c3::C3, c2_beamstop::C2, c3_beamstop::C3)
    #Renormalize c2
    c2_corrected = c2 ./ renormalize_correlation(c2_beamstop)

    #Renormalize c3
    c3_corrected = c3 ./ renormalize_correlation(c3_beamstop)

    #Filter c2
    c2_corrected_filtered = add_Gaussian_filter(symmetrize_correlation(c2_corrected), 1.0)

    #Filter c3
    c3_corrected_filtered = add_Gaussian_filter(symmetrize_correlation(c3_corrected), 1.0)
    return c2_corrected_filtered, c3_corrected_filtered
end
