export
    C1, C2, C3,
    C2Shared, C3Shared,
    AbstractBasisType,
    BasisType,
    twoPhotons,
    retrieveSolution,
    alpharange,
    complexBasis,
    FullCorrelation_parallized,
    energy,
    flab,
    alpharange


typealias C1 Vector{Float64}
typealias C1Shared SharedArray{Float64, 1}
typealias C2 Array{Float64, 3}
typealias C2Shared SharedArray{Float64, 3}
typealias C3 Array{Float64, 5}
typealias C3Shared SharedArray{Float64, 5}

abstract AbstractBasisType

"""Datatype for precalculated three photon correlation basis function"""
type BasisType <: AbstractBasisType
    basis::SharedArray{Float64,2}
    basisindices::SharedArray{Int64,2}
    basislen::Int64
    N::Int64
    L::Int64
    LMAX::Int64
    lrange::StepRange
    ctr::Dict
    rtc::Dict
end

"""Calculates the two photon correlation from spherical harmonics coefficients"""
function twoPhotons(volume::SphericalHarmonicsVolume, basis::BasisType, K::Int64, minimal::Bool=false, normalize::Bool=false)
    c = zeros(Float64,basis.N,K,K)

    for k1=1:K
        for k2=1:(minimal ? k1 : K)
            slice = zeros(Float64, basis.N)
            for l = basis.lrange
                fac = 0
                for m = -l:l
                    fac += getc(volume, k1, l, m) * conj(getc(volume, k2, l, m))
                end
                slice += fac * Float64[ Plm(l,0,alpha) for alpha = alpharange(basis.N)]
            end
            c[:,k2,k1] = real(slice)
        end
    end
    return normalize ? c/sumabs(c) : c
end

"""Calculates the difference between 2 two-photon correlations"""
function diff_c2(c2_a::C2,c2_b::C2)
    sabs = 0
    for k1=1:K
        for k2=1:K
            sabs += sumabs(c2_a[:,k2,k1] - c2_b[:,k2,k1])/sumabs(c2_a[:,k2,k1])
        end
    end
    return sabs/K^2
end

"""Alpharange of theory"""
function alpharange(N)
    da = pi/N
    linspace(da/2,pi-da/2,N)
end

"""Alpharange of histograms"""
function alpharange_alt(N)
    da = pi/N
    linspace(0,pi-da,N)
end

"""Retrieves a set of spherical harmonics coefficeints"""
function retrieveSolution(c2::C2, L::Int64, LMAX::Int64, KMAX::Int64, qmax::Float64)
    N,K,_ = size(c2)
    println("Extracting solution with K=$K and L=$L.")

    A = Float64[ Plm(l,0,alpha) for alpha = alpharange(N), l = 0:L]
    AI = pinv(A)

    #Create empty Spherical Harmonics volume
    intensity = SphericalHarmonicsVolume(LMAX, KMAX, qmax)
    eigenvecs = Dict()
    eigenvals = Dict()
    Gmatrices = Dict()

    for l = 0:2:L
        G = zeros(K, K)
        for k1 = 1:K
            for k2 = 1:k1
                slice = c2[:,k2,k1]

                #Symmetrizing helps for intensities, but is wrong for complex Fourier results
                slice = 0.5*(slice + reverse(slice)) #symmetrize

                fac = AI*slice
                val = fac[l+1]
                G[k1,k2] = val
                G[k2,k1] = val
            end
        end
        #Diagonalize the matrix
        eigenval, eigenvectors = eig(G) #, permute=false, scale=false doesnt do anything here

        #Calculate the vectors
        eigenvalmatrix = diagm(sqrt(max(eigenval, 0.0)))
        eigenvecs[l] = eigenvectors'
        eigenvals[l] = eigenvalmatrix
        Gmatrices[l] = G
        m = eigenvalmatrix*eigenvectors'
        for k = 1:K
            vec = m[K-(2*l):K,k]
            cvec_set(intensity,k,l, vec)
        end
    end

    intensity = real_to_comp(intensity)

    sign = negativityCheck(intensity) > 0.33 ? -1.0 : 1.0
    for k = 1:intensity.KMAX intensity.coeff[k] *= sign end

    return intensity
end

"Part of the three photon correlation basis functions"
function flab(l1::Int64, l2::Int64, l3::Int64, a::Float64, b::Float64)
    fl = 0.0
    ph = pi/2

    if !(abs(l1-l2) <= l3 <= (l1+l2)) return 0.0 end

    dl1 = mod(l1,2) == 0 ? 2 : 1
    dl2 = mod(l2,2) == 0 ? 2 : 1
    dl3 = mod(l3,2) == 0 ? 2 : 1

    for m1p=-l1:dl1:l1
        for m2p=-l2:dl2:l2
            for m3p=-l3:dl3:l3

                fac = wigner(l1,m1p, l2,m2p, l3,-m3p)

                if abs(fac) > 100*eps()
                    cont = Ylm(l1,m1p, 0.0,ph) * Ylm(l2,m2p, a,ph)* Ylm(l3,-m3p, b,ph)
                    fl += fac* cont
                end
            end
        end
    end
    return Float64(real(fl)) #this really is always real
end

"Precalcualtes the basis functions of the three photon correlations"
function complexBasis(L::Int64, N::Int64, LMAX::Int64, forIntensity=true)
    basisindiceslist = Int64[]
    basislist = Float64[]

    lrange = forIntensity ? (0:2:L) : (0:L)
    ctr = Dict(l => Umat(l) for l=lrange)
    rtc = Dict(l => ctranspose(ctr[l]) for l=lrange)

    num = 0
    for l1 = lrange
        for l2 = lrange
            for l3 = lrange

                m = Float64[ flab(l1,l2,l3, a,b) for a=alpharange(N), b=alpharange(N)]
                if norm(m) > 100*eps()

                    for m1 = -l1:l1
                        for m2 = -l2:l2
                            for m3 = -l3:l3

                                w = wigner(l1,m1, l2,m2, l3,-m3)

                                if abs(w) > 100*eps()

                                    append!(basisindiceslist, Int64[seanindex(l1,m1,LMAX), seanindex(l2,m2,LMAX),seanindex(l3,-m3,LMAX), l1,m1, l2,m2, l3,-m3])

                                    append!(basislist, reshape(m*w, N*N))
                                    num += 1
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    basis = SharedArray(Float64,(N*N,num))
    basis[:,:] = reshape(basislist, N*N, num)

    basisindices = SharedArray(Int64, (num,3+6))
    basisindices[:,:] = transpose(reshape(basisindiceslist,3+6, num))

    basislen = Base.size(basisindices,1)
    println("Calculated complex basis with N=$N L=$L (LMAX=$LMAX): $basislen basislen")

    return BasisType(basis, basisindices, basislen, N, L, LMAX, lrange, ctr, rtc)
end

complexBasis_choice = complexBasis

function calculate_coefficient_matrix(intensity::SphericalHarmonicsVolume, basis::BasisType, K::Int64)
    coeff = intensity.coeff
    faclist = zeros(Float64, basis.basislen,Integer(K*(K+1)*(K+2)/6))
    i = 1
    for k1 = 1:K
        ck1 = coeff[k1]
        for k2 = 1:k1
            ck2 = coeff[k2]
            for k3 = 1:k2
                ck3 = coeff[k3]
                faclist[:,i] = Float64[real(ck1[basis.basisindices[i,1]]*ck2[basis.basisindices[i,2]]*ck3[basis.basisindices[i,3]]) for i=1:basis.basislen]
                i = i + 1
            end
        end
    end
    return faclist
end

"Calculates the three photon correlation from spherical harmonics coefficients in a paralllized way."
function FullCorrelation_parallized(intensity::SphericalHarmonicsVolume, basis::BasisType, K::Int64=8, minimal::Bool=true, normalize::Bool=false, return_raw::Bool=false)
    c = SharedArray(Float64,(basis.N,basis.N,K,K,K))

    fac_matrix = calculate_coefficient_matrix(intensity, basis, K)
    res = basis.basis * fac_matrix

    if return_raw return res end

    t = zeros(Float64, basis.N, basis.N, K, K, K)
    i = 1
    for k1=1:K for k2=1:k1 for k3=1:k2
        t[:,:, k3, k2, k1] = reshape(res[:, i], basis.N, basis.N)
        i += 1
    end end end
    return normalize ? t / sumabs(t) : t
end

"""Calculate correlation slices for a list of (k1,k2,k3)"""
function calculateCorrelationSlices!(c::C3Shared, coeff::Array{Array{Complex{Float64}}}, basis::BasisType, k1::Int64, k2::Int64, k3::Int64)
    basisvec = basis.basis

    ck1,ck2,ck3 = coeff[k1],coeff[k2],coeff[k3]
    faclist = Float64[real(ck1[basis.basisindices[i,1]]*ck2[basis.basisindices[i,2]]*ck3[basis.basisindices[i,3]]) for i=1:basis.basislen]

    c[:,:,k3,k2,k1] = reshape(basisvec*faclist,basis.N,basis.N)
end

"Calculates the three photon correlation from spherical harmonics coefficients in a paralllized way."
function FullCorrelation_parallized_piecewise(intensity::SphericalHarmonicsVolume, basis::BasisType, K::Int64=8, minimal::Bool=true, normalize::Bool=false)
    kcombinations = Tuple{Int64,Int64,Int64}[(k1,k2,k3) for k1 = 1:K for k2 = 1:(minimal ? k1 : K) for k3 = 1:(minimal ? k2 : K)]

    c = SharedArray(Float64,(basis.N,basis.N,K,K,K))

    # Here I had a GC problem:
    # https://github.com/JuliaLang/julia/issues/15155
    # https://github.com/JuliaLang/julia/issues/15415
    pmap((combo)-> calculateCorrelationSlices!( c, intensity.coeff, basis, combo[1], combo[2], combo[3]), kcombinations)
    res = sdata(c)
    return normalize ? res / sumabs(res) : res
end

#--------------------------------------------------

"""Central energy calculation function"""
function energy(intensity::SphericalHarmonicsVolume, basis::AbstractBasisType, c3ref::C3, K::Int64, measure::String="Bayes")
    c3 = FullCorrelation_parallized(intensity, basis, K, true, true, true)
    c3 = max(c3, 1e-30) #Filter out results where a negative c3 value is expected. This happens in particular when starting structures are derived from sparse (histogrammed) photon correlations.

    res = 0.0
    if measure == "Bayes"
        i = 1
        for k1 = 1:K
            for k2 = 1:k1
                for k3 = 1:k2
                    p = reshape(c3[:, i], basis.N, basis.N)
                    i += 1
                    q = c3ref[:,:,k3,k2,k1]
                    p = p / sumabs(p)
                    q = q / sumabs(q)
                    @fastmath res += -1.0*sum(q.*log(p))*tripletFactor(k1,k2,k3)
                end
            end
        end
    end
    return res
end

#------------------------------------------------------------------------

"""Calculate complete three-photon correlation
This is just for testing because it's a slow calculation."""
    function FullC3(volume::SphericalHarmonicsVolume, L::Int64, K::Int64, N::Int64, LMAX::Int64)

        coeff = volume.coeff
        c3 = zeros(Float64, N, N, K, K, K)

        for k1 = 1:K
            ck1 = coeff[k1]
            for k2 = 1:k1
                ck2 = coeff[k2]

                for k3 = 1:k2
                    ck3 = coeff[k3]

                    slice = zeros(Complex{Float64}, N, N)
                    for l1 = 0:2:L
                        for l2 = 0:2:L
                            for l3 = 0:2:L
                                m = Float64[ flab(l1,l2,l3, a,b) for a=alpharange(N), b=alpharange(N)]
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
                    c3[:,:, k3,k2,k1] = real(slice)
                end
            end
        end
        return c3 / sumabs(c3)
    end

    """From a two-photon correlation with k1>k2 restriction, calculates the full version."""
    function complete_two_photon_correlation(c2::C2)
        N,K,_= Base.size(c2)
        Float64[k1 >= k2 ? c2[a,k2,k1] : c2[a,k1,k2] for a = 1:N, k1=1:K, k2=1:K]
    end

    """From a three-photon correlation with k1>k2>k3 restriction, calculates the full version.
    This is used to complete the histogrammed three-photon correlation before integration."""
    function complete_three_photon_correlation(c3::C3)
        N,_,K,_,_ = Base.size(c3)
        c3new = deepcopy(c3)
        for k1 = 1:K
            for k2 = 1:K
                for k3 = 1:K
                    if     k1 >= k3 && k3 >= k2 c3new[:,:,k3, k2, k1] = c3[:,:,k2,k3,k1]'

                    elseif k3 >= k1 && k1 >= k2
                        b = c3[:,:,k2,k1,k3]
                        c3new[:,:,k3, k2, k1] = [b[j,mod(j-i,N)+1] for i = 1:N, j=1:N]
                    elseif k3 >= k2 && k2 >= k1
                        b = c3[:,:,k1,k2,k3]
                        c3new[:,:,k3, k2, k1] = [b[mod(j-i,N)+1,j] for i = 1:N, j=1:N]

                    elseif k2 >= k3 && k3 >= k1
                        b = c3[:,:,k1,k3,k2]
                        c3new[:,:,k3, k2, k1] = [b[mod(i-j,N)+1,i] for i = 1:N, j=1:N]
                    elseif k2 >= k1 && k1 >= k3
                        b = c3[:,:,k3,k1,k2]
                        c3new[:,:,k3, k2, k1] = [b[i,mod(i-j,N)+1] for i = 1:N, j=1:N]
                    end
                end
            end
        end
        return c3new
    end

    # #d t / d U_l,m,mp
    # function derivative_slice(coeff,coeff0,k1,k2,k3, l, m, mp)
    #     c = zeros(Complex{Float64}, N*N)

    #     for i=1:Base.size(basisindices,1)
    #         b = basisindices[i,:]
    #         l1,m1,l2,m2,l3,m3
    #         l1,m1,l2,m2,l3,m3 = b[1],b[2],b[3],b[4],b[5],b[6]
    #         # println("$l1 $l2 $l3 $m1 $m2 $m3")
    #         f1 = getc(coeff[k1], l1, m1)
    #         f2 = getc(coeff[k2], l2, m2)
    #         f3 = getc(coeff[k3], l3, -m3)
    #         bas = basis[:,i]

    #         if l == l1 && m == m1
    #             c+= getc(coeff0[k1], l, mp) * f2 * f3 * basis[:,i]
    #         end
    #         if l == l2 && m == m2
    #             c+= getc(coeff0[k2], l, mp) * f1 * f3 * basis[:,i]
    #         end
    #         if l == l3 && m == -m3
    #             c+= getc(coeff0[k3], l, mp) * f1 * f2 * basis[:,i]
    #         end
    #     end
    #     return c
    # end

    # function derivative(coeff, coeff0, l, m, mp, ksize = KMAX)
    # #     c = Complex{Float64}[]
    #     for k1 = 1:ksize
    #         for k2 = 1:ksize
    #             for k3 =1:ksize
    #                 append!(c, derivative_slice(coeff, coeff0, k1, k2, k3, l, m, mp))
    #             end
    #         end
    #     end
    #     return c
    # end

    # function delta_U(l,coeff, coeff0, K)

    #     global tnew = FullCorrelation(coeff, K)
    #     diff = tref - tnew
    #     n = norm(diff)/norm(tref)

    #     d = zeros(Complex{Float64},2*l+1, 2*l+1)

    #     for m = -l:l
    #         for mp = -l:l
    #             dtdU = derivative(coeff, coeff0, l,m,mp, K)
    #             # dtdU = dtdU / dtdU:norm()

    #             delta = - 2.0 * transpose(diff) * dtdU
    #             d[m+l+1,mp+l+1] = delta[1]
    #         end
    #     end

    #     return d,n
    # end

    # function dtest()
    #     K = 7
    #     L = 4
    #     list = Float64[]
    #     complexBasis(L, 16)
    #     global new = deleteTerms(deepcopy(rCoeff), K, L)
    #     # global new = deepcopy(intensCoeff)
    #     global tref = FullCorrelation(intensCoeff, K)
    #     for k = 1:K
    #         setc(new[k], 0, 0, 0.0)
    #     end

    #     global orig = deepcopy(new)

    #     #Initialize with random unitary matrices
    #     for l = 4:2:L
    #         # global rot = convert(Array{Complex{Float64}},eye(2*l+1)*random_rotation_step(2*l+1, pi/180.0 * 1.0))
    #         # println(rot)
    #         global rot = complex(random_rotation(2*l+1))
    #         # rot = complex(eye(2*l+1))
    #         for k = 1:K cvec_set(new[k],l, rot*cvec_get(orig[k],l)) end
    #     end
    #     modCube = getCube(new, cubesize, K, L)
    #     origCube = getCube(orig, cubesize, K, L)
    #     saveCube(modCube, cubesize, "mirror.mrc")
    #     saveCube(origCube, cubesize, "orig.mrc")
    #     a = reshape(modCube,cubesize^3)
    #     b = reshape(origCube,cubesize^3)
    #     println(dot(a,b) / (norm(a)*norm(b)))
    #     # rot = convert_real_to_complex_matrix(rot, L)
    #     # new = real_to_comp(new)
    #     # global nrot = deepcopy(rot)
    #     # global bla
    #     # du, = delta_U(4, new, intensCoeff, K) #/norm(tref)
    #     # du = du / norm(du)
    #     # for i = 1:1e5
    #     #     du,n = delta_U(4,new,intensCoeff, K)
    #     #     du = du/ norm(du)
    #     #     # push!(list,n)

    #     #     nrot = ( real(nrot-10.0*n*du))
    #     #     println("n=$n det=",det(nrot)," eyediff = ", norm(eye(2*4+1)-nrot))
    #     #     for k = 1:K cvec_set(new[k],L, complex(nrot)*cvec_get(intensCoeff[k],L)) end
    #     #     # plot([1:length(list)], list)
    #     # end
    #     # plotCorrelationSlice(FullCorrelation(new,3),1,2,3,3)
    # end

    #----------------------------------------------------------------

    function integrateShell_2pc_alt(intensity::SphericalHarmonicsVolume, N::Int64, K::Int64, iterations::Int64=1e7)
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

    function integrateShell_3pc(intensity::SphericalHarmonicsVolume, N::Int64, K::Int64, iterations::Int64=1e7)
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

    function integrateShell_2pc(intensity::SphericalHarmonicsVolume, N::Int64, K::Int64, iterations::Int64=1e7)
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

    function integrateShell_3pc_alt(intensity::SphericalHarmonicsVolume, N::Int64, K::Int64, iterations::Int64=1e7)
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
