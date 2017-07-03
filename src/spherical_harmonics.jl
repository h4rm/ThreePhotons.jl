abstract Volume

#General type of a volume decomposed into spherical shells
abstract SphericalVolume <: Volume

"""A volume described by shell-wise spherical harmonics expansions"""
type SphericalHarmonicsVolume <: SphericalVolume
  #spherical harmonics coefficients
  coeff::Array{Array{Complex{Float64}}}

  #bandlimit
  LMAX::Int64

  #radial parameters
  KMAX::Int64
  rmax::Float64

  #radial interpolation
  radial_interp::Bool

  function SphericalHarmonicsVolume(coeff, LMAX::Int64, KMAX::Int64, rmax::Float64, radial_interp::Bool=true)
    new(coeff, LMAX, KMAX, rmax, radial_interp)
  end

  function SphericalHarmonicsVolume(LMAX::Int64, KMAX::Int64, rmax::Float64, radial_interp::Bool=true)
    coeff = [zeros(Complex{Float64}, num_coeff(LMAX)) for k=1:KMAX]
    new(coeff, LMAX, KMAX, rmax, radial_interp)
  end
end

"""A volume described by uniform samples of spherical shells"""
type SurfaceVolume <: SphericalVolume
  #spherical harmonics coefficients
  surf::Array{Array{Complex{Float64}}}

  #bandlimit
  LMAX::Int64

  #radial parameters
  KMAX::Int64
  rmax::Float64

  #radial interpolation
  radial_interp::Bool

  function SurfaceVolume(surf, LMAX::Int64, KMAX::Int64, rmax::Float64, radial_interp::Bool=true)
    new(surf, LMAX, KMAX, rmax, radial_interp)
  end

  function SurfaceVolume(LMAX::Int64, KMAX::Int64, rmax::Float64, radial_interp::Bool=true)
    surf = [zeros(Complex{Float64}, num_val(LMAX)) for k=1:KMAX]
    new(surf, LMAX, KMAX, rmax, radial_interp)
  end
end

#############################################################################
# General Spherical Harmonics functions
#############################################################################

"Wigner 3j symbol from the GSL library"
function wigner(l1::Int64, m1::Int64, l2::Int64, m2::Int64, l3::Int64, m3::Int64)
    return ccall((:gsl_sf_coupling_3j, "libgsl"), Float64, (Int32, Int32, Int32, Int32, Int32, Int32), 2*l1,2*l2,2*l3, 2*m1, 2*m2, 2*m3)
end
#Tested with https://www.wolframalpha.com/input/?i=ThreeJSymbol%5B%7B8,8%7D,+%7B8,-7%7D,+%7B8,-1%7D%5D

"Wrapper for the normalized Legendre polynomials"
function sphPlm(l::Int64,m::Int64,theta::Float64)
    if m>= 0
        return ccall((:gsl_sf_legendre_sphPlm, "libgsl"), Float64, (Int32, Int32, Float64), l,m,cos(theta))
    else
        ma = abs(m)
        return result = (-1.0)^ma * ccall((:gsl_sf_legendre_sphPlm, "libgsl"), Float64, (Int32, Int32, Float64), l,ma,cos(theta))
    end
end

"Wrapper for the assocaited legendre polynomials"
function Plm(l::Int64,m::Int64,theta::Float64)
    if m>= 0
        return ccall((:gsl_sf_legendre_Plm, "libgsl"), Float64, (Int32, Int32, Float64), l,m,cos(theta))
    else
        ma = abs(m)
        return result = (-1.0)^ma*Float64(factorial(BigInt(l-ma))/factorial(BigInt(l+ma))) * ccall((:gsl_sf_legendre_Plm, "libgsl"), Float64, (Int32, Int32, Float64), l,ma,cos(theta))
    end
end

"Spherical harmonics functions"
function Ylm(l::Int64,m::Int64,phi::Float64,theta::Float64)
    return sphPlm(l,m,theta) * exp(1im*m*phi)
end

#############################################################################
# Spherical Harmonics Expansion functions
#############################################################################

function get_size(LMAX::Int64) return 2*LMAX end
function num_val(LMAX::Int64) return (2*LMAX)^2 end
function num_coeff(LMAX::Int64) return LMAX^2 end

function dr(volume::SphericalVolume) return volume.rmax/volume.KMAX end
function dq(volume::SphericalVolume) return volume.rmax/volume.KMAX end
function qmax(volume::SphericalVolume) return pi/dr(volume) end
function qmax(KMAX::Int64, rmax::Float64) return pi/(rmax/KMAX) end

function +(a::SurfaceVolume, b::SurfaceVolume) return SurfaceVolume(a.surf + b.surf, a.LMAX, a.KMAX, a.rmax) end
function -(a::SurfaceVolume, b::SurfaceVolume) return SurfaceVolume(a.surf - b.surf, a.LMAX, a.KMAX, a.rmax) end
function *(a::SurfaceVolume, b::Number) return SurfaceVolume(a.surf * b, a.LMAX, a.KMAX, a.rmax) end
function *(a::Number, b::SurfaceVolume) return SurfaceVolume(b.surf * a, b.LMAX, b.KMAX, b.rmax) end
function /(a::SurfaceVolume, b::Number) return SurfaceVolume(a.surf / b, a.LMAX, a.KMAX, a.rmax) end
function sumabs(a::SurfaceVolume) return sumabs(map(sumabs, a.surf)) end
function length(a::SurfaceVolume) return length(a.surf) end

"""Averages the structure of a range of results for a given run."""
function average(list::Array{SphericalHarmonicsVolume})
    return getSphericalHarmonicsVolume(reduce(+, map(getSurfaceVolume, list)) / length(list))
end

"""Disables the radial interpolation of the spherical harmonics volume"""
function disable_radial_interpolation(volume::SphericalHarmonicsVolume)
  volume.radial_interp = false
end

"""calculate the coefficients from surface"""
function transform(surface::Vector{Complex{Float64}}, LMAX::Int64)
    rcoeff = zeros(Float64, num_coeff(LMAX))
    icoeff = zeros(Float64, num_coeff(LMAX))

    ccall( (:transform_sample, "sh/libsh.so"), Ptr{Void},
    (Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}, Int32, Int32),
    real(surface), imag(surface), rcoeff, icoeff, LMAX, 1)

    return complex(rcoeff,icoeff)
end

"""Calculates surface from coefficients"""
function backtransform(coeff::Vector{Complex{Float64}}, LMAX::Int64)
    rdata = zeros(Float64, num_val(LMAX))
    idata = zeros(Float64, num_val(LMAX))

    ccall((:transform_sample, "sh/libsh.so"), Ptr{Void},
    (Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64,}, Int32, Int32),
    real(coeff), imag(coeff), rdata, idata, LMAX, 0)

    return complex(rdata, idata)
end

"""SphericalHarmonicsVolume to SurfaceVolume"""
function getSurfaceVolume(volume::SphericalHarmonicsVolume)
  surf = map((x) -> backtransform(x, volume.LMAX), volume.coeff)
  return SurfaceVolume(surf, volume.LMAX, volume.KMAX, volume.rmax, volume.radial_interp)
end

"""SurfaceVolume to SphericalHarmonicsVolume"""
function getSphericalHarmonicsVolume(volume::SurfaceVolume)
  coeff = map((x) -> transform(x, volume.LMAX), volume.surf)
  return SphericalHarmonicsVolume(coeff, volume.LMAX, volume.KMAX, volume.rmax, volume.radial_interp)
end

"""Helper function to get index in linear-array for (l,m)-doublet"""
function seanindex(l::Int64,m::Int64,LMAX::Int64)
    bigL = LMAX - 1
    if m >= 0
        return round(Int64, ( m * ( bigL + 1 ) - ( ( m * (m - 1) ) /2 ) + ( l - m ) ) +1 )
    else
        return round(Int64,  ( ( ( bigL * ( bigL + 3 ) ) /2 ) + 1 + ( ( bigL + m ) * ( bigL + m + 1 ) / 2 ) + ( l - abs( m ) ) ) + 1 )
    end
end

"""Helper function: gets specific spherical harmonics coefficient!"""
function getc(volume::SphericalHarmonicsVolume,k::Int64,l::Int64,m::Int64)
    return volume.coeff[k][seanindex(l,m,volume.LMAX)]
end
"""Helper function: sets specific spherical harmonics coefficient"""
function setc(volume::SphericalHarmonicsVolume,k::Int64,l::Int64,m::Int64,val::Complex{Float64})
    volume.coeff[k][seanindex(l,m,volume.LMAX)] = val
end

"""Helper function: sets specific spherical harmonics coefficient"""
function setc(volume::SphericalHarmonicsVolume,k::Int64,l::Int64,m::Int64,val::Float64)
    setc(volume, k, l, m, Complex(val))
end

"""Gets the whole 2*l+1 dimensional coefficients vector"""
function cvec_get(volume::SphericalHarmonicsVolume,k,l)
    v = Array(Complex{Float64},2*l+1)
    for m=-l:l
        v[m+l+1] = getc(volume,k,l,m)
    end
    return v
end
"""Sets the whole 2*l+1 dimensional coefficient vector"""
function cvec_set(volume::SphericalHarmonicsVolume,k, l, v)
    for m=-l:l
        setc(volume,k, l, m, v[m+l+1])
    end
end

"""helper function that translates an index to two angles, according to the size
returns PHI, THETA"""
function i_to_ang(index, LMAX)
    size = get_size(LMAX)
    iphi = mod(index-1,size)
    itheta = floor((index-1)/size)
    phi = iphi / (size) * 2.0*pi
    theta = pi*(2*itheta + 1)/(2*size)
    return (phi,theta)
end

"""Calculates index in linear surface array for given angles"""
function ang_to_i(phi::Float64,theta::Float64, LMAX)
    size = get_size(LMAX)
    phi = mod(phi, 2.0*pi)
    theta = mod(theta, pi)
    ip = round(Int64, phi/(2.0*pi)*size)
    it = round(Int64, ((theta*2*size)/pi - 1)/ 2)
    i = clamp(it*size+ip + 1,1,size^2)
    return Int64(i)
end

"""From a surface function, calculate the spherical harmonics coefficients up to LMAX"""
function getSurface(func, LMAX::Int64)
    surface = Array(Complex{Float64}, num_val(LMAX))
    for i = 1:num_val(LMAX)
        phi,theta = i_to_ang(i,LMAX)
        surface[i] = func(phi,theta)
    end
    return surface
end

"""Calculates the SH coefficients for a density"""
function getDensity(atomlist, LMAX::Int64, KMAX::Int64, rmax::Float64)
    surflist = Array{Complex{Float64},1}[]
    dr = rmax/KMAX
    for k = 1:KMAX
        r = k*dr
        push!(surflist, getSurface(densityShell(atomlist, k*dr), LMAX))
    end
    surfVolume = SurfaceVolume(surflist, LMAX, KMAX, rmax)

    return getSphericalHarmonicsVolume(surfVolume)
end

"""wraps the density as a function only of the angle with fixed radius"""
function densityShell(atomlist, r)
   return function(phi,theta) return getDensity_cartesian(atomlist, sphericalToCartesian(r,phi,theta)) end
end

"""Calculates the Fourier shell"""
function fourierShell(atomlist, k)
    return function(phi,theta) return getFourier_cartesian(atomlist,sphericalToCartesian(k,phi,theta)) end
end
"""Gets intensity shell"""
function intensityShell(atomlist, k)
    return function(phi,theta) return getIntensity_cartesian(atomlist,sphericalToCartesian(k,phi,theta)) end
end

# """Calculates the spherical harmonics coefficients from a cubic representation"""
# function cubic_to_spherical(cube)
#     cubic_shell = function(r)
#       return function(phi, theta)
#         x = clamp(round(Int64,sphericalToCartesian(r,phi,theta)) + round(Int64,[(cubesize-1)/2 + 1, (cubesize-1)/2 + 1, (cubesize-1)/2 + 1] ), 1, cubesize)
#         return cube[x[1],x[2],x[3]]
#       end
#     end
#     coeff = []
#     for i = 1:KMAX
#         r = i
#         c = get_coefficients(cubic_shell(r))
#         push!(coeff,c)
#     end
#     return coeff
# end

"""Deletes upper terms of a spherical harmonics representation
deleteTerms(coeff,kuct,L, LMAX)"""
function deleteTerms(volume::SphericalHarmonicsVolume, K::Int64, L::Int64)
    newvol = deepcopy(volume)
    for k = 1:newvol.KMAX
        for l = 0:newvol.LMAX-1
            for m = -l:l
                if k > K || l > L
                    setc(newvol,k,l,m,0.0)
                end
            end
        end
    end
    return newvol
end

# Hankel transform

"""Numerical integration of \$f \\cdot r dr\$ """
function frdr(f::Vector{Float64}, r::Vector{Float64})
    l = length(r)
    res = Array(Float64,l)

    @inbounds res[1] = sqrt(0.5) * (r[1] + r[2])^2
    for j = 2:l-1
         @inbounds res[j] = (r[j+1] - r[j-1]) * (r[j-1] + 2*r[j] + r[j+1])
    end
    @inbounds res[l] = 4*r[l]^2 - (r[l-1] + r[l])^2

    return pi/4 * (res.^(1.5) .*f)
end

"""Calculates the Hankel transform
Source:
# *http://www.mathworks.de/matlabcentral/fileexchange/13371-hankeltransform
# *http://documents.epfl.ch/users/l/le/leuteneg/www/MATLABToolbox/HankelTransform.html"""
function hat!(res::Vector{Float64}, h::Vector{Float64}, I::Matrix{Float64}, forward::Bool=true)
    l = length(h)
    if forward
      r = collect(1.0:l)
      res[:] = I*frdr(h, r)
    else
      k = pi/l * collect(1.0:l)
      res[:] = transpose(I)*frdr(h, k)/(2.0*pi)^2
    end
end

"""Precalcualtes Imatrix for numerical integration of Hankel transform"""
function precalcImatrices(LMAX::Int64, KMAX::Int64)

  Imatrix = function(len, n)
      r = Float64[1:len;]
      k = pi/len * Float64[1:len;]
      return besselj(n, k*transpose(r))
  end

  #Calculate the matrices, if not already available
  if !isdefined(:Ilist) || length(Ilist) != LMAX || size(Ilist[1])[1] != KMAX
    global Ilist = []
    for l = 0:LMAX-1
        push!(Ilist, Imatrix(KMAX, l+0.5))
    end
  end

end

"""On a set of spherical harmonics coefficients, this performs a hankel transform with a given max value of rho (rmax or qmax)"""
function sphBesselTransform(volume::SphericalHarmonicsVolume, forward::Bool)
    newvol = deepcopy(volume)

    #Precalculate the basis
    precalcImatrices(volume.LMAX, volume.KMAX)
    samp_r, samp_c = zeros(Float64, volume.KMAX), zeros(Float64, volume.KMAX)

    for l = 0:volume.LMAX-1
        ipart = Ilist[l+1]
        for m = -l:l
            rho = Complex{Float64}[getc(volume,k,l,m) for k=1:volume.KMAX]
            hat!(samp_r, real(rho), ipart, forward)
            hat!(samp_c, imag(rho), ipart, forward)

            for k=1:volume.KMAX
                fac = forward ? Complex{Float64}((1im)^l *sqrt(1/k)) : Complex{Float64}((-1im)^l *sqrt(1/k))
                setc(newvol, k, l, m, fac*complex(samp_r[k],samp_c[k]))
            end
        end
    end
    newvol.rmax = qmax(volume) #this also gives rmax for backtransform
    return newvol
end

"""Takes density SH and performs hankel transform to get amplitude SH"""
function forward(volume::SphericalHarmonicsVolume)
   return sphBesselTransform(volume, true)
end

"""Takes density SH and performs hankel transform to get amplitude SH"""
function backward(volume::SphericalHarmonicsVolume)
   return sphBesselTransform(volume, false)
end

"Converts a density surface array into its Fourier counterpart"
function forward(volume::SurfaceVolume)
	volume = getSphericalHarmonicsVolume(volume)
	volume = forward(volume)
	return getSurfaceVolume(volume)
end

"Converts a fourier surface array into its real space counterpart"
function backward(volume::SurfaceVolume)
  volume = getSphericalHarmonicsVolume(volume)
	volume = backward(volume)
	return getSurfaceVolume(volume)
end

"""Indirectly calculate intensity coefficients by squaring on surface"""
function absoluteSquare(volume::SphericalHarmonicsVolume)
    newvol = deepcopy(volume)
    newvol.coeff = map((x) -> transform( complex(abs(backtransform(x, volume.LMAX)).^2), volume.LMAX), volume.coeff)
    return newvol
end

"Calculates the Intensity Shell Correlation for a set of coefficients"
function shell_correlation(volume1::SphericalHarmonicsVolume, volume2::SphericalHarmonicsVolume, K::Int64=0, cor_method=cor_FSC)
  surf1 = getSurfaceVolume(volume1)
  surf2 = getSurfaceVolume(volume2)
  @assert length(surf1) == length(surf2)
  K = (K == 0 ? length(surf1) : K)

  return Float64[cor_method(surf1.surf[k], surf2.surf[k]) for k=1:K]
end

function shell_correlation_ISC(volume1::SphericalHarmonicsVolume, volume2::SphericalHarmonicsVolume, K::Int64=0)
  shell_correlation(volume1, volume2, K, cor_ISC)
end

function shell_correlation_FSC(volume1::SphericalHarmonicsVolume, volume2::SphericalHarmonicsVolume, K::Int64=0)
  shell_correlation(volume1, volume2, K, cor_FSC)
end

function ISC(density::SphericalHarmonicsVolume, intensity::SphericalHarmonicsVolume, K::Int64=0)
  return shell_correlation_ISC(intensity, absoluteSquare(forward(density)), K)
end

function FSC(density::SphericalHarmonicsVolume, fourier::SphericalHarmonicsVolume, K::Int64=0)
  return shell_correlation_FSC(fourier, forward(density), K)
end

"Calculates the similarity between two sets of coefficients via Fourier Shell Correlation"
function similarity(volume1::SphericalHarmonicsVolume, volume2::SphericalHarmonicsVolume, K::Int64)
  return sum(shell_correlation_ISC(volume1, volume2, K))/K
end

"Calculates the similarity between two sets of coefficients via Fourier Shell Correlation"
function similarity(volume1::SphericalHarmonicsVolume, volume2::SphericalHarmonicsVolume)
  return similarity(volume1, volume2, volume1.KMAX)
end

#Complex to real and real to complex

"""Calculates the U matrix for order l that translates complex to real coefficients"""
function Umat(l)
    l = float(l)
    Gamma = (m)-> m>0 ? 1 : 0

    Delta = (a,b)-> a==b ? 1 : 0

    U = function(u,m)
        return ( Delta(m,0)*Delta(u,0) + 1/sqrt(2)*( Gamma(u)*Delta(m,u) + Gamma(-u)*1im*(-1)^m*Delta(m,u)+
                Gamma(-u)*(-1im)*Delta(m,-u)+Gamma(u)*(-1)^m*Delta(m,-u)))
    end

    return Complex{Float64}[ U(x,y) for x = -l:l, y= -l:l]
end

"""input array of [k,l] entries with [-m ... m] spherical harmonics coefficients in a vector"""
function comp_to_real(volume::SphericalHarmonicsVolume)
    newvol = deepcopy(volume)
    for l = 1:volume.LMAX-1
        mat = Umat(l)
        for k = 1:volume.KMAX
            cvec_set(newvol, k,l, mat*cvec_get(volume, k, l))
        end
    end
    return newvol
end

"""Calculates real spherical harmonics coefficients from complex"""
function real_to_comp(volume::SphericalHarmonicsVolume)
    newvol = deepcopy(volume)
    for l = 1:volume.LMAX-1
        mat = ctranspose(Umat(l))
        for k = 1:volume.KMAX
            cvec_set(newvol,k,l, mat*cvec_get(volume, k, l))
        end
    end
    return newvol
end

"""Given a set of surfaces, calculates the cartesian value via interpolation"""
function getVolumeInterpolated(volume::SurfaceVolume, vec::Vector{Float64})
    q,p,t = cartesianToSpherical(vec)
    size = get_size(volume.LMAX)
    dp = 2.0*pi / size
    dt = pi / size
    surf = volume.surf
    k = q/dr(volume)

    k0 = Int64(clamp(round(Int64, round(k)), 1, volume.KMAX))
    k1 = k0

    if volume.radial_interp
      k0 = Int64(clamp(round(Int64, floor(k)), 1, volume.KMAX))
      k1 = Int64(clamp(round(Int64, ceil(k)), 1, volume.KMAX))
    end

    #phi_k = 2*pi*k/2B,
    p0 = floor(p/dp)*dp
    p1 = ceil(p/dp)*dp

    #theta_j = (2j + 1)*pi/4B
    # t0 = clamp(floor(t/dt-0.5)*dt,0,pi)
    # t1 = clamp(ceil(t/dt-0.5)*dt,0,pi)
    t0 = max(floor( ((t*2*size)/pi - 1)/ 2) *dt + dt/2, 0)
    t1 = min(ceil( ((t*2*size)/pi - 1)/ 2) *dt+ dt/2, pi)

    a = k-k0
    b = (p-p0)/dp
    c = (t-t0)/dt

    V000 = (1-a) * (1-b) * (1-c)
    V100 = a *(1-b) *(1-c)
    V010 = (1-a)* b* (1-c)
    V001 = (1-a)* (1-b)* c
    V101 = a*(1-b)* c
    V011 = (1-a)* b* c
    V110 = a* b* (1-c)
    V111 = a* b* c

    return      V000 * surf[k0][ang_to_i(p0,t0,volume.LMAX)] +
                V100 * surf[k1][ang_to_i(p0,t0,volume.LMAX)] +
                V010 * surf[k0][ang_to_i(p1,t0,volume.LMAX)] +
                V001 * surf[k0][ang_to_i(p0,t1,volume.LMAX)] +
                V101 * surf[k1][ang_to_i(p0,t1,volume.LMAX)] +
                V011 * surf[k0][ang_to_i(p1,t1,volume.LMAX)] +
                V110 * surf[k1][ang_to_i(p1,t0,volume.LMAX)] +
                V111 * surf[k1][ang_to_i(p1,t1,volume.LMAX)]
end

"""Calculates a cube from spherical shells"""
function getCube(volume::SurfaceVolume)
  cubesize = 2*volume.KMAX+1
  getCube(volume, cubesize)
end

"""Calculates a cube from spherical harmonics coefficients"""
function getCube(volume::SphericalHarmonicsVolume)
  return getCube(getSurfaceVolume(volume))
end

"""Calculates a cube from spherical shells"""
function getCube(volume::SurfaceVolume, cubesize::Int64)
  r = linspace(-volume.rmax, volume.rmax, cubesize)
  return CubeVolume(Complex{Float64}[getVolumeInterpolated(volume,[x,y,z]) for x=r, y=r, z=r], cubesize, volume.rmax)
end

"""Calculates a cube from spherical harmonics coefficients"""
function getCube(volume::SphericalHarmonicsVolume, cubesize::Int64)
  return getCube(getSurfaceVolume(volume), cubesize)
end

"""Saving the cube calculated from a spherical harmonics volume"""
function saveCube(volume::SphericalHarmonicsVolume, filename)
  saveCube(getCube(volume), filename)
end

##########################################################
# rotations
##########################################################


"""Big number factorial function"""
function myfactorial(x)
  return factorial(BigInt(x))
end

"""working calculation of the matrix d elements taken from:
Source: unknown, maybe Ludger?"""
function d_work(l::Int64, m1::Int64, m2::Int64, beta::Float64)
    t1 = BigInt(l+m1)
    t2 = BigInt(l-m2)
    t3 = BigInt(m2 - m1)
    t4 = BigInt(2*l + m1 - m2)
    vmin = (-t3 < 0) ? 0 : -t3
    vmax = (t2 < t1) ? t2 : t1

    factor = (myfactorial(l+m1)*myfactorial(l-m1)*myfactorial(l+m2)*myfactorial(l-m2))^(1/2)
    val = 0

    for s = vmin:vmax
        val += cos(beta/2)^(2*l+m1-m2-2*s) * sin(beta/2)^(m2-m1+s*2) * (-1.0)^(s+m2-m1)/(myfactorial(l+m1-s)*myfactorial(s)*myfactorial(m2-m1+s)*myfactorial(l-m2-s))
    end
    return Complex{Float64}(val*factor)
end

function Dlms(l::Int64, m::Int64, s::Int64, theta::Float64, phi::Float64, gamma::Float64)
    return exp(-1im*s*gamma-1im*m*phi) * d_work(l, m, s, theta)
end

function rotationMatrix(theta::Float64, phi::Float64, gamma::Float64, l::Int64)
    return Complex{Float64}[ Dlms(l, m, s, theta, phi, gamma) for m=-l:l,s=-l:l]
end

"""Rotates a volume via rotating the coefficients"""
function rotateStructure(volume::SphericalHarmonicsVolume, theta::Float64, phi::Float64, gamma::Float64, K::Int64, lrange::Range)
	newvol = deepcopy(volume)

	for l in lrange
		rot = rotationMatrix(theta, phi, gamma, l)
		for k = 1:K
			cvec_set(newvol,k,l, rot*cvec_get(volume, k,l) )
		end
	end
	return newvol
end
"""Overloaded version of the rotateStruture function"""
function rotateStructure(volume::SphericalHarmonicsVolume, theta::Float64, phi::Float64, gamma::Float64, K::Int64, L::Int64)
  return rotateStructure(volume, theta, phi, gamma, K, 1:L)
end

"""Calculates the function at k,phi,theta for a set of spherical harmonic coefficient"""
function get_value(volume::SphericalHarmonicsVolume,k, phi, theta)
    val = 0
    for l = 0:volume.KMAX-1
        for m=-l:l
            val += getc(volume,k,l,m)*(Ylm(l,m,phi, theta))
        end
    end
    return val
end

function negativityCheck(volume::SphericalHarmonicsVolume)
    volume = getSurfaceVolume(volume)
    mapreduce((x) -> max(sumabs(min(real(x),0))/sumabs(x),0), +, volume.surf)
end
