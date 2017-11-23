##########################################################################
# Spherical Harmonics phase retrieval
##########################################################################

"Contraint in real space for spherical harmonics RAAR"
function A(density::SurfaceVolume, Koff::Integer)
	newdensity = deepcopy(density)

	#enforce support and positivity
	for k = 1:density.KMAX
		newdensity.surf[k] = k > Koff ?  zeros(num_val(density.LMAX)) : max(real(density.surf[k]), 0.0)
	end
	return newdensity
end

"applies projections in Fourier space"
function B(density::SurfaceVolume, amplitudes::SurfaceVolume)
	# Fourier Transformation
	fourier = forward(density)

	# Extract the phases from the Fourier transformation
	for k = 1:fourier.KMAX
		fourier.surf[k] = amplitudes.surf[k] .* exp(1im*angle(fourier.surf[k]))
	end

	# Inverse Fourier Transform
	return backward(fourier)
end

"Spherical Harmonics RAAR implementation"
function sRAAR(intensity::SphericalHarmonicsVolume, iterations::Integer, beta0::Float64 = 0.75, beta_max::Float64 = 0.75, tau::Float64 = 100.0, plotting::Bool=false; cutoff_factor::Float64=0.5)

	Koff = ceil(Int64, cutoff_factor*intensity.KMAX)

	amplitudes = getSurfaceVolume(intensity)
	amplitudes.surf = map(sqrt,  map((x)-> max(x, 0.0),  real(amplitudes.surf)) )

	# #Calculating the phases of a ball
	# ball = deepcopy(intensity)
	# for k = 1:ball.KMAX
	#     ball.coeff[k] = zeros(Complex{Float64}, Base.size(ball.coeff[k]))
	#     if k <= Integer(ball.KMAX/2) setc(ball, k, 0, 0, 1.0) end
	# end
	# ball_phases = getSurfaceVolume(forward(ball))

	# Setting up the starting amplitudes with random initial phases
	fourier = deepcopy(amplitudes)
	for k = 1:fourier.KMAX
		phases = pi * 2 * rand(num_val(fourier.LMAX))
		# phases = angle(ball_phases.surf[k])
		fourier.surf[k] = amplitudes.surf[k] .* exp(1im * phases)
	end

	density = backward(fourier)
	tmp1 = 2.0*B(density, amplitudes) - density #reflection
	diffa = deepcopy(tmp1)
	tmp2 = deepcopy(tmp1)

	# progress = Progress(iterations, 1, "Phase retrieval...", 50)
	diff_dens_list = Float64[]
	diff_amp_list = Float64[]

	for n = 1:iterations
		#Alternate reflections
		beta = exp((-n/tau)^3.0)*beta0 + (1.0-exp((-n/tau)^3.0))*beta_max
		tmp3 = A(tmp1, Koff)
		tmp_u = 0.5*(beta*(2.0*tmp3 -tmp1) + (1.0-beta)*tmp1 + density)
		tmp2 = B(tmp_u, amplitudes)

		#Difference tracking in real space
		diff_dens = sumabs(density - tmp_u)/sumabs(density)
		push!(diff_dens_list, diff_dens)

		diff_amp = 0.0
		#Difference tracking in Fourier space
		# d_fourier = forward(tmp_u)
		# diff_amp = sumabs([sumabs(amplitudes.surf[k]) > 1e-16 ? sumabs( abs(d_fourier.surf[k]) - amplitudes.surf[k])/sumabs(amplitudes.surf[k]) : 0.0 for k = 1:density.KMAX])
		# push!(diff_amp_list, diff_amp)

		#Last step
		density = deepcopy(tmp_u)
		tmp1 = 2.0*tmp2 - tmp_u

		#show progress
		# progress.desc = "$(diff_dens) ($(diff_amp)): "
		# next!(progress)
	end
	if plotting
		plot(collect(1:iterations), diff_dens_list, label="diff dens")
		plot(collect(1:iterations), diff_amp_list, label="diff amp")
	end

	density = A(tmp2, Koff)
	return getSphericalHarmonicsVolume(density)
end

##########################################################################
# Cubic phase retrieval
##########################################################################

"Creates a cube of 1.0s and 0.0s with the non-zero region in the center of the cube"
function createCenterCube(cubesize, ratio)
	core = ceil(cubesize/ratio)
	edge = floor( (cubesize - core)/2 )
	center = cubesize/2
	println("Creating support cube with cubesize:$cubesize,ratio:$ratio,edge:$edge,core:$core")
	c = [x > edge && x <= cubesize-edge &&
	y > edge && y <= cubesize-edge &&
	z > edge && z <= cubesize-edge ? 1.0 : 0.0 for x=1:cubesize, y=1:cubesize, z=1:cubesize]
	# c = [ norm([x,y,z] - [center, center, center]) < cubesize/4 ? 1.0 : 0.0 for x=1:cubesize, y=1:cubesize, z=1:cubesize  ]
	# saveCube(c, cubesize, "supportCube.mrc")
	return c
end
#
#
# "Creates a cube with a gaussian kernel of defined sigma"
# function gaussian_kernel(size, sigma)
#     sigma = max(sigma, 0.0)
#     local center = Float64[ceil(size/2),ceil(size/2),ceil(size/2)]
#     kern = [exp(-norm([x,y,z]-center)^2/sigma) for x=1:size, y=1:size, z=1:size]
# end
#
# "Applies a gaussian kernel to a data cube via FFT"
# function gaussian_filter(data, kernel)
#     fkern = fft(ifftshift(kernel))
#     fdata = fft((data))
#     return real((ifft(fkern.*fdata)))
# #     return ones(10,10,10)
# end
#
# "Symmetrizes the cube and multiplies a factor 1.2 to the center voxel as preparation for phase retrieval (this is emprical)"
# function symmetrize_cube(cube)
# 	size = Base.size(cube)[1]
# 	mid = ceil(Int,cubesize/2)
# 	newcube = [(cube[x,y,z] + cube[(size-x)+1,(size-y)+1,(size-z)+1])*0.5 for x = 1:size, y=1:size, z=1:size]
# 	newcube[mid,mid,mid] = newcube[mid,mid,mid]*1.15
# 	return newcube
# end
#
# "Cubic RAAR phase retrieval"
# function craar(coeff,iter=501, b1 = 0.85, b2 = 0.99, tau = 350, K=KMAX, cubesize=cubesize, outputfile="density.mrc")
# 	supportCube = createCenterCube(cubesize, 2)
#
# 	local newIntensCube
# 	if length(coeff) == KMAX
# 		newIntensCube = getCube(coeff, cubesize, K)
# 	else
# 		newIntensCube = coeff
# 	end
#
# 	# newIntensCube = symmetrize_cube(newIntensCube)
# 	gkernel = gaussian_kernel(cubesize, 1.0)
# 	# newIntensCube = gaussian_filter(newIntensCube, gkernel)
#
# 	# newIntensCube = deepcopy(intensCube)
# 	# newIntensCube[floor(cubesize/2+1),floor(cubesize/2+1),floor(cubesize/2+1)] = 0.0#intensCube[cubesize/2+1,cubesize/2+1,cubesize/2+1]
# 	# newIntensCube[cubesize/2+1,cubesize/2+1,cubesize/2+1] = reduce(max, newIntensCube)
# 	amplitudes = sqrt(max(real(newIntensCube), 0.0))
# 	intensStart = amplitudes.* exp( 2.0*pi*rand(cubesize,cubesize,cubesize)*1im)
#
# 	amplitudes = ifftshift(amplitudes)
# 	intensStart = ifftshift(intensStart)
# 	supportCube = ifftshift(supportCube)
#
# 	res = RAAR(iter, b1, b2, tau, amplitudes, intensStart,supportCube, gkernel)
# 	res = fftshift(res)
# 	saveCube(res, cubesize, outputfile)
# 	# res = [res[cubesize-i+1,j,k] for i=1:cubesize,j=1:cubesize,k=1:cubesize]
# 	# saveCube(res, cubesize,"res_m.mrc")
# 	# say("Done")
#
# end
#
# "Magnitude projection for RAAR"
# function MagProj(u, M)
#
# 	epsilon = 3e-16;
#
# 	modsq_u = abs(u).^2
# 	denom2 = modsq_u + epsilon
# 	denom = sqrt(denom2)
# 	r_eps = modsq_u ./ denom - M
# 	dr_eps = (denom2 + epsilon)./(denom2 .* denom)
# 	unew = (1- dr_eps.*r_eps) .* u
# 	return unew
# end
#
# # function RA(data, gkernel, support)
# # 	filtered = gaussian_filter(data, gkernel)
# # 	threshold = maximum(filtered) * 0.1
# # 	new_support = map((x) -> x > threshold ? 1.0: 0.0, filtered)
# # 	u = data.*support.*new_support
# # 	return max(u, 0.0)
# # end
#
# "Real space constraints"
# function RA(data, gkernel, support)
# 	u = data.*support
# 	return max(real(u), 0.0)
# end
#
# "Actual RAAR implementation"
# function RAAR(maxit, beta0, beta_max, tau, amplitudes, intens, support, gkernel)
#
# 	# g = 0.0
# 	beta = beta0
# 	u = ifft(intens)
# 	m_proj = MagProj(intens,amplitudes)
# 	newu = ifft(m_proj)
#
# 	tmp1 = real(2*newu - u)
#
# 	@showprogress 1 "Phasing..." 50 for iter=1:maxit
# 		beta = exp((-iter/tau)^3.0)*beta0 + (1-exp((-iter/tau)^3.0))*beta_max
#
# 		r_sp = 2*RA(tmp1, gkernel, support)  - tmp1
# 		tmp_u = 0.5*(beta*r_sp + (1-beta)*tmp1 +u)
# 		u = real(tmp_u)
# 		intens = fft(u)
# 		m_proj = MagProj(intens,amplitudes)
# 		newu = ifft(m_proj)
# 		tmp1 = real(2*newu - u)
# 	end
# 	return RA(real(newu), gkernel, support)
# end
