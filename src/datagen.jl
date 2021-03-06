export
    Noise,
    GaussianNoise,
    pointsPerOrientation,
    tripletFactor,
    doubletFactor,
    histogramCorrelationsInPicture_alltoall,
    get_noise_volume,
    compton_scattering,
    compton_scattering_q,
    get_compton_noise_volume,
    generateHistogram,
    loadHistograms,
    countTriplets,
    countDoublets,
    detector_to_Ewald_sphere,
    alpha_star_inverse,
    alpha_star,
    ImageVolume,
    getVolumeInterpolated


"""Abstract noise type, to be extended in the future"""
abstract type Noise end

"""Radially symmetric Gaussian noise with sigma(width) and gamma(height)"""
type GaussianNoise <: Noise
    gamma::Float64
    sigma::Float64
    photons::Int64
end

function detector_to_Ewald_sphere(vec::Vector{Float64}, lambda::Float64=1.0)
    q = norm(vec[1:2])
    theta = acos(q*lambda/(4*pi))
    c = cos(theta) #Note that theta should be pi/2 for lambda -> 0
    s = sin(theta)
    return Float64[s*vec[1:2]; c*q]
end

"""Corrects for Ewald sphere curvature"""
function alpha_star_inverse(alpha_star::Float64, k1::Int64, k2::Int64, dq::Float64, lambda::Float64)
    theta1 = acos(k1*dq*lambda/(4*pi))
    theta2 = acos(k2*dq*lambda/(4*pi))
    return acos((cos(alpha_star) - cos(theta1)*cos(theta2))/(sin(theta1)*sin(theta2)))
end

"""Corrects for Ewald sphere curvature"""
function alpha_star(alpha::Float64, k1::Int64, k2::Int64, dq::Float64, lambda::Float64)
    theta1 = acos(k1*dq*lambda/(4*pi))
    theta2 = acos(k2*dq*lambda/(4*pi))
    return acos(sin(theta1)*sin(theta2)*cos(alpha) + cos(theta1)*cos(theta2))
end

"""Interface that allows to treat images as dissections through the 3D Fourier space"""
type ImageVolume <: Volume
    img::Array{Float64,2}
    size::Tuple{Int64,Int64}
    rmax::Float64

    function ImageVolume(img::Array{Float64,2}, rmax::Float64)
        new(img, Base.size(img), rmax)
    end

end

"""Accesing ImageVolume per position"""
function getVolumeInterpolated(img::ImageVolume, pos::Vector{Float64})
    center = Float64[img.size[1]/2.0, img.size[2]/2.0] .+ 0.5
    half_width = center .- 1.0
    shifted_pixel_pos = clamp.(round.(Int64, (pos[1:2]./img.rmax).*half_width + center), 1, minimum(img.size))
    return img.img[shifted_pixel_pos[1], shifted_pixel_pos[2]]
end

function maximum(img::ImageVolume)
    return Base.maximum(img.img)
end

"""
Returns a scattering image for one random orientation given a cubic Intensity volume
@param volume:
@params qcut:
@params number_incident_photons
@params rot:
@params gamma: contribution of noise as ratio to intensity
@params sigma: width of noise
"""
function pointsPerOrientation(volume::Volume, qcut::Float64, envelope_sigma::Float64, number_incident_photons::Int64; incident_photon_variance::Int64 = 0, rot::Matrix{Float64}=random_rotation(3), lambda::Float64=1.0, beamstop_width::Float64=0.0, print_warning::Bool=true)
    maxIntensity = real(maximum(volume))
    if typeof(volume) == SurfaceVolume
        maxIntensity *= 1.3
    end

    photons_per_image = number_incident_photons
    #Introduce a variable number of photons per image
    if incident_photon_variance > 0
        photon_dist = Distributions.Normal(number_incident_photons, incident_photon_variance)
        photons_per_image = round(Int64,Distributions.rand(photon_dist, 1)[1])
    end

    #Generate gaussian distributed incident photons
    dist = Distributions.Normal(0.0, envelope_sigma)
    incidents = Distributions.rand(dist, 3, photons_per_image)

    #Delete z dimension
    incidents[3,:] = 0.0

    #Make array of points
    incidents = Vector{Float64}[ incidents[:,i] for i = 1:photons_per_image]

    #Filter out photons that are in the beamstop
    filter!((p)-> abs(p[1]) > beamstop_width && abs(p[2]) > beamstop_width, incidents)
    #e.g. beamstop_width=qcut/20.0

    #Define envolpe function g(x)
    dist_3D = Distributions.MvNormal(zeros(3), [envelope_sigma,envelope_sigma,envelope_sigma])
    gx = (x) -> pdf(dist_3D, x)
    max_gauss = gx([0.0, 0.0, 0.0])

    #counts how many photons were above the envelope Gaussian
    overflow_count = 0

    accepted = Vector{Float64}[]
    for i = 1:length(incidents)
        p = incidents[i]

        if norm(p) <= qcut
            rp = rot*detector_to_Ewald_sphere(p, lambda)
            intensity_value = real(getVolumeInterpolated(volume, rp))
            ratio = intensity_value/(maxIntensity/max_gauss*gx(rp))
            if ratio > 1.0
                overflow_count += 1
                if print_warning
                    println(STDERR, "Probability ratio is larger than 1.0: ratio:$ratio, norm:$(norm(rp)), incident_photons:$(number_incident_photons), sigma:$(envelope_sigma), qcut:$(qcut).")
                end
            end
            if rand() < ratio
                push!(accepted, p)
            end
        end
    end
    println("Overflow percentage: $(overflow_count/number_incident_photons*100.0)%, scattered photons: $(length(accepted))")
    return accepted,rot
end

"""Uniform random points except beamstop"""
function beamstop_shot(qcut::Float64, photons_per_image::Int64, beamstop_width::Float64=0.0)

    incidents = 2*qcut*(rand(3, photons_per_image) - 0.5)

    #Delete z dimension
    incidents[3,:] = 0.0

    #Make array of points
    incidents = Vector{Float64}[ incidents[:,i] for i = 1:photons_per_image]

    #Filter out photons that are in the beamstop
    filter!((p)-> abs(p[1]) > beamstop_width && abs(p[2]) > beamstop_width, incidents)
    filter!((p)-> norm(p) <= qcut, incidents)

    return incidents
end


"""Calculates, how often a triple of k1,k2,k3 appears in the three photon correlation"""
function tripletFactor(k1::Int64,k2::Int64,k3::Int64)
    if k1 == k2 == k3
        return 1.0
    elseif (k1 == k2 != k3) || (k1 == k3 != k2) || (k2 == k3 != k1)
        return 3.0
    elseif (k1 != k2) && (k1 != k3) &&  (k2 != k3)
        return 6.0
    end
end

"""Calculates, how often a pair of k1,k2 appears in the two photon correlation"""
function doubletFactor(k1::Int64,k2::Int64)
    if k1 != k2
        return 2.0
    else
        return 1.0
    end
end

# """This calculates the 3 photon correlation for a sparse image"""
# function histogramCorrelationsInPicture_random(picture::Array{Array{Float64,1},1}, c2::C2, c3::C3, dq::Float64, N::Int64)
#   da = pi/N
#   l = length(picture)
#   if l < 3 return end
#   indices = collect(1:l)
#
#   for i = 1:Int64(l*(l+1)/2)
#     #Choose random triplet and sort according to k1>k2
#     shuffle!(indices)
#     doublet = Array{Float64}[picture[indices[j]] for j = 1:2 ]
#     sort!(doublet,lt=function(a,b) norm(a)>norm(b) end)
#     p1, p2 = doublet[1], doublet[2]
#
#     #Calculate k1, k2, alpha
#     k1, k2 = round(Int64,norm(p1)/dq), round(Int64,norm(p2)/dq)
#     alpha = angle_between_simple(p1,p2)
#     ai = Int64(mod(floor(Int64, alpha/da),N)+1)
#
#     @inbounds c2[ai,k2,k1] += 1.0 / (doubletFactor(k1,k2)*k1*k2)
#     @inbounds c2[N-ai+1,k2,k1] += 1.0 / (doubletFactor(k1,k2)*k1*k2)
#   end
#
#   for i = 1:Int64(l*(l+1)*(l+2)/6)
#     #Sort them according to k1 > k2 > k3
#     shuffle!(indices)
#     triplet = Array{Float64,1}[picture[indices[j]] for j = 1:3 ]
#     sort!(triplet,lt=function(a,b) norm(a)>norm(b) end)
#     p1,p2,p3 = triplet[1],triplet[2],triplet[3]
#
#     #Calculate k1, k2, k3, alpha, beta
#     k1,k2,k3 = round(Int64,norm(p1)/dq),round(Int64,norm(p2)/dq),round(Int64,norm(p3)/dq)
#     alpha,beta = mod(angle_between(p1,p2), pi),mod(angle_between(p1,p3), pi)
#     ai,bi = Int64(mod(floor(Int64, alpha/da),N)+1), Int64(mod(floor(Int64, beta/da),N)+1)
#
#     #histogram
#     @inbounds c3[ai,bi,k3,k2,k1] += 1.0 / (tripletFactor(k1,k2,k3)*k1*k2*k3)
#     @inbounds c3[N-ai+1,N-bi+1,k3,k2,k1] += 1.0 / (tripletFactor(k1,k2,k3)*k1*k2*k3)
#   end
# end


"""This calculates the 3 photon correlation for a sparse image"""
function histogramCorrelationsInPicture_alltoall(picture::Vector{Vector{Float64}}, c1::C1, c2::C2, c3::C3, dq::Float64, N::Int64, K2::Int64, K3::Int64, lambda::Float64)
    da = pi/N

    #Define picture filter function
    picturefilter = function(p::Vector{Float64})
        k = round(Int64,norm(p)/dq)
        return k >= 1 && k <= K2
    end

    filter!(picturefilter,picture)
    l = length(picture)
    if l < 3 return end
    picture = sort(picture,lt=function(a,b) norm(a)<norm(b) end)

    for i = 1:l
        p1 = picture[i]
        k1 = round(Int64,norm(p1)/dq)

        c1[k1] += Float64(1.0 / k1)

        for j = 1:(i-1)

            p2 = picture[j]
            k2 = round(Int64,norm(p2)/dq)

            alpha = angle_between(p1,p2)
            alpha_s = alpha
            if alpha > pi
                alpha_s = 2*pi - alpha
            end
            ai = Int64(mod(floor(Int64, alpha_s/da),N)+1)

            @fastmath val2 = Float64(1.0 / (doubletFactor(k1,k2)*k1*k2))
            @inbounds c2[ai,k2,k1] += val2

            for k = 1:(j-1) #this implies k2 >= k3

                p3 = picture[k]
                k3 = round(Int64,norm(p3)/dq)

                if k1 <= K3 && k2 <= K3 && k3 <= K3

                    beta = angle_between(p1,p3)
                    if alpha > pi
                        beta = 2*pi - beta
                    end
                    bi = Int64(mod(floor(Int64, beta/da),2*N)+1)

                    @fastmath val3 = Float64(1.0 / (tripletFactor(k1,k2,k3)*k1*k2*k3))
                    @inbounds c3[ai,bi,k3,k2,k1] += val3
                end
            end
        end
    end
end

"""Noise volume in CubeVolume"""
function get_noise_volume(intensity::CubeVolume, sigma::Float64)
    volume = deepcopy(intensity)
    #Let's create a normalized noise volume
    r = linspace(-volume.rmax, volume.rmax, volume.cubesize)
    noise_volume = Float64[ gaussian_distribution([x,y,z], [0.0, 0.0, 0.0], [sigma,sigma,sigma]) for x=r, y=r, z=r]
    noise_volume /= sum(abs, noise_volume)
    volume.cube = noise_volume
    return volume
end

"""Noise volume in SphericalHarmonicsVolume"""
function get_noise_volume(intensity::SurfaceVolume, sigma::Float64)
    intensity_sh = getSphericalHarmonicsVolume(intensity)
    intensity_sh = deleteTerms(intensity_sh, 0, 0)
    for k = 1:intensity.KMAX
        q = k*dr(intensity)
        setc(intensity_sh,k,0,0, gaussian_distribution(q, 0.0, sigma))
    end
    return getSurfaceVolume(intensity_sh)
end

#http://rcwww.kek.jp/research/shield/photon_r.pdf
#http://ac.els-cdn.com/S0969806X15300621/1-s2.0-S0969806X15300621-main.pdf?_tid=1dcbacf4-3a16-11e7-87fd-00000aab0f6c&acdnat=1494925412_b9ac6c901a449780b30c1076c1a57ccc
# """Calculates the compton scattering distribution"""
function compton_scattering(theta::Float64)
    r0_in_m = 2.817940e-15 #classical electron radius [m]
    r0 = r0_in_m * 1.0e10 #classical electron radius [A]
    mec2 = 0.5109989461#electron rest mass energy [MeV] : m_e * c^2
    hnu0 = 0.005#energy of incident photon [MeV] : 5keV
    k0 = hnu0 / mec2
    hnu1 = hnu0 / (1 + k0 * (1 - cos(theta)))#energy of outbound photon
    k1 = hnu1 / mec2
    #Klein and Nishina differential cross section formula
    #From "Compton scattering I: Angular distribution and polarization degree" paper
    #     return 0.25 * r0^2 * (k1/k0)^2 * (k1/k0 + k0/k1 - 2 + 4*sin(theta)^2)
    #From wikipedia
    return 0.5 * r0^2 * (k1/k0)^2 * (k1/k0 + k0/k1 - sin(theta)^2)
end

function compton_scattering_q(q::Float64)
    h = 4.135667516e-15# [eV*s]
    c = 299792458 #[m/s]
    E = 5000 #[eV]
    lambda_in_m = h*c / E #wavelength [m]
    lambda = lambda_in_m * 1.0e10 #wavelength [A]
    val = q * lambda / (4*pi)
    if abs(val) > 1.0
        return 0.0
    else
        theta = asin(val)
        return compton_scattering(theta) * (4*pi/lambda * cos(theta))
    end
end

function get_compton_noise_volume(intensity::CubeVolume)
    volume = deepcopy(intensity)
    #Let's create a normalized noise volume
    r = linspace(-volume.rmax, volume.rmax, volume.cubesize)
    noise_volume = Float64[ compton_scattering_q(norm([x,y,z])) for x=r, y=r, z=r]
    noise_volume /= sum(abs, noise_volume)
    volume.cube = noise_volume
    return volume
end

"""
Generates pictures and histograms the triplets in these pictures without reusing photons
@params intensity:
@params qcut:
@params K:
@params N:
@params maxpictures:
@params number_incident_photons:
@params numprocesses:
@params three:
@params two:
@params gamma:
@params sigma:
"""
function generateHistogram(intensity::Volume; qcut::Float64=1.0, K2::Int64=38, K3::Int64=26, N::Int64=32, max_pictures::Int64=Int64(0), number_incident_photons::Int64=1000, incident_photon_variance::Int64 = 0, numprocesses::Int64=1, file::String="histo.dat", noise::Noise=GaussianNoise(0.0, 1.0, false), batchsize::Int64 = 1000, histogramMethod=histogramCorrelationsInPicture_alltoall, number_particles::Int64=1, lambda::Float64=1.0, beamstop_width::Float64=0.0, beamstop_only::Bool=false)

    #cutoff parameters
    dq = qcut/K2

    c1_full = zeros(Float64, K2)
    c2_full = zeros(Float64,(N, K2,K2))
    c3_full = zeros(Float64,(N,2*N,K3,K3,K3))
    params = Dict("num_pictures"=>0, "num_incident_photons"=>number_incident_photons, "noise"=>noise, "qcut"=>qcut, "K2"=>K2, "K3"=>K3, "N"=>N, "qcut"=>qcut, "dq"=>dq)

    #Load preexisting histogram file so we can continue from where we left off
    if isfile(file)
        println("Continue from preexisting correlation in $file.")
        (params, c2_full, c3_full, c1_full) = deserializeFromFile(file)
    end

    #Calculate noise volume if necessary
    noise_volume = get_noise_volume(intensity, noise.sigma)

    #Count current triplets
    current_triplets = countTriplets(c3_full)

    println("Starting to collect images")
    flush(STDOUT)

    #Make batches of pictures so we can save often
    while max_pictures > 0 && params["num_pictures"] < max_pictures

        c1_part,c2_part,c3_part = @sync @parallel ( (a,b) -> (a[1]+b[1], a[2]+b[2], a[3]+b[3])) for i = 1:numprocesses
            c1 = zeros(Float64,K2)
            c2 = zeros(Float64,N,K2,K2)
            c3 = zeros(Float64,N,2*N,K3,K3,K3)

            for j=1:batchsize
                photon_list = Vector{Float64}[]

                if beamstop_only == false
                    for n=1:number_particles
                        single_molecule,rot = pointsPerOrientation(intensity, qcut, qcut/3.0, number_incident_photons, incident_photon_variance=incident_photon_variance, lambda=lambda, beamstop_width=beamstop_width)

                        if noise.gamma > 0.0
                            noise_photons,_ = pointsPerOrientation(noise_volume,qcut, noise.sigma*1.05, noise.photons, incident_photon_variance=0, rot=rot, lambda=lambda, beamstop_width=beamstop_width)
                            append!(single_molecule, noise_photons)
                        end
                        append!(photon_list, single_molecule)
                    end
                else
                    #beamstop only
                    append!(photon_list, beamstop_shot(qcut, number_incident_photons, beamstop_width))
                end

                histogramMethod(photon_list, c1, c2, c3, dq, N, K2, K3, lambda)
            end
            (c1, c2, c3)
        end
        c1_full += c1_part
        c2_full += c2_part
        c3_full += c3_part

        params["num_pictures"] += numprocesses*batchsize #In numprocesses, we created 'batchsize' pictures
        current_triplets = countTriplets(c3_full)
        serializeToFile(file, (params, c2_full, c3_full, c1_full))

        println("triplets: $current_triplets\tpictures: $(params["num_pictures"]) / $(max_pictures) ($(params["num_pictures"]/max_pictures*100.0) perc.)")
        flush(STDOUT)
    end
end

"Loading a histogram from a file and store two/three photon histograms in global variables."
function loadHistograms(K2::Int64, K3::Int64, file::String, load_c1::Bool=true)

    params, c2_full, c3_full, c1_full = 0,0,0,0
    #Try loading two different structure types
    if load_c1
        params, c2_full, c3_full, c1_full = deserializeFromFile(file)
    else
        params, c2_full, c3_full = deserializeFromFile(file)
    end

    println("Loaded $(countDoublets(c2_full)) doublets and $(countTriplets(c3_full)) triplets from $file generated from $(params["num_pictures"]) pictures.")

    #Doubling old histograms
    NA,NB,K3,_,_ = Base.size(c3_full)
    if NA == NB
        println("Doubling old histograms.")
        new_c3 = zeros(NA, 2*NB,K3,K3,K3)
        for k1=1:K3
            for k2=1:k1
                for k3=1:k2
                    new_c3[:,1:NB,k3,k2,k1] = c3_full[:,:,k3,k2,k1]
                    new_c3[:,NB+1:2*NB,k3,k2,k1] = c3_full[:,:,k3,k2,k1]
                end
            end
        end
        c3_full = new_c3
    end

    c1_full = max.(c1_full, 1e-30)
    c2_full = max.(c2_full, 1e-30)
    c3_full = max.(c3_full, 1e-30)

    #Cut down to designated K
    c2 = c2_full[:,1:K2,1:K2]
    c2 = c2 / sum(abs, c2)
    c3 = c3_full[:,:,1:K3,1:K3,1:K3]
    c3 = c3 / sum(abs, c3)

    if load_c1
        return c2_full, c2, c3_full, c3, c1_full
    else
        return c2_full, c2, c3_full, c3
    end
end

"""Calculates the total number of triplets in a histogram"""
function countTriplets(c3::C3)
    s = size(c3)
    replace_nan!(c3)
    replace_inf!(c3)
    return sum(abs, c3[a,b,k3,k2,k1]*Float64(tripletFactor(k1,k2,k3)*k1*k2*k3) for a=1:s[1],b=1:s[2],k3=1:s[3],k2=1:s[4],k1=1:s[5])
end

"""Calculates the total number of doublets in a histogram"""
function countDoublets(c2::C2)
    s = size(c2)
    replace_nan!(c2)
    replace_inf!(c2)
    return sum(abs, c2[a,k2,k1]*doubletFactor(k1,k2)*k1*k2 for a=1:s[1],k2=1:s[2],k1=1:s[3])
end
