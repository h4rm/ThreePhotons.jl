export 	images_10p,
	calculate_expected_triplets,
	calculate_required_images,
	calculate_images_ppi,
	calculate_incident_photons,
	minutues_to_measure,
	cor_ISC,
	cor_FSC,
	angle_between,
	angle_between_simple,
	euler,
	get_euler_angles,
	sphericalToCartesian,
	cartesianToSpherical,
	gaussian,
	gaussian_distribution,
	normal_distribution,
	replace_NaN!,
	random_rotation,
	random_rotation_step,
	SON_parametrized,
	SON_parametrized_derivative,
	random_doublet,
	random_triplet,
	serializeToFile,
	deserializeFromFile,
    gaussian_filter,
	generate_gaussian_kernel,
	calculate_ppi

images_10p = Int64[1280000,5120000,20480000,81920000,327680000,3276800000]
calculate_expected_triplets = (images,ppi::Int64) -> images * (ppi * (ppi-1) * (ppi-2) / 6.0)
calculate_required_images = (triplets,ppi::Int64) -> triplets / (ppi * (ppi-1) * (ppi-2) /6.0)
calculate_images_ppi = (ppi::Int64) -> map((img)->ceil(Int64,calculate_required_images(calculate_expected_triplets(img, 10), ppi)/8)*8, images_10p)
calculate_incident_photons = (ppi::Int64, incident_10p::Int64=550) -> round(Int64, (ppi / 10.0) * incident_10p)
calculate_ppi = (incident_photons::Int64, incident_10p::Int64=550) -> round(Int64, (incident_photons * 10.0) / incident_10p)

minutues_to_measure = (images, pulses_per_second::Int64=27000, hitrate::Float64=1.0) -> images / (pulses_per_second*60) / hitrate

###############################################################
#                     MATH
###############################################################

"""Pearson correlation between two vectors"""
function cor_ISC(a::Vector{Complex{Float64}}, b::Vector{Complex{Float64}})
    return cor(real(a), real(b))
end

"""Fourier shell correlation between two complex vectors"""
function cor_FSC(a::Vector{Complex{Float64}}, b::Vector{Complex{Float64}})
    na = norm(a)
    nb = norm(b)
    if na > eps() && nb > eps()
        return abs(sum(a.*conj(b))/(na*nb))
    else
        return 0.0
    end
end

"Angle between two vectors ranging from 0 to 2.0*pi"
function angle_between(p1::Vector{Float64}, p2::Vector{Float64})
    val = clamp(Base.dot(p1,p2)/ (norm(p1) * norm(p2)), -1.0, 1.0)
    if abs(val) > 0.0
        angle = acos(val)
        if cross(p1, p2)[3] > 0
            angle =  2.0*pi - angle
        end
        return angle
    else
        return 0.0
    end
end

"Angle between two vectors ranging from 0 to pi"
function angle_between_simple(p1::Vector{Float64}, p2::Vector{Float64})
    angle = acos( clamp(Base.dot(p1,p2)/ (norm(p1) * norm(p2)), -1.0, 1.0))
    return angle
end

"""
Calculate a three dimensional rotation matrix from the euler angles.

@param a: alpha, angle between the x-axis and the line of nodes
@param b: beta, angle between the z axis of the different coordinate systems
@param c: gamma, angle between the line of nodes and the X-axis
"""
function euler(a::Float64, b::Float64, c::Float64)

    ca, cb, cc = cos(a), cos(b), cos(c)
    sa, sb, sc = sin(a), sin(b), sin(c)

    @fastmath return Float64[[ cc * cb * ca - sc * sa, cc * cb * sa + sc * ca, -cc * sb] [-sc * cb * ca - cc * sa, -sc * cb * sa + cc * ca, sc * sb] [ sb * ca, sb * sa, cb ]]
end

"Calculates the Euler angles of a 3D rotation matrix"
function get_euler_angles(m)
    a = mod(atan2(m[3, 2], m[3, 1]), 2.0*pi)
    b = mod(atan2((m[3, 1] + m[3, 2]) / (cos(a) + sin(a)), m[3, 3]), 2.0*pi)
    c = mod(atan2(m[2, 3], -m[1, 3]), 2.0*pi)

    return a, b, c
end

"""Converts spherical coordinates to cartesian coordinates"""
function sphericalToCartesian(r::Float64, phi::Float64, theta::Float64)
    x = r*cos(phi)*sin(theta)
    y = r*sin(phi)*sin(theta)
    z = r*cos(theta)
    return Float64[x,y,z]
end

"""Converts a vector to spherical coordinates
Returns r, phi,theta
"""
function cartesianToSpherical(vec::Vector{Float64})
    r::Float64 = norm(vec)
    phi::Float64 = atan2(vec[2],vec[1])
    if phi < 0.0 phi = phi + 2.0*pi end
    theta::Float64 = r > 0.0 ? acos(clamp(vec[3]/r, -1.0, 1.0)) : 0.0
    return r,phi,theta
end

"Creates a gaussian distributed number with sigma and mu"
function gaussian(mu::Float64, sigma::Float64)
    return randn()*sigma + mu
end

"""Calcualtes the gaussian distribution wiht sigma and mu"""
function gaussian_distribution(x::Float64, mu::Float64, sigma::Float64)
    @assert sigma != 0.0
    s2 = 2.0 * sigma^2
    1/sqrt(pi*s2)*exp(-norm(x-mu)^2/s2)
end

"""Calcualtes the gaussian distribution wiht sigma and mu"""
function gaussian_distribution(x::Vector{Float64}, mu::Vector{Float64}, sigma::Vector{Float64})
    @assert sigma != 0.0
    return pdf(Distributions.MvNormal(mu, sigma), x)
end

"""Calcualtes the normal distribution wiht sigma=1.0 and mu=0.0"""
function normal_distribution(x::Vector)
    gaussian_distribution(x, 0.0, 1.0)
end

"""In a numerical Vector, replaces NaN with zero"""
function replace_NaN!(vec::Vector)
    for i = 1:length(vec)
        vec[i] = isnan(vec[i]) ? 0.0 : vec[i]
    end
end

"""Calculates a random rotation in `dim`
Source: http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf"""
function random_rotation(dim::Int64)
    Q,R = qr(randn(dim, dim))
    r = diag(R)
    r = r ./ abs.(r)
    L = diagm(r)
    res = Q*L
    if(det(res) > 0)
        return res
    else
        random_rotation(dim)
    end
end

"""returns a random rotation with a slight angle beta=pi/8
Taken from "Random rotations: characters and random walk on SO(N)"""
function random_rotation_step(b, beta=pi/8)
    R = eye(b)
    R[1,1] = cos(beta)
    R[1,2] = -sin(beta)
    R[2,1] = sin(beta)
    R[2,2] = cos(beta)

    C = random_rotation(b)
    B = inv(C)*R*C
    return B
end

"Calculates a rotation matrix in `dim` dimensions with dim*(dim-1)/2 rotation angles"
function SON_parametrized(dim, angles)
    R = eye(dim)
    assert(dim*(dim-1)/2 == length(angles))
    num = 1
    for i = 1:dim
        for j = i+1:dim
            M = eye(dim)
            c = cos(angles[num])
            s = sin(angles[num])
            M[i,i] = c
            M[j,j] = c
            M[i,j] = s
            M[j,i] = -s
            R *= M
            num += 1
        end
    end
    return R
end

"Calculates the derivative of R by the kth parameter"
function SON_parametrized_derivative(dim, angles, k)
    R = eye(dim)
    num = 1
    for i = 1:dim
        for j = i+1:dim
            c = cos(angles[num])
            s = sin(angles[num])
            if(num != k)
                M = eye(dim)
                M[i,i] = c
                M[j,j] = c
                M[i,j] = s
                M[j,i] = -s
                R *= M
            else
                M = zeros(dim,dim)
                M[i,i] = -s
                M[j,j] = -s
                M[i,j] = c
                M[j,i] = -c
                R *= M
            end
            num += 1
        end
    end
    return R
end

"""Creates random doublet"""
function random_doublet(K::Int64=8, half::Bool=true)
    k1 = rand(1:K)
    k2 = rand(1:K)
    if half
        if k1 >= k2
            return (k1, k2)
        else
            return random_doublet(K, half)
        end
    else
        return (k1,k2)
    end
end

"""Creates random triplet"""
function random_triplet(K=8, half::Bool=true)
    k1 = rand(1:K)
    k2 = rand(1:K)
    k3 = rand(1:K)
    if half
        if k1 >= k2 && k2 >= k3
            return (k1,k2,k3)
        else return random_triplet(K, half)
        end
    else
        return (k1,k2,k3)
    end
end

###############################################################
#                     FILES
###############################################################
using FileIO

"""Serializes the data to a file"""
function serializeToFile(filename::String, data)
    # out = open(filename, "w+")
    # serialize(out, data)
    # close(out)
	save("$(filename).jld", "data", data)
end

"""Deserializes a julia data object from a file using JLD.jl and fixing backward compatibility"""
function deserializeFromFile(filename::String)
	# try
	# 	return deserializeFromFile_raw(filename)
	# catch
		try
	    	data = load(filename)
			return data["data"]
		catch err
			if isa(err, FileIO.UnknownFormat{FileIO.File{FileIO.DataFormat{:UNKNOWN}}}) || isa(err, FileIO.File{FileIO.DataFormat{:UNKNOWN}})
				println("Trying to rewrite $filename into JLD format and reload")
				run(`julia5 -e "using ThreePhotons5; using JLD; data = deserializeFromFile(\"$(filename)\"); save(\"$(filename).jld\", \"data\", data)"`)
				try
					data = load("$(filename).jld")
					return data["data"]
				catch
					println("Something failed loading $filename")
				end
			else
				error("Failed loading with error: $(err)")
			end
		end
	# end
end

function deserializeFromFile_raw(filename::String)
    out = open(filename, "r")
    data = deserialize(out)
    close(out)
    return data
end

###############################################################
#                     MISC
###############################################################

"""Say something via MacOS say function"""
function say(str::String)
    voices = ["Agnes",
    "Albert",
    "Alex",
    "Alice",
    "Alva",
    "Amelie",
    "Anna",
    "Bad",
    "Bahh",
    "Bells",
    "Boing",
    "Bruce",
    "Bubbles",
    "Carmit",
    "Cellos",
    "Damayanti",
    "Daniel",
    "Deranged",
    "Diego",
    "Ellen",
    "Fiona",
    "Fred",
    "Good",
    "Hysterical",
    "Ioana",
    "Joana",
    "Junior",
    "Kanya",
    "Karen",
    "Kathy",
    "Kyoko",
    "Laura",
    "Lekha",
    "Luciana",
    "Maged",
    "Mariska",
    "Mei-Jia",
    "Melina",
    "Milena",
    "Moira",
    "Monica",
    "Nora",
    "Paulina",
    "Pipe",
    "Princess",
    "Ralph",
    "Samantha",
    "Sara",
    "Satu",
    "Sin-ji",
    "Tessa",
    "Thomas",
    "Ting-Ting",
    "Trinoids",
    "Veena",
    "Vicki",
    "Victoria",
    "Whisper",
    "Xander",
    "Yelda",
    "Yuna",
    "Zarvox",
    "Zosia",
    "Zuzana"]
    try
        voice = voices[rand(1:length(voices))]
        # println(voice)
        run(`say $str --voice $voice`)
    end
end

function whisper(str)
    try
        run(`say $str --voice 'whisper'`)
    end
end

ϕ(x::Float64, σ::Float64) = exp(-x^2 / (2 * σ^2)) / (σ * sqrt(2 * π))
ϕ(x::Vector{Float64}, σ::Float64) = exp(-norm(x)^2 / (2 * σ^2)) / (σ * sqrt(2 * π))

function generate_gaussian_kernel(σ::Float64, kernel_size::Vector{Int64})
    center = map((dim) -> mod(dim,2) == 0 ? (dim/2)+1 : ceil(dim/2), kernel_size)
	if length(kernel_size) == 1
    	return Float64[ϕ(x-center[1],σ) for x=1:kernel_size[1]]
	elseif length(kernel_size) == 2
		return Float64[ϕ(Float64[x,y]-center,σ) for x=1:kernel_size[1], y=1:kernel_size[2]]
	end
end

function generate_gaussian_kernel(σ::Float64, kernel_size::Tuple)
	generate_gaussian_kernel(σ, Int64[t for t in kernel_size])
end

function generate_gaussian_kernel(σ::Float64, kernel_size::Int64)
	return generate_gaussian_kernel(σ, [kernel_size])
end

function gaussian_filter(data::Vector{Float64}, σ::Float64=1.0)
    kernel = generate_gaussian_kernel(σ, Base.size(data))
    fkern = fft(ifftshift(kernel))
    fdata = fft((data))
    return real((ifft(fkern.*fdata)))
end

function gaussian_filter(data::Matrix{Float64}, σ::Float64=1.0)
    kernel = generate_gaussian_kernel(σ, Base.size(data))
    fkern = fft(ifftshift(kernel))
    fdata = fft((data))
    return real((ifft(fkern.*fdata)))
end
