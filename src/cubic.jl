#Abstract volume description
abstract type Volume end

"""Describing a volume with cubic voxels"""
type CubeVolume <: Volume
    cube::Array{Complex{Float64},3}
    cubesize::Int64
    rmax::Float64

    function CubeVolume(cube::Array{Complex{Float64},3}, cubesize::Int64, rmax::Float64)
        new(cube, cubesize, rmax)
    end
    function CubeVolume(cube::Array{Float64,3}, cubesize::Int64, rmax::Float64)
        new(complex(cube), cubesize, rmax)
    end
    function CubeVolume(cubesize::Int64, rmax::Float64)
        new(zeros(Complex{Float64}, cubesize, cubesize, cubesize), cubesize, rmax)
    end
end

"""Get voxel size of Cube Volume"""
function dr(cube::CubeVolume) return 2*cube.rmax/(cube.cubesize-1) end
function qmax(cube::CubeVolume) return pi/dr(cube) end

function +(a::CubeVolume, b::CubeVolume) return CubeVolume(a.cube + b.cube, a.cubesize, a.rmax) end
function -(a::CubeVolume, b::CubeVolume) return CubeVolume(a.cube - b.cube, a.cubesize, a.rmax) end
function *(a::CubeVolume, b::Number) return CubeVolume(a.cube * b, a.cubesize, a.rmax) end
function *(a::Number, b::CubeVolume) return CubeVolume(b.cube * a, b.cubesize, b.rmax) end
function /(a::CubeVolume, b::Number) return CubeVolume(a.cube / b, a.cubesize, a.rmax) end
function sum(abs,a::CubeVolume) return sum(abs, a.cube) end
function abs(a::CubeVolume) return CubeVolume(abs(a.cube), a.cubesize, a.rmax) end
function real(a::CubeVolume) return Cubevolume(real(a.cube), a.cubesize, a.rmax) end
function maximum(a::CubeVolume) return Base.maximum(real(a.cube)) end

"""Calculates the density cube from atomlist coordinates"""
function getDensityCube(atomlist, cubesize::Int64, rmax::Float64)
    r = linspace(-rmax, rmax, cubesize)
    cube = Complex{Float64}[getDensity_cartesian(atomlist, [x,y,z]) for x=r, y=r, z=r]
    return CubeVolume(cube,cubesize, rmax)
end

"""Average a list of cubes"""
function average(list::Array{CubeVolume})
    return reduce(+, list) / length(list)
end

"""Center of mass calculation for fitting in real space."""
function calculate_center_of_mass(volume::CubeVolume)
    r = linspace(-volume.rmax, volume.rmax, volume.cubesize)
    positions = Vector{Float64}[Float64[x,y,z] for x=r, y=r, z=r ]
    cube = real(volume.cube)
    convert(Vector{Float64},sum(positions .* cube)/sum(cube))
end

"""Translates the content of a cube into a certain direction"""
function translate_cube(volume::CubeVolume, direction::Vector{Float64})
    r = linspace(-volume.rmax, volume.rmax, volume.cubesize)
    return CubeVolume(Complex{Float64}[getVolumeInterpolated(volume, Float64[x,y,z]+direction) for x=r, y=r, z=r], volume.cubesize, volume.rmax)
end

"""Fits a cube into another cube"""
function center_cube(volume::CubeVolume)
    translated_vol = translate_cube(volume, calculate_center_of_mass(volume) )
    return translated_vol
end

"""Translates the volume to the center of mass of the reference"""
function displacement_fitting(volume::CubeVolume, reference::CubeVolume)
    translated_vol = translate_cube(volume, -(calculate_center_of_mass(reference) - calculate_center_of_mass(volume)) )
    return translated_vol,similarity(translated_vol, reference)
end

# """Calculates the Intensity cube from direct expression"""
# function getIntensityCube(cubesize)
#     r = linspace(-qmax, qmax, cubesize)
#     cube = Complex{Float64}[getIntensity_cartesian([x,y,z]) for x=r, y=r, z=r]
#     return cube
# end

"""cubic forward fourier transform"""
function forward(volume::CubeVolume)
    ftcube = fftshift(fft(ifftshift(volume.cube)))
    return CubeVolume(ftcube, volume.cubesize, qmax(volume))
end

"""Calculates the aboslute square of volume"""
function absoluteSquare(volume::CubeVolume)
    return CubeVolume(abs.(volume.cube).^2, volume.cubesize, volume.rmax)
end

"""Mirrors the cube along one axis"""
function mirrorCube(volume::CubeVolume)
    cubesize = volume.cubesize
    mx,my,mz = deepcopy(volume),deepcopy(volume),deepcopy(volume)

    # mx.cube = Complex{Float64}[volume.cube[cubesize-x+1,y,z] for x=1:cubesize,y=1:cubesize,z=1:cubesize]
    # my.cube = Complex{Float64}[volume.cube[x,cubesize-y+1,z] for x=1:cubesize,y=1:cubesize,z=1:cubesize]
    mz.cube = Complex{Float64}[volume.cube[x,y,cubesize-z+1] for x=1:cubesize,y=1:cubesize,z=1:cubesize]

    # return [mx,my,mz]
    return mz
end

"Calculates the similarity between two cubes"
function similarity(volume1::CubeVolume, volume2::CubeVolume)
    vec1 = real(reshape(volume1.cube, volume1.cubesize^3))
    vec2 = real(reshape(volume2.cube, volume2.cubesize^3))
    vec1 = vec1 / norm(vec1)
    vec2 = vec2 / norm(vec2)
    return cor(vec1, vec2)
end

"""Saves a cubic structure to a .mrc file"""
function saveCube(volume::CubeVolume, filename::String)
    cubesize = volume.cubesize
    rmax = volume.rmax
    cube = convert(Array{Float32}, real(volume.cube))
    cube = cube / sum(abs,cube)
    offset = floor(cubesize/2)

    println("Saving cube to $filename with cubesize=$cubesize")
    out = open(filename, "w+")
    write(out, Int32[cubesize, cubesize, cubesize]) #columns, rows, sections
    write(out, Int32[2]) #mode
    write(out, Int32[-offset,-offset,-offset]) #starting point
    write(out, Int32[cubesize,cubesize,cubesize]) #grid size
    write(out, Float32[2.0*rmax, 2.0*rmax, 2.0*rmax]) #cell size
    write(out, Float32[90.0, 90.0, 90.0]) #cell angles
    write(out, Int32[1,2,3]) #mapc mapr, maps
    write(out, Float32[minimum(cube), maximum(cube),mean(cube)])
    write(out, Int32[0, 0]) #ispg next
    write(out, zeros(Int8, 1024-96))

    write(out, cube)
    close(out)
end

"""Loads a cube from a .mrc file"""
function loadCube(filename, layers=0)
    out = open(filename, "r")
    cubesizes = convert(Array{Int64},read(out, Int32, 3)) #columns, rows, sections
    mode = read(out, Int32) #mode
    startingpoint = read(out, Int32, 3) #starting point
    gridsize = read(out, Int32, 3) #grid size
    cellsize = read(out, Float32, 3) #cell size
    cellangles = read(out, Float32, 3) #cell angles
    mapparameters = read(out, Int32, 3) #mapc mapr, maps
    maxmimean = read(out, Float32, 3)
    ispgnext = read(out, Int32, 2) #ispg next
    stuff = read(out, Int8, 1024-96)
    length = layers > 0 ? layers : cubesizes[3]
    data = read(out, Float32, cubesizes[1]*cubesizes[2]*length)
    close(out)
    return CubeVolume(reshape(convert(Array{Complex{Float64}}, data),cubesizes[1], cubesizes[2], length), cubesizes[1], convert(Float64,cellsize[1]/2.0) )
end

"""Gets the value within a cube via trilinear interpolation"""
function getVolumeInterpolated(volume::CubeVolume, vec::Vector{Float64})
    cubesize = volume.cubesize
    rmax = volume.rmax
    deltar = dr(volume)
    cube = volume.cube

    x = clamp((vec[1] + rmax)/deltar + 1.0, 1.0, Float64(cubesize))
    y = clamp((vec[2] + rmax)/deltar + 1.0, 1.0, Float64(cubesize))
    z = clamp((vec[3] + rmax)/deltar + 1.0, 1.0, Float64(cubesize))

    x0,x1 = floor(Int64, x), ceil(Int64, x)
    y0,y1 = floor(Int64, y), ceil(Int64, y)
    z0,z1 = floor(Int64, z), ceil(Int64, z)

    a = Float64(x-x0)
    b = Float64(y-y0)
    c = Float64(z-z0)

    V000 = (1.0-a)*(1.0-b)*(1.0-c)
    V100 = a *(1.0-b)*(1.0-c)
    V010 = (1.0-a)* b*(1.0-c)
    V001 = (1.0-a)*(1.0-b)*c
    V101 = a*(1.0-b)*c
    V011 = (1.0-a)*b*c
    V110 = a* b*(1.0-c)
    V111 = a* b* c

    @inbounds @fastmath return V000 * cube[x0,y0,z0] + V100 * cube[x1,y0,z0] + V010 * cube[x0,y1,z0] + V001 * cube[x0,y0,z1] + V101 * cube[x1,y0,z1] + V011 * cube[x0,y1,z1] + V110 * cube[x1,y1,z0] + V111 * cube[x1,y1,z1]
end

"""Rotates a cubic volume"""
function rotateStructure(volume::CubeVolume, theta::Float64, phi::Float64, gamma::Float64)
    newvol = deepcopy(volume)
    r = linspace(-volume.rmax, volume.rmax, volume.cubesize)
    newvol.cube = Complex{Float64}[getVolumeInterpolated(volume, euler(theta, phi, gamma)*Float64[x,y,z]) for x=r, y=r, z=r]
    return newvol
end

"""Calculates the shell-wise correlation of two cubes"""
function shell_correlation(volume1::CubeVolume, volume2::CubeVolume, K::Int64=0, cor_method=cor_FSC; N::Int64=48)
    K = (K == 0 ? floor(Int64,volume1.cubesize/2) : K)
    fsc = zeros(Float64,K)

    #let's sample the sphere evenly
    for k = 1:K
        r = k*dr(volume1)
        surfA = Complex{Float64}[ getVolumeInterpolated(volume1, sphericalToCartesian(r, phi, theta)) for phi=linspace(0,2.0*pi,N), theta=acos.(linspace(-1.0,1.0,N)) ]
        surfB = Complex{Float64}[ getVolumeInterpolated(volume2, sphericalToCartesian(r, phi, theta)) for phi=linspace(0,2.0*pi,N), theta=acos.(linspace(-1.0,1.0,N)) ]

        surfA = reshape(surfA,N^2)
        surfB = reshape(surfB,N^2)
        fsc[k] = cor_method(surfA, surfB)
    end

    return fsc
end

function shell_correlation_ISC(volume1::CubeVolume, volume2::CubeVolume, K::Int64=0)
    shell_correlation(volume1, volume2, K, cor_ISC)
end

function shell_correlation_FSC(volume1::CubeVolume, volume2::CubeVolume, K::Int64=0)
    shell_correlation(volume1, volume2, K, cor_FSC)
end

"""Calculate the Intensity Shell Correlation between two cubes"""
function ISC(density::CubeVolume, intensity::CubeVolume, K::Int64=0)
    return shell_correlation_ISC(intensity, absoluteSquare(forward(density)), K)
end

"""Calculate the Fourier Shell Correlation between two cubes"""
function FSC(density::CubeVolume, fourier::CubeVolume, K::Int64=0)
    return shell_correlation_FSC(fourier, forward(density), K)
end
