#Atom factors for specific atom types
type atomfactor
    R::Float64
    B::Float64
end

#Constant directory with atom information
atomfactorlist = Dict( "C" => atomfactor(0.7, 6.0),
                    "N" => atomfactor(0.65, 7.0),
                    "S" => atomfactor(1.0, 16.0),
                    "O" => atomfactor(0.6, 8.0),
                "default" => atomfactor(0.7, 6.0),
                    "P" => atomfactor(0.65, 15.0),
                    "H" => atomfactor(0.3, 1.0),
                    "D" => atomfactor(0.3, 1.0)
)

"""Atom type with position and type"""
type atomType
    pos::Array{Float64}
    t::String
    group::String
end

"""Loads the atom coordinates from a pdb file and returns centered list with atomic coordinates"""
function loadPDB(filename)
    file = open(filename, "r")
    list = atomType[]
    #print(readstring(file))
    for line in readlines(file)
        if line[1:4] == "ATOM"
            x = float(line[31:38])
            y = float(line[39:46])
            z = float(line[47:54])
            t = line[78:78]
            group = line[18:20]
            push!(list,atomType(Float64[x,y,z],t, group))
        end
    end
    close(file)
    return centerStructure(list)
end

"""Centers list of atoms to the center of mass"""
function centerStructure(atomlist)
    center = zeros(3)
    for atom in atomlist
        center = center + atom.pos
    end

    center = center / length(atomlist)

    # shift by center
    for atom in atomlist
        atom.pos = atom.pos - center
    end
    return atomlist
end

"""Calculates density for a given point in spherical coordinates"""
function getDensity_cartesian(atomlist, pos::Array{Float64,1})
    result = 0.0
    for atom in atomlist
        rsq = norm(pos-atom.pos)^2
        if rsq < 5 && atom.t != "1"
            atomtype = atomfactorlist[atom.t]
            B = atomtype.B
            R = atomtype.R
            result += B/((R*sqrt(2.0*pi))^3) *exp(-rsq/(2.0*R^2))
        end
    end
    return result
end

"""Calculates intensity for a given point in spherical coordinates"""
function getIntensity_cartesian(atomlist, pos::Array{Float64,1})

    k2 = norm(pos)^2
    sumc = 0
    sums = 0
    for atom in atomlist
        R = atomfactorlist[atom.t].R
        B = atomfactorlist[atom.t].B
        expk = exp( -0.5 * R*R*k2)
        sumc += B* cos( dot(pos, atom.pos)) * expk
        sums += B* sin( dot(pos, atom.pos)) * expk
    end
    return sumc^2 + sums^2
end

"""Calculates intensity for a given point in spherical coordinates"""
function getFourier_cartesian(atomlist, pos::Array{Float64,1})

    k2 = norm(pos)^2
    sumc = 0
    sums = 0
    for atom in atomlist
        R = atomfactorlist[atom.t].R
        B = atomfactorlist[atom.t].B
        expk = exp( -0.5 * R*R*k2)
        sumc += B* cos( dot(pos, atom.pos)) * expk
        sums += B* sin( dot(pos, atom.pos)) * expk
    end
    return sumc + 1i*sums
end

##################################################################
# Structure handling
##################################################################

"""Loads a pdb file and creates a spherical harmonics description of the electron denisty, the Fourier transform and the Intensity (absolute square FT)"""
function createSphericalHarmonicsStructure(moleculepdb, LMAX::Integer,KMAX::Integer,rmax::Float64)
    println("Initializing SH structure from $moleculepdb with LMAX=$LMAX, KMAX=$KMAX, rmax=$rmax.")
    atomlist = loadPDB(moleculepdb)

    density = getDensity(atomlist, LMAX, KMAX, rmax)
    fourierTransform = forward(density)
    intensity = absoluteSquare(fourierTransform)

    return density, fourierTransform, intensity
end

"""Loads a pdb file and creates a cubic description of the electron denisty, the Fourier transform and the Intensity (absolute square FT)"""
function createCubicStructure(moleculepdb, cubesize::Integer, rmax::Float64)
  println("Initializing cubic structure from $moleculepdb with cubesize=$cubesize, rmax=$rmax.")
  atomlist = loadPDB(moleculepdb)

  #Calculate optional cubic representation
  densCube = getDensityCube(atomlist, cubesize, rmax)
  fourierCube = forward(densCube)
  intensCube = absoluteSquare(fourierCube)

  return densCube, fourierCube, intensCube
end

"""Expands a cubic volume in a spherical harmonics expansion"""
function cubeToSphericalHarmonics(volume::CubeVolume, KMAX::Int64, LMAX::Int64)
  rmax = volume.rmax
  surflist = Array{Complex{Float64},1}[]
  dr = rmax/KMAX
  for k = 1:KMAX
      r = k*dr
      push!(surflist, getSurface((phi,theta)->getVolumeInterpolated(volume, sphericalToCartesian(r,phi,theta)), LMAX))
  end
  surfVolume = SurfaceVolume(surflist, LMAX, KMAX, rmax)

  return getSphericalHarmonicsVolume(surfVolume)

end
