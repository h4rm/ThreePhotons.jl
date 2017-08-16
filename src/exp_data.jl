export image_to_photons
export precompute_distances_and_angles, calculate_correlations_in_image

"""Transforms a sparse block-image into a list of photons that can be used to calculate correlations"""
function image_to_photons(img::Matrix{UInt32}, qmax::Float64=1.0)
    photons = Vector{Float64}[]
    xsize, ysize = Base.size(img)
    xhalf = floor(xsize/2.0)
    yhalf = floor(ysize/2.0)
    dq = qmax / xhalf
#     println("Sorting $(sumabs(img)) photons.")
    for x=1:xsize
        for y=1:ysize
            if img[x,y] > 0 && (x-xhalf,y-yhalf) != (0,0)
                for p in repeated(dq*[x-xhalf,y-yhalf,0], 1)#Int64(img[x,y])
                    push!(photons, p)
                end
            end
        end
    end
    return photons
end

############################################################################
## EXPERIMENTAL DATA
############################################################################

function precompute_distances_and_angles(image_width::Int64, N::Int64)
    distances = zeros(Int64, image_width, image_width)
    angles = zeros(Int64, image_width, image_width, image_width, image_width)

    range = 1:image_width
    center = Float64[round(Int64, image_width/2),round(Int64, image_width/2), 0.0]
    da = pi/N
    for x1 = range
        for y1 = range
            p1 = Float64[x1,y1, 0.0] - center
            @inbounds distances[x1,y1] = round(Int64,norm(p1 - center))
            for x2 = range
                for y2 = range
                    p2 = Float64[x2,y2, 0.0] - center
                    angle = mod(angle_between(p1,p2),pi)
                    @inbounds angles[x1, y1, x2, y2] = clamp(floor(Int64, angle/da) + 1, 1, N)
                end
            end
        end
    end
    return distances, angles
end

"""Calculates the two- and three-photon correlation from dense pixelized images"""
function calculate_correlations_in_image(image_list::Array{Array{Float64,2},1}, K2::Int64, K3::Int64, N::Int64=32, filename::String="correlations.dat")

    (sx,sy) = Base.size(image_list[1])
    (cx,cy) = (round(Int64, sx/2),round(Int64, sx/2))
    range = 1:sx
    center = Float64[cx,cy, 0.0]
    da = pi/N

    distances,angles = precompute_distances_and_angles(sx, N)

    c2_full,c3_full = @sync @parallel ( (a,b) -> (a[1]+b[1], a[2]+b[2])) for i=1:length(image_list)
        image = image_list[i]
        println("Processing image #$(i)")
        c2_local = zeros(Float64, N, K2, K2)
        c3_local = zeros(Float64, N, N, K3, K3, K3)
        for x1 = range
            for y1 = range
                @inbounds k1 = distances[x1,y1]

                for x2 = range
                    for y2 = range
                        @inbounds k2 = distances[x2,y2]

                        if k1 >= k2 && k1 != 0 && k2 != 0 && k1 <= K2 && k2 <= K2
                            @inbounds ai = angles[x1,y1,x2,y2]
                            @inbounds c2_local[ai,k1,k2] += real(image[x1,y1]*image[x2,y2]) * doubletFactor(k1,k2) * 1/(k1*k2)

                            for x3 = range
                                for y3 = range
                                    @inbounds k3 = distances[x3,y3]

                                    if k2 >= k3 && k3 != 0 && k1 <= K3 && k2<= K3 && k3 <= K3
                                        bi = angles[x1,y1,x3,y3]
                                        @inbounds c3_local[ai,bi,k1,k2,k3] += real(image[x1,y1]*image[x2,y2]*image[x3,y3]) * tripletFactor(k1,k2,k3) * 1/(k1*k2*k3)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        flush(STDOUT)
        (c2_local, c3_local)
    end
    serializeToFile(filename, ( Dict("num_pictures"=>0, "num_incident_photons"=>0, "qcut"=>1.0, "K2"=>K2, "K3"=>K3, "N"=>N, "dq"=>0.1), sdata(c2_full),sdata(c3_full)))
end
