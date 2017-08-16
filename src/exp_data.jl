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

"""Calculates the two- and three-photon correlation from dense pixelized images"""
function calculate_correlations_in_image(image_list::Array{Array{Float64,2},1}, ksize::Int64, N::Int64=32, filename::String="correlations.dat")
    c2_full = zeros(Float64, N, ksize, ksize)
    c3_full = zeros(Float64, N, N, ksize, ksize, ksize)
    #@showprogress
    @time for (i,image) in enumerate(image_list)
        println("Processing $(i)th image.")
        (sx,sy) = Base.size(image)
        (cx,cy) = (round(Int64, sx/2),round(Int64, sx/2))
        range = 1:sx
        center = Float64[cx,cy, 0.0]
        da = pi/N
        c2_part,c3_part = @sync @parallel ( (a,b) -> (a[1]+b[1], a[2]+b[2])) for x1 = range

            c2_local = zeros(Float64, N, ksize, ksize)
            c3_local = zeros(Float64, N, N, ksize, ksize, ksize)
            println("X1: $x1")
            for y1 = range
                for x2 = range
                    for y2 = range
                        p1 = Float64[x1,y1, 0.0] - center
                        p2 = Float64[x2,y2, 0.0] - center

                        k1 = round(Int64,norm(p1))
                        k2 = round(Int64,norm(p2))

                        if k1 >= k2 && k1 != 0 && k2 != 0 && k1 <= cx && k2 <= cx

                            alpha = mod(angle_between(p1,p2),pi)
                            ai = clamp(floor(Int64, alpha/da) + 1, 1, N)
                            @inbounds c2_local[ai,k1,k2] += real(image[x1,y1]*image[x2,y2]) * doubletFactor(k1,k2) * 1/(k1*k2)

                            for x3 = range
                                for y3 = range
                                    p3 = Float64[x3,y3,0.0] - center
                                    k3 = round(Int64,norm(p3))
                                    if k2 >= k3 && k3 != 0 && k3 <= cx

                                        beta = mod(angle_between(p1,p3),pi)
                                        bi = clamp(floor(Int64, beta/da) + 1, 1, N)
                                        @inbounds c3_local[ai,bi,k1,k2,k3] += real(image[x1,y1]*image[x2,y2]*image[x3,y3]) * tripletFactor(k1,k2,k3) * 1/(k1*k2*k3)

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

        c2_full += c2_part
        c3_full += c3_part
    end
    serializeToFile(filename, ( Dict("num_pictures"=>0, "num_incident_photons"=>0, "qcut"=>1.0, "K"=>ksize, "N"=>N, "dq"=>0.1), sdata(c2_full),sdata(c3_full)))
end
