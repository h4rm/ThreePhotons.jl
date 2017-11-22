export image_to_photons
export precompute_distances_and_angles, calculate_correlations_in_image, symmetrize_image

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
    center = ceil(Float64[image_width/2.0,image_width/2.0, 0.0])
    da = pi/N
    for x1 = range
        for y1 = range
            p1 = Float64[x1,y1, 0.0] - center
            @inbounds distances[x1,y1] = round(Int64,norm(p1))
            for x2 = range
                for y2 = range
                    p2 = Float64[x2,y2, 0.0] - center

                    #calculate angles
                    alpha = angle_between(p1,p2)
                    @inbounds angles[x1, y1, x2, y2] = clamp(floor(Int64, alpha/da) + 1, 1, 2*N)
                end
            end
        end
    end
    return distances, angles
end

"""Calculates the two- and three-photon correlation from dense pixelized images
This version is depricated because it produces pixel artifacts
"""
function calculate_correlations_in_image_old(image_list::Array{Array{Float64,2},1}, K2::Int64, K3::Int64, N::Int64=32, filename::String="histo.dat")

    (sx,sy) = Base.size(image_list[1])
    range = 1:sx
    da = pi/N

    distances,angles = precompute_distances_and_angles(sx, N)

    c1_full = zeros(Float64, K2)
    c2_full = zeros(Float64, N, K2, K2)
    c3_full = zeros(Float64, N, 2*N, K3, K3, K3)
    number_analyzed_images = 0

    for j=1:ceil(Int64, length(image_list)/nworkers())
        println("Processing batch $j")
        #Processing next batch of nworkers() images
        c1_part,c2_part,c3_part = @sync @parallel ( (a,b) -> (a[1]+b[1], a[2]+b[2], a[3]+b[3])) for i=((j-1)*nworkers()+1):clamp(j*nworkers()+1, 1, length(image_list))
            image = image_list[i]
            println("Processing image #$(i)")
            c1_local = zeros(Float64, K2)
            c2_local = zeros(Float64, N, K2, K2)
            c3_local = zeros(Float64, N, 2*N, K3, K3, K3)

            c1_counts = zeros(Float64, K2)
            c2_counts = zeros(Float64, N, K2, K2)
            c3_counts = zeros(Float64, N, 2*N, K3, K3, K3)
            for x1 in range
                for y1 in range
                    @inbounds k1 = distances[x1,y1]

                    if k1 > 0 && k1 <= K2
                        @inbounds c1_local[k1] += real(image[x1,y1])*(1/k1)
                        @inbounds c1_counts[k1] += 1.0

                        for x2 in range
                            for y2 in range
                                @inbounds k2 = distances[x2,y2]

                                if k2 <= k1  && k2 > 0 && k2 <= K2
                                    @inbounds ai = angles[x1,y1,x2,y2]
                                    ais = ai
                                    if ai > N ais = 2*N-ais+1 end
                                    @fastmath val2 = real(image[x1,y1]*image[x2,y2]) * doubletFactor(k1,k2) * 1/(k1*k2)
                                    @inbounds c2_local[ais,k2,k1] += val2
                                    @inbounds c2_counts[ais,k2,k1] += 1.0

                                    if k1 <= K3 && k2<= K3
                                        for x3 in range
                                            for y3 in range
                                                @inbounds k3 = distances[x3,y3]

                                                if k3 <= k2 && k3 > 0 && k3 <= K3
                                                    @inbounds bi = angles[x1,y1,x3,y3]
                                                    bis = bi
                                                    if ai > N bis = 2*N-bis+1 end

                                                    @fastmath val3 = real(image[x1,y1]*image[x2,y2]*image[x3,y3]) * tripletFactor(k1,k2,k3) * 1/(k1*k2*k3)
                                                    @inbounds c3_local[ais,bis,k3,k2,k1] += val3
                                                    @inbounds c3_counts[ais,bis,k3,k2,k1] += 1.0
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            flush(STDOUT)
            (c1_local ./ max(c1_counts, 1.0), c2_local ./ max(c2_counts, 1.0), c3_local ./ max(c3_counts, 1.0))
        end
        c1_full += c1_part
        c2_full += c2_part
        c3_full += c3_part
        number_analyzed_images += nworkers()

        serializeToFile(filename, ( Dict("num_pictures"=>number_analyzed_images, "num_incident_photons"=>0, "qcut"=>1.0, "K2"=>K2, "K3"=>K3, "N"=>N, "dq"=>0.1), c2_full, c3_full, c1_full))
        flush(STDOUT)
    end
end

"""Symmetrizes any image"""
function symmetrize_image(img::Matrix{Float64})
    img_sym = deepcopy(img)
    center = [ceil(Base.size(img)[1]/2.0), ceil(Base.size(img)[2]/2.0)]
    for x=1:Base.size(img)[1]
        for y=1:x
            pos = [x,y]
            mpos = round(Int64,-1.0*([x,y] - center)+center)
            m = 0.5*(img_sym[x,y] + img_sym[mpos[1],mpos[2]])
            img_sym[x,y] = m
            img_sym[mpos[1],mpos[2]] = m
        end
    end
    return img_sym
end

"""Calculates the two- and three-photon correlation from dense pixelized images"""
function calculate_correlations_in_image(image_list::Array{Array{Float64,2},1}, K2::Int64, K3::Int64, N::Int64=32, filename::String="histo.dat", symmetrize::Bool=false)

    da = pi/N
    (image_width,sy) = Base.size(image_list[1])
    center = ceil(Float64[image_width/2.0,image_width/2.0])
    rotations = [[cos(phi) sin(phi); -sin(phi) cos(phi)] for phi = 0:da:2*pi]
    up = Float64[0.0, 1.0]

    c1_full = zeros(Float64, K2)
    c2_full = zeros(Float64, N, K2, K2)
    c3_full = zeros(Float64, N, 2*N, K3, K3, K3)
    number_analyzed_images = 0

    for j=1:ceil(Int64, length(image_list)/nworkers())
        println("Processing batch $j")
        #Processing next batch of nworkers() images
        c1_part,c2_part,c3_part = @sync @parallel ( (a,b) -> (a[1]+b[1], a[2]+b[2], a[3]+b[3])) for i=((j-1)*nworkers()+1):clamp(j*nworkers()+1, 1, length(image_list))
            image = image_list[i]
            if symmetrize == true
                image = symmetrize_image(image)
            end
            println("Processing image #$(i)")
            c1_local = zeros(Float64, K2)
            c2_local = zeros(Float64, N, K2, K2)
            c3_local = zeros(Float64, N, 2*N, K3, K3, K3)

            for k1 = 1:K2
                write(STDOUT,"$k1 ")
                p1 = k1*up

                sum1 = 0.0
                for rot in rotations
                    rp1 = round(Int64,rot*p1 + center)
                    sum1 += image[rp1[1], rp1[2]]
                end
                c1_local[k1] = sum1 * 1/(k1)

                for k2 = 1:k1
                    for alpha = 1:2*N
                        a = Float64(da*alpha-da/2)
                        p2 = k2*[cos(a) sin(a); -sin(a) cos(a)]*up
                        ais = alpha
                        if alpha > N ais = 2*N-ais+1 end

                        sum2 = 0.0
                        for rot in rotations
                            rp1 = round(Int64,rot*p1 + center)
                            rp2 = round(Int64,rot*p2 + center)
                            sum2 += image[rp1[1], rp1[2]]*image[rp2[1], rp2[2]]
                        end
                        c2_local[ais,k2,k1] = sum2 * 1/(k1*k2) * doubletFactor(k1,k2)

                        for k3=1:k2
                            for beta=1:2*N

                                if k1 <= K3 && k2 <= K3 && k3 <= K3
                                    b = Float64(da*beta-da/2)
                                    p3 = k3*[cos(b) sin(b); -sin(b) cos(b)]*up
                                    bis = beta
                                    if alpha > N bis = 2*N-bis+1 end

                                    sum3 = 0.0
                                    for rot in rotations
                                        rp1 = round(Int64,rot*p1 + center)
                                        rp2 = round(Int64,rot*p2 + center)
                                        rp3 = round(Int64,rot*p3 + center)
                                        sum3 += image[rp1[1], rp1[2]]*image[rp2[1], rp2[2]]*image[rp3[1], rp3[2]]
                                    end
                                    c3_local[ais,bis,k3,k2,k1] = sum3 * 1/(k1*k2*k3) * tripletFactor(k1,k2,k3)
                                end
                            end
                        end
                    end
                end
            end
            (c1_local, c2_local, c3_local)
        end
        c1_full += c1_part
        c2_full += c2_part
        c3_full += c3_part
        number_analyzed_images += nworkers()

        serializeToFile(filename, ( Dict("num_pictures"=>number_analyzed_images, "num_incident_photons"=>0, "qcut"=>1.0, "K2"=>K2, "K3"=>K3, "N"=>N, "dq"=>0.1), c2_full, c3_full, c1_full))
        flush(STDOUT)
    end
end
