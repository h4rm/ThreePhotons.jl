export image_to_photons
export precompute_distances_and_angles, calculate_correlations_in_image_integral, symmetrize_image, calculate_correlations_in_image_using_single_photons

"""Transforms a sparse block-image into a list of photons that can be used to calculate correlations"""
function image_to_photons(img::Matrix{UInt32}, qmax::Float64=1.0)
    photons = Vector{Float64}[]
    xsize, ysize = Base.size(img)
    xhalf = floor(xsize/2.0)
    yhalf = floor(ysize/2.0)
    dq = qmax / xhalf
    #     println("Sorting $(sum(abs, img)) photons.")
    for x=1:xsize
        for y=1:ysize
            if img[x,y] > 0 && (x-xhalf,y-yhalf) != (0,0)
                for p in Base.Iterators.repeated(dq*[x-xhalf,y-yhalf,0], 1)#Int64(img[x,y])
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


"""Symmetrizes any image"""
function symmetrize_image(img::Matrix{Float64})
    # img_sym = deepcopy(img)
    # center = [ceil(Base.size(img)[1]/2.0), ceil(Base.size(img)[2]/2.0)]
    # for x=1:Base.size(img)[1]
    #     for y=1:x
    #         pos = [x,y]
    #         mpos = round(Int64,-1.0*([x,y] - center)+center)
    #         m = 0.5*(img_sym[x,y] + img_sym[mpos[1],mpos[2]])
    #         img_sym[x,y] = m
    #         img_sym[mpos[1],mpos[2]] = m
    #     end
    # end
    # return img_sym
    return 0.5*(img+rot180(img))
end

"""Calculates the two- and three-photon correlation from dense pixelized images"""
function calculate_correlations_in_image_integral(image_list::Array{Array{Float64,2},1}, K2::Int64, K3::Int64, N::Int64=32, filename::String="histo.dat", symmetrize::Bool=false)

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

"""Calculates the two- and three-photon correlation from dense pixelized images"""
function calculate_correlations_in_image_using_single_photons(image_list::Array{Array{Float64,2},1}, K2::Int64, K3::Int64, N::Int64, filename::String, overall_maximum::Float64, photons_per_image::Int64, symmetrize::Bool=false)

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

            #Calculating maximum intensity and expected photon counts
            integral = sum(abs,image)
            ppi = integral/overall_maximum * photons_per_image

            image_volume = ImageVolume(image, float(Base.size(image)[1]))

            photon_list,_ = pointsPerOrientation(image_volume, K2, float(K2), photons_per_image, rot=eye(3), lambda=0.0, beamstop_width=0.0, print_warning=false)

            histogramMethod(photon_list, c1, c2, c3, 1.0, N, K2, K3, lambda)

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
