"""Plots random 2 photon correlation slices"""
function plot_random_2photon_slices(c2list::Dict; normalization=false, list=[random_doublet(25) for i = 1:20])
    for i = 1:length(list)
        k1 = list[i][1]
        k2 = list[i][2]

        figure(i)
        for name in keys(c2list)
          slice = real(c2list[name][:,k2,k1])

          if normalization slice /= sumabs(slice) end
          plot(slice,"-o", label=name)
        end

        title(latexstring("2 Photon Correlation \$k_1=$k1, k_2=$k2\$"))
        xlabel(L"$\alpha$ [rad]")
        ylabel(latexstring("c_2($k1,$k2)") )
        legend()
    end
end

"""Plots random 3 photon correlation slices"""
function plot_random_3photon_slices(c3list::Dict; list=[random_triplet(25) for i = 1:20], integrate=false,  normalization=false)

  for i = 1:length(list)
      k1 = list[i][1]
      k2 = list[i][2]
      k3 = list[i][3]

      figure(2*i)
      for (j,name) in enumerate(keys(c3list))
        slice = real(c3list[name][:,:,k3,k2,k1])
        if normalization slice /= sumabs(slice) end

        #Plot slices
        set_cmap("hot")
        subplot(1,length(keys(c3list)),j)
        imshow(slice, extent=[0, pi, 0, pi])
        clim(minimum(slice), maximum(slice))
        title(latexstring("$name\n\$k_1=$k1\,k_2=$k2\,k_3=$k3\$"))
        xlabel(L"\alpha")
        ylabel(L"\beta")
        grid(false)
      end

      figure(2*i+1)
      for (j,name) in enumerate(keys(c3list))
        slice = real(c3list[name][:,:,k3,k2,k1])
        if normalization slice /= sumabs(slice) end

        #plot diagonales
        plot( diag(slice), "-o", label=name)
      end

      xlabel(L"\alpha=\beta [rad]")
      ylabel("t")
      title(latexstring("3 Photon Correlation \$k_1=$k1 k_2=$k2 k_3=$k3\$"))
      legend()

    end
end

"""Plots an array of 2D-points as 2D scatter plot"""
function plotPoints2D(points)
  scatter( [points[i][1] for i=1:length(points)], [points[i][2] for i=1:length(points)])
  # xlim(-qmax, qmax)
  # ylim(-qmax, qmax)
end

"""Plots an array of 3D-points as 3D scatter plot"""
function plotPointCloud(points)
  # points = []
  # for i = 1:4000
  #     p,rot = pointsPerOrientation(200)
  #     p = map((x)->rot*x, p)
  #     append!(points, p)
  # end

  # colors = [qmax/sumabs(points[i]) for i = 1:length(points)]
  clf()
  scatter3D( [points[i][1] for i=1:length(points)], [ points[i][2]for i=1:length(points)], [ points[i][3]for i=1:length(points)], s=0.1)
end

"""From an intensity cube, generates scattering data and plots it"""
#noise::Noise=GaussianNoise(0.0, 0.0, false)
function plot_scattering_image(intensity::CubeVolume; number_incident_photons::Integer=10000, qcut_ratio::Float64=0.83, point_size::Float64=50.0, colorfill="red", coloredge="black")

    local qmax = intensity.rmax*qcut_ratio
    # noisy_intensity = get_noisy_intensity(intensity, noise)
    noisy_intensity = intensity
    p,rot = pointsPerOrientation(noisy_intensity, qmax,qmax/3.0, number_incident_photons)

    #Plot underlying intensity
    r = linspace(-qmax, qmax, noisy_intensity.cubesize)
    myslice = Float64[getVolumeInterpolated(noisy_intensity, rot*[-x,y,0]) for x=r,y=r]
    myslice = (myslice).^(0.2)
    # myslice = log(max(myslice, 1e-6))
    # myslice = log(myslice)
    fig = imshow(myslice, interpolation="hermite", extent=[-qmax, qmax, -qmax, qmax], cmap="Blues")

    #Plot scattering photons
    scatter( [p[i][2] for i=1:length(p)], [ p[i][1] for i=1:length(p)], c=colorfill, s=point_size, alpha=1.0, edgecolors=coloredge)#

    #plot center
    scatter([0.0], [0.0], marker="+", alpha=1.0, s=100.0, color="red")

    fig[:axes][:get_xaxis]()[:set_visible](false)
    fig[:axes][:get_yaxis]()[:set_visible](false)
end

function compare_c2_grid(a::C2, b::C2)
    _,_,K = Base.size(a)
    ac = complete_two_photon_correlation(a)
    bc = complete_two_photon_correlation(b)

    ratios = Float64[100.0*mean(ac[:,k2,k1]-bc[:,k2,k1]).^2/sumabs(ac[:,k2,k1]) for k1 = 1:K, k2=1:K]
    imshow(ratios, interpolation="none", extent=[1,K,K,1])
    title("Average Deviation [%]")
    colorbar()
    println("Difference: $(100.0*sumabs((ac-bc).^2)/sumabs(ac))")
end

"""Compares a list of histograms visually"""
function compare_histogram_with_theory(histograms::Dict=Dict(), N::Int64=32, K::Int64=16, normalization::Bool=false, ctype::String="c3", volumes::Dict=Dict(), L::Int64=10)
    histolist_c3 = Dict()
    histolist_c2 = Dict()
    for (name,file) in histograms
        _,c2,_,c3 = loadHistograms(K,file)
        histolist_c3[name] = c3
        histolist_c2[name] = c2
    end
    basis = complexBasis(L,N,25)

    if ctype=="c3"
        c3list = Dict( name=>FullCorrelation_parallized(volume, basis, K, true, true) for (name,volume) in volumes)
        plot_random_3photon_slices(merge(histolist_c3,c3list); list=[random_triplet(K) for i = 1:20], normalization=normalization)
    else
        c2list = Dict( name=>twoPhotons(volume, basis, K, true, true) for (name,volume) in volumes)
        plot_random_2photon_slices(merge(histolist_c2,c2list), normalization=normalization, list=[random_doublet(K) for i = 1:20])
    end
end
