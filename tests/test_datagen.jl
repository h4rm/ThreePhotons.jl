#Data generation Cubic
batchsize = Integer(1e4)
generateHistogram(intensityCube; qcut=intensityCube.rmax, K=KMAX, N=N, max_pictures=batchsize,max_triplets=0, batchsize = batchsize, number_incident_photons=calculate_incident_photons(10), numprocesses=1, file="/tmp/correlations.dat", noise=GaussianNoise(0.1, 0.5, 14))
loadHistograms(8, "/tmp/correlations.dat")

#Data generation spherical harmonics
volume = getSurfaceVolume(intensity)
volume.radial_interp = false
try rm("/tmp/correlations.dat") end
@profile generateHistogram(volume; qcut=intensity.rmax, K=KMAX, N=N, max_pictures=batchsize,max_triplets=0, batchsize = batchsize, number_incident_photons=calculate_incident_photons(10), numprocesses=nworkers(), file="/tmp/correlations.dat", noise=GaussianNoise(0.1, intensity.rmax/4.0, 10))
loadHistograms(8, "/tmp/correlations.dat")
