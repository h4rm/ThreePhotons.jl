module ThreePhotons
using Distributions
using Optim
using ProgressMeter

export
	#structures
	atomfactor,
	atomfactorlist,
	atomType,
	loadPDB,
	centerStructure,
	getDensity_cartesian,
	getIntensity_cartesian,
	getFourier_cartesian,
	createSphericalHarmonicsStructure,
	createCubicStructure,
	cubeToSphericalHarmonics,

	#cubic volumes
	CubeVolume,
	dr,
	qmax,
	getDensityCube,
	average,
	calculate_center_of_mass,
	translate_cube,
	center_cube,
	displacement_fitting,
	forward,
	absoluteSquare,
	mirrorCube,
	similarity,
	saveCube,
	loadCube,
	getVolumeInterpolated,
	rotateStructure,
	shell_correlation,
	shell_correlation_ISC,
	shell_correlation_FSC,
	ISC,
	FSC,

	#spherical harmonics
	Volume,
	SphericalVolume,
	SphericalHarmonicsVolume,
	SurfaceVolume,
	wigner,
	sphPlm,
	Plm,
	Ylm,
	get_size,
	num_val,
	num_coeff,
	dr,
	dq,
	qmax,
	average,
	disable_radial_interpolation,
	transform,
	backtransform,
	getSurfaceVolume,
	getSphericalHarmonicsVolume,
	seanindex,
	getc,
	setc,
	cvec_get,
	cvec_set,
	i_to_ang,
	ang_to_i,
	getSurface,
	getDensity,
	densityShell,
	fourierShell,
	intensityShell,
	deleteTerms,
	sphBesselTransform,
	forward,
	backward,
	absoluteSquare,
	shell_correlation,
	shell_correlation_ISC,
	shell_correlation_FSC,
	ISC,
	FSC,
	similarity,
	Umat,
	comp_to_real,
	real_to_comp,
	getVolumeInterpolated,
	getCube,
	saveCube,
	d_work,
	Dlms,
	rotationMatrix,
	rotateStructure,
	get_value,
	negativityCheck,

	#structure determination
	randomStartStructure,
	randomPertubation,
	rotation_search,
	get_MC_step,
	calculate_temperature,
	saveState,
	regular_action,
	rotate_all_at_once,
	rotate_hierarchical,
	increase_stepsize,
	decrease_stepsize,
	evaluate_step,
	checkRotationSearch,

	#phasing
	sRAAR,

	#data processing
	fitStructures_full,
	fitStructures_random,
	getAllEvenRotations,
	sRAAR_bestfit,
	calculateSC,
	postprocess_run,
	get_qrange,
	calculate_cutoff,
	calculate_maximum_resolution,
	load_runs,
	analyse_ensemble,
	load_fitted_intensity,
	load_density,
	average_intensities,
	calculate_optimal,
	load_parameter_list,
	get_parameter_entry,
	get_optimal,
	load_sc_vs_triplets,
	substract_gamma,
	radial_difference_with_intensity,
	denoise_structure,

	#plotting
	plot_random_2photon_slices,
	plot_random_3photon_slices,
	plotPoints2D,
	plotPointCloud,
	plot_scattering_image,
	compare_c2_grid,
	compare_histogram_with_theory,
	initialize_CUDA

import Base.real, Base.abs, Base.sumabs, Base.+, Base.-, Base.*, Base./, Base.length

include("utilities.jl")
include("spherical_harmonics.jl")
include("cubic.jl")
include("structure.jl")
include("correlations.jl")
include("datagen.jl")
include("phases.jl")
include("data_processing.jl")
include("determination.jl")
include("exp_data.jl")
# include("cuda.jl")

"""
A Julia package for the structure determination from single molecule X-ray scattering experiments at very low photon counts.
"""

end #module
