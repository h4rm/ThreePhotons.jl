using Compat
using ThreePhotons

include("runs.jl")

directory_list = ["final_res_vs_pictures", "paper_res_vs_pictures"]

is_run = (root) -> (true in map((x)->contains(root,x), directory_list))

"""Checks if there is core dump of a job"""
function check_crash(delete::Bool=false)
	deletestring = (delete ? " -delete" : "")
	file = "core.*"
	run(`find . -name $file$deletestring`)
end

"""Prints the result of a run given by `dir`"""
function print_result(dir::String)
	det = "$dir/det.out"
	termout = "$dir/terminal-out"
	statefilepath = "$dir/state.dat"

	if isfile(statefilepath)
		state = Dict()
		try
			state = deserializeFromFile(statefilepath)
		catch
			println("Error loading $statefilepath")
		end

		if haskey(state, "state")
			if state["state"] != "finished"
				print_with_color(:red, "$dir\t$(state["state"])\n")
			elseif state["state"] == "finished"
				print_with_color(:green, "$dir\t$(state["state"])\n")
			end
		end

		if isfile(det)
			tail = readstring(`tail -n 1 $det`)
			println("\t$tail")
		end
		if isfile(termout)
			out = readstring(`tail -n200 $termout`)

			intensity_correlation = match(r"Best intensity correlation: (.+)", out)
			if intensity_correlation != nothing && length(intensity_correlation[1]) < 1000
				print_with_color(:green, "\tBest intensity correlation:  $(intensity_correlation[1])\n")
			end

			density_correlation = match(r"Best density correlation: (.+)", out)
			if density_correlation != nothing && length(density_correlation[1]) < 1000
				print_with_color(:green, "\tBest density correlation:  $(density_correlation[1])\n")
			end

			resolution = match(r"Best resolution: (.+)", out)
			if resolution != nothing && length(resolution[1]) < 1000
				print_with_color(:green, "\tBest resolution:  $(resolution[1])\n")
			end
		end
	else
		print_with_color(:grey, "$dir:\tNo state file\n")
	end

end

"""Plots a summary of a certain job"""
function monitor(dir::String, gnuplot::Bool=false)
	gnuplot = []
	for subdir in readdir(dir)
		println("---------------------------------------------------------")
		if isdir("$dir/$subdir")
			print_result("$dir$subdir")
		end
	end

	if ENV["LOGNAME"] != "bardenn" && gnuplot == true
		list = join(gnuplot,",")
		script = """set title 'Energy vs. 5000 last steps';plot $list"""
		run(`gnuplot -p -e "$script" `)
	end
end


"""Iterates over all subdirectories in parallel and displays 'nonfinished' runs"""
function check_finished(root_dir::String = "$ENV_root")
  for (root, dirs, files) in Compat.walkdir(root_dir)

    if ("terminal-out" in files)# && is_run(root)
			println("--------------------------------------")
      print_result(root)
    end
  end
end

"""Goes over a list of defined directories and checks if job is finished, if not - resubmits"""
function resubmit_unfinished_jobs(root_dir::String = "$ENV_root")

	for (root, dirs, files) in Compat.walkdir(root_dir)

    if ("job.sh" in files) #&& is_run(root)
			state = Dict()
			try state = deserializeFromFile("$root/state.dat") end

			if !haskey(state, "state") || state["state"] != "finished"
				# println("Submit job $root/job.sh")
				submit_job("$root/job.sh")
			end
    end
  end
end

function recalculate_sc(root_dir::String, K::Int64)
	densityCube,fourierCube,intensityCube = createCubicStructure("$(ENV["DETERMINATION_PATH"])/structures/crambin.pdb", 2*K+1, float(K))

	for (root, dirs, files) in Compat.walkdir(root_dir)

    if "density.mrc" in files
			try
				state = deserializeFromFile("$root/state.dat")
				isc,old_fsc,isc_nofitting = state["sc"]
				fsc = FSC(loadCube("$root/density.mrc"), fourierCube)
				replace_NaN!(fsc)
				state["sc"] = (isc,fsc,isc_nofitting)
				serializeToFile("$root/state.dat", state)

				println("Succesfully reprocessed $root")
				println("\tOld resolution: ", calculate_maximum_resolution(old_fsc, dr(fourierCube)))
				println("\tNew resolution: ", calculate_maximum_resolution(fsc, dr(fourierCube)))

			catch
				println("Failed to process $root")
			end
    end
  end
end

if length(ARGS)>0 && length(ARGS[1])>0
	check_finished(ARGS[1])
	exit()
end
