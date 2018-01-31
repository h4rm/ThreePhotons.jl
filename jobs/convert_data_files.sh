#This script converts old julia0.5 serialized files to JLD and then loads the JLD in julia0.6 and converts them back to julia0.6 format
julia5 -e "using ThreePhotons5; using JLD; data = deserializeFromFile(\"/home/bardenn/Documents/projects/reconstruction/code/data/output_owl/$1\"); save(\"/home/bardenn/Documents/projects/reconstruction/code/data/output_owl/$1.jld\", \"data\", data);"
#rm(\"/home/bardenn/Documents/projects/reconstruction/code/data/output_owl/$1\")

julia -e "using ThreePhotons; using JLD; data = deserializeFromFile(\"/home/bardenn/Documents/projects/reconstruction/code/data/output_owl/$1.jld\"); serializeToFile(\"/home/bardenn/Documents/projects/reconstruction/code/data/output_owl/$1\", data);"
