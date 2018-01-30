julia5 -e "using ThreePhotons5; using JLD; data = deserializeFromFile(\"/home/bardenn/Documents/projects/reconstruction/code/data/output_owl/$1\"); save(\"/home/bardenn/Documents/projects/reconstruction/code/data/output_owl/$1.jld\", \"data\", data);
rm(\"/home/bardenn/Documents/projects/reconstruction/code/data/output_owl/$1\")"
