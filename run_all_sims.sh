threads=7

julia --threads $threads stochtree.jl nsites=75 only_simulate=true ntrees=28 


mkdir trees

mv *.csv trees
mv *.fasta trees
mv *.tre trees

mkdir 500_iteration_simulation
mkdir 1000_iteration_simulation
cp -r trees 500_iteration_simulation
cp -r trees 1000_iteration_simulation

julia --threads $threads stochtree.jl tree_directory=trees directory=500_iteration_simulation iters=500
julia --threads $threads stochtree.jl tree_directory=trees directory=1000_iteration_simulation iters=1000 

rm -rf trees
mkdir trees


mkdir dependent_simulation


julia --threads $threads stochtree.jl nsites=75 nsites_big=300 only_simulate=true ntrees=28 simulation_type=dependent

mv *.csv trees
mv *.fasta trees
mv *.tre trees

cp -r trees dependent_simulation

julia --threads $threads stochtree.jl tree_directory=trees directory=dependent_simulation iters=1000 simulation_type=dependent

julia stochtree.jl simulation_type=time_benchmark