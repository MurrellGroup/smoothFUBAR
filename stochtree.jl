using CodonMolecularEvolution, MolecularEvolution, FASTX, Distributions, EvalMetrics, Plots, DataFrames, CSV



argument_dictionary = Dict(
    "directory" => ".",
    "positive_sites" => "1",
    "ntrees" => "30",
    "ntaxa" => "50",
    "nsites" => "500",
    "nsites_big" => "2500",
    "K" => "30",
    "iters" => "1000",
    "tree_directory" => "",
    "only_simulate" => "false",
    "simulation_type" => "independent")

if length(ARGS) == 0

    println("You have not entered any arguments. Do you wish to proceed with default values: ")
    for keyvalue in argument_dictionary
        println(keyvalue)
    end
    println("Yes/No? [Y/N]")
    yes_no_prompt = readline()

    if yes_no_prompt == "N"
        println("To modify default values, pass command line arguments of the form [VARIABLE_NAME]=[DESIRED_VALUE]")
        exit()
    end

end

for argument in ARGS
    key, value = (split(argument, "=")[1], split(argument, "=")[2])
    argument_dictionary[key] = value
end

cd(argument_dictionary["directory"])




function simulate_single_positive_site(nsites, σα, rate, multiplier, pos_sites)

    negative_alphavec::Vector{Float64} = []
    negative_betavec::Vector{Float64} = []

    while length(negative_alphavec) < nsites

        alpha = exp(randn() * σα)
        beta = rand(Exponential(rate))

        if beta < alpha
            push!(negative_alphavec, alpha)
            push!(negative_betavec, beta)
        end

    end

    positive_alphavec::Vector{Float64} = [exp(randn() * σα)]
    positive_betavec::Vector{Float64} = [rand(Exponential(rate))]
    while sum(positive_betavec .> positive_alphavec) < pos_sites || length(positive_alphavec) < nsites
        alpha = exp(randn() * σα)
        beta = rand(Exponential(rate))
        if length(positive_alphavec) < nsites
            push!(positive_alphavec, alpha)
            push!(positive_betavec, beta)
        else
            insertion_point = rand(1:nsites)
            positive_alphavec[insertion_point] = alpha
            positive_betavec[insertion_point] = beta
        end
    end

    return negative_alphavec, negative_betavec, positive_alphavec, positive_betavec

end

function simulate_dependent_trees(nsites_small, nsites_big, σα, rate, multiplier, pos_sites)

    negative_small_alphavec, negative_small_betavec, positive_small_alphavec, positive_small_betavec = simulate_single_positive_site(nsites_small, σα, rate, multiplier, pos_sites)

    positive_big_alphavec = copy(positive_small_alphavec)
    positive_big_betavec = copy(positive_small_betavec)

    while length(positive_big_alphavec) < nsites_big
        alpha = exp(randn() * σα)
        beta = rand(Exponential(rate))

        if beta < alpha
            push!(positive_big_alphavec, alpha)
            push!(positive_big_betavec, beta)
        end
    end

    return positive_small_alphavec, positive_small_betavec, positive_big_alphavec, positive_big_betavec

end

function generate_independent_trees(argument_dictionary)
    simulated_negative_trees = []
    simulated_positive_trees = []

    positive_sites = parse(Int64, argument_dictionary["positive_sites"])

    ntrees = parse(Int64, argument_dictionary["ntrees"])
    ntaxa = parse(Int64, argument_dictionary["ntaxa"])
    nsites = parse(Int64, argument_dictionary["nsites"])
    

    Threads.@threads for i in 1:ntrees

        local_negative_tree = standard_tree_sim(ntaxa)

        negative_alphavec, negative_betavec, positive_alphavec, positive_betavec = simulate_single_positive_site(nsites, 0.2, 0.3, 2, positive_sites)

        negative_alphabeta_dataframe = DataFrame(alpha = negative_alphavec, beta = negative_betavec)
        positive_alphabeta_dataframe = DataFrame(alpha = positive_alphavec, beta = positive_betavec)

        CSV.write(string(i)*"_negative_alphabeta.csv", negative_alphabeta_dataframe)
        CSV.write(string(i)*"_positive_alphabeta.csv", positive_alphabeta_dataframe)


        negative_nucs, negative_seqnams, negative_tre = sim_alphabeta_seqs(negative_alphavec, negative_betavec,
            local_negative_tree, outpath=string(i) * "_negative_simdata_FUBAR",
            CodonMolecularEvolution.demo_nucmat, CodonMolecularEvolution.demo_f3x4) #Note the nucmat and f3x4, not from real data

        push!(simulated_negative_trees, (negative_seqnams, negative_nucs, negative_tre, negative_alphavec, negative_betavec))


        local_positive_tree = standard_tree_sim(ntaxa)

        single_strong_index = rand(1:length(positive_alphavec))

        positive_betavec[single_strong_index] = 2 * positive_alphavec[single_strong_index]

        positive_nucs, positive_seqnams, positive_tre = sim_alphabeta_seqs(positive_alphavec, positive_betavec,
            local_positive_tree, outpath=string(i) * "_positive_simdata_FUBAR",
            CodonMolecularEvolution.demo_nucmat, CodonMolecularEvolution.demo_f3x4) #Note the nucmat and f3x4, not from real data

        push!(simulated_positive_trees, (positive_seqnams, positive_nucs, positive_tre, positive_alphavec, positive_betavec))

    end
    return simulated_negative_trees, simulated_positive_trees
end

function read_independent_trees(argument_dictionary)

    tree_directory = argument_dictionary["tree_directory"]

    tree_files = readdir(tree_directory)

    tree_count = Int64(length(tree_files) / 6)

    simulated_positive_trees = []
    simulated_negative_trees = []

    for i in 1:tree_count

        positive_seqnams, positive_nucs = read_fasta(tree_directory*"/"*string(i)*"_positive_simdata_FUBAR.fasta")
        positive_tre = read_newick_tree(tree_directory*"/"*string(i)*"_positive_simdata_FUBAR.tre")

        negative_seqnams, negative_nucs = read_fasta(tree_directory*"/"*string(i)*"_negative_simdata_FUBAR.fasta")
        negative_tre = read_newick_tree(tree_directory*"/"*string(i)*"_negative_simdata_FUBAR.tre")

        positive_alphabeta_frame = CSV.read(tree_directory*"/"*string(i)*"_positive_alphabeta.csv", DataFrame)
        negative_alphabeta_frame = CSV.read(tree_directory*"/"*string(i)*"_negative_alphabeta.csv", DataFrame)

        positive_alphavec, positive_betavec = (positive_alphabeta_frame.alpha, positive_alphabeta_frame.beta)
        negative_alphavec, negative_betavec = (negative_alphabeta_frame.alpha, negative_alphabeta_frame.beta)

        push!(simulated_positive_trees, (positive_seqnams, positive_nucs, positive_tre, positive_alphavec, positive_betavec))
        push!(simulated_negative_trees, (negative_seqnams, negative_nucs, negative_tre, negative_alphavec, negative_betavec))

    end
    
    return simulated_negative_trees, simulated_positive_trees


end

function run_independent_tree_simulation(argument_dictionary)

    simulated_negative_trees = []
    simulated_positive_trees = []

    if argument_dictionary["tree_directory"] == ""
        simulated_negative_trees, simulated_positive_trees = generate_independent_trees(argument_dictionary)
    else 
        simulated_negative_trees, simulated_positive_trees = read_independent_trees(argument_dictionary)
    end
    ntrees = length(simulated_negative_trees)
    K = parse(Int64, argument_dictionary["K"])
    iters = parse(Int64, argument_dictionary["iters"])

    targets = zeros(2 * ntrees)
    scores = zeros(2 * ntrees)

    println("Running with: " * string(Threads.nthreads()) * " threads.")
    println(length(simulated_negative_trees))

    result_tuples = []

    Threads.@threads for i in eachindex(simulated_negative_trees)

        negative_seqnams, negative_nucs, negative_tre, negative_alphavec, negative_betavec = simulated_negative_trees[i]
        positive_seqnams, positive_nucs, positive_tre, positive_alphavec, positive_betavec = simulated_positive_trees[i]

        # Perform the init_2_grid calculation only once 

        f_grid_negative = CodonMolecularEvolution.alphabetagrid(negative_seqnams, negative_nucs, newick(negative_tre))
        f_grid_positive = CodonMolecularEvolution.alphabetagrid(positive_seqnams, positive_nucs, newick(positive_tre))


        dispatch = CodonMolecularEvolution.FUBARweightedpos()

        smooth_negative_result_df, smooth_negative_tuple = smoothFUBAR(dispatch, f_grid_negative, string(i) * "_negative_smooth"; K=K, HMC_samples=iters, plots=false)
        smooth_positive_result_df, smooth_positive_tuple = smoothFUBAR(dispatch, f_grid_positive, string(i) * "_positive_smooth"; K=K, HMC_samples=iters, plots=false)



        dirichlet_negative_result_df, dirichlet_negative_θ = FUBAR(f_grid_negative, string(i) * "_negative_dirichlet", plots=false)
        dirichlet_positive_result_df, dirichlet_positive_θ = FUBAR(f_grid_positive, string(i) * "_positive_dirichlet", plots=false)



        scores[2*i-1] = smooth_negative_tuple.global_posterior_probability
        scores[2*i] = smooth_positive_tuple.global_posterior_probability
        targets[2*i] = 1
        push!(result_tuples, (smooth_negative_result_df, smooth_positive_result_df, dirichlet_negative_result_df, dirichlet_positive_result_df, negative_alphavec, negative_betavec, positive_alphavec, positive_betavec))

        write(string(i) * "_negative_mixiness.txt", string(smooth_negative_tuple.mixing))
        write(string(i) * "_positive_mixiness.txt", string(smooth_positive_tuple.mixing))

    end


    alignment_wide_result_dataframe = DataFrame(assigned_probability=scores, ground_truth=targets)
    CSV.write("alignment_wide_results.csv", alignment_wide_result_dataframe)

    prplot(targets, scores)
    savefig("alignment_wide_prc.svg")
    rocplot(targets, scores)
    savefig("alignment_wide_roc.svg")

    positive_sites_smooth = []
    positive_sites_dirichlet = []
    negative_sites_smooth = []
    negative_sites_dirichlet = []

    for i in eachindex(result_tuples)

        smooth_negative_result_df, smooth_positive_result_df, dirichlet_negative_result_df, dirichlet_positive_result_df, negative_alphavec, negative_betavec, positive_alphavec, positive_betavec = result_tuples[i]

        for i in eachindex(positive_alphavec)

            if positive_alphavec[i] < positive_betavec[i]
                push!(positive_sites_smooth, smooth_positive_result_df.positive_posterior[i])
                push!(positive_sites_dirichlet, dirichlet_positive_result_df.positive_posterior[i])
            else
                push!(negative_sites_smooth, smooth_positive_result_df.positive_posterior[i])
                push!(negative_sites_dirichlet, dirichlet_positive_result_df.positive_posterior[i])
            end

        end



    end

    rectangle(w, h, x, y) = Shape(x .+ [0, w, w, 0], y .+ [0, 0, h, h])

    rff_outperformance = sum((positive_sites_smooth .> 0.5) .* (positive_sites_dirichlet .< 0.5)) / length(positive_sites_smooth)
    dirichlet_outperformance = sum((positive_sites_smooth .< 0.5) .* (positive_sites_dirichlet .> 0.5)) / length(positive_sites_smooth)
    scatter(positive_sites_smooth, positive_sites_dirichlet, label="Positive site posteriors", size=(1000, 1000), color="red", legend=:outertopright, markerstrokewidth=0, alpha=0.75)

    # plot!(rectangle(0.5,0.5,0.5,0.0), opacity = 0.3, label = "Region of RFF outperformance, share="*string(round(rff_outperformance*100,sigdigits=2))*"%",color="green")
    # plot!(rectangle(0.5,0.5,0,0.5), opacity = 0.3, label = "Region of Dirichlet outperformance, share ="*string(round(dirichlet_outperformance*100,sigdigits=2))*"%",color="red")



    rff_outperformance = sum((negative_sites_smooth .< 0.5) .* (negative_sites_dirichlet .> 0.5)) / length(negative_sites_smooth)
    dirichlet_outperformance = sum((negative_sites_smooth .> 0.5) .* (negative_sites_dirichlet .< 0.5)) / length(negative_sites_smooth)
    scatter!(negative_sites_smooth, negative_sites_dirichlet, label="Negative site posteriors", size=(1000, 1000), color="blue", legend=:outertopright, markerstrokewidth=0, alpha=0.75)
    title!("Scatter plot of Dirichlet prediction vs RFF prediction")
    xlabel!("Smooth prediction")
    ylabel!("Dirichlet prediction")
    # plot!(rectangle(0.5,0.5,0.0,0.5), opacity = 0.3, label = "Region of RFF outperformance, share="*string(round(rff_outperformance*100,sigdigits=2))*"%",color="green")
    # plot!(rectangle(0.5,0.5,0.5,0.0), opacity = 0.3, label = "Region of Dirichlet outperformance, share ="*string(round(dirichlet_outperformance*100,sigdigits=2))*"%",color="red")
    plot!([0, 1], [0, 1], line=(1, :red), label="y=x", markerstrokewidth=0)
    savefig("scatter.svg")

    smooth_scores = [positive_sites_smooth; negative_sites_smooth]
    smooth_targets = Int64.([exp.(zeros(length(positive_sites_smooth))); zeros(length(negative_sites_smooth))])
    dirichlet_scores = [positive_sites_dirichlet; negative_sites_dirichlet]
    dirichlet_targets = Int64.([exp.(zeros(length(positive_sites_dirichlet))); zeros(length(negative_sites_dirichlet))])

    rocplot(smooth_targets, smooth_scores, label="RFF", fillalpha=0)
    rocplot!(dirichlet_targets, dirichlet_scores, label="Dirichlet", fillalpha=0)
    savefig("single_site_roc.svg")
    prplot(smooth_targets, smooth_scores, label="RFF", fillalpha=0)
    prplot!(dirichlet_targets, dirichlet_scores, label="Dirichlet", fillalpha=0)
    savefig("single_site_prc.svg")


    histogram(positive_sites_smooth, color="red", bins=0:0.05:1, alpha=0.6, xlabel="P(β > α | data)", label="Positive class")
    histogram!(negative_sites_smooth, color="blue", bins=0:0.05:1, alpha=0.6, label="Negative class")
    title!("RFF classification histogram")
    savefig("rff_histogram.svg")

    histogram(positive_sites_dirichlet, color="red", bins=0:0.05:1, alpha=0.6, xlabel="P(β > α | data)", label="Positive class")
    histogram!(negative_sites_dirichlet, color="blue", bins=0:0.05:1, alpha=0.6, label="Negative class")
    title!("Dirichlet classification histogram")
    savefig("dirichlet_histogram.svg")


end

function generate_dependent_trees(argument_dictionary)
    simulated_small_trees = []
    simulated_big_trees = []

    positive_sites = parse(Int64, argument_dictionary["positive_sites"])

    ntrees = parse(Int64, argument_dictionary["ntrees"])
    ntaxa = parse(Int64, argument_dictionary["ntaxa"])
    nsites = parse(Int64, argument_dictionary["nsites"])
    nsites_big = parse(Int64, argument_dictionary["nsites_big"])
    K = parse(Int64, argument_dictionary["K"])
    iters = parse(Int64, argument_dictionary["iters"])

    Threads.@threads for i in 1:ntrees

        small_tree = standard_tree_sim(ntaxa)
        big_tree = standard_tree_sim(ntaxa)

        small_alphavec, small_betavec, big_alphavec, big_betavec = simulate_dependent_trees(nsites, nsites_big, 0.2,0.3,1,positive_sites)

        small_nucs, small_seqnams, small_tre = sim_alphabeta_seqs(small_alphavec, small_betavec,
            small_tree, outpath=string(i) * "_small_simdata_FUBAR",
            CodonMolecularEvolution.demo_nucmat, CodonMolecularEvolution.demo_f3x4) #Note the nucmat and f3x4, not from real data

        push!(simulated_small_trees, (small_seqnams, small_nucs, small_tre, small_alphavec, small_betavec))


        big_nucs, big_seqnams, big_tre = sim_alphabeta_seqs(big_alphavec, big_betavec,
            big_tree, outpath=string(i) * "_big_simdata_FUBAR",
            CodonMolecularEvolution.demo_nucmat, CodonMolecularEvolution.demo_f3x4) #Note the nucmat and f3x4, not from real data

        push!(simulated_big_trees, (big_seqnams, big_nucs, big_tre, big_alphavec, big_betavec))

        small_alphabeta_dataframe = DataFrame(alpha = small_alphavec, beta = small_betavec)
        small_alphabeta_dataframe = DataFrame(alpha = big_alphavec, beta = big_betavec)

        CSV.write(string(i)*"_small_alphabeta.csv", small_alphabeta_dataframe)
        CSV.write(string(i)*"_big_alphabeta.csv", small_alphabeta_dataframe)



    end

    return simulated_small_trees, simulated_big_trees

end

function read_dependent_trees(argument_dictionary)

    simulated_small_trees = []
    simulated_big_trees = []

    tree_directory = argument_dictionary["tree_directory"]

    tree_files = readdir(tree_directory)

    tree_count = Int64(length(tree_files) / 6)

    for i in 1:tree_count

        big_seqnams, big_nucs = read_fasta(tree_directory*"/"*string(i)*"_big_simdata_FUBAR.fasta")
        big_tre = read_newick_tree(tree_directory*"/"*string(i)*"_big_simdata_FUBAR.tre")

        small_seqnams, small_nucs = read_fasta(tree_directory*"/"*string(i)*"_small_simdata_FUBAR.fasta")
        small_tre = read_newick_tree(tree_directory*"/"*string(i)*"_small_simdata_FUBAR.tre")

        big_alphabeta_frame = CSV.read(tree_directory*"/"*string(i)*"_big_alphabeta.csv", DataFrame)
        small_alphabeta_frame = CSV.read(tree_directory*"/"*string(i)*"_small_alphabeta.csv", DataFrame)

        big_alphavec, big_betavec = (big_alphabeta_frame.alpha, big_alphabeta_frame.beta)
        small_alphavec, small_betavec = (small_alphabeta_frame.alpha, small_alphabeta_frame.beta)

        push!(simulated_big_trees, (big_seqnams, big_nucs, big_tre, big_alphavec, big_betavec))
        push!(simulated_small_trees, (small_seqnams, small_nucs, small_tre, small_alphavec, small_betavec))

    end

    return simulated_small_trees, simulated_big_trees


end

function run_dependent_tree_simulation(argument_dictionary)

    nsites = parse(Int64, argument_dictionary["nsites"])
    nsites_big = parse(Int64, argument_dictionary["nsites_big"])
    K = parse(Int64, argument_dictionary["K"])
    iters = parse(Int64, argument_dictionary["iters"])

    simulated_small_trees = []
    simulated_big_trees = []
    if argument_dictionary["tree_directory"] == ""
        simulated_small_trees, simulated_big_trees = generate_dependent_trees(argument_dictionary)
    else
        simulated_small_trees, simulated_big_trees = read_dependent_trees(argument_dictionary)
    end
   

    small_alignment_wide = []
    big_alignment_wide = []

    series_annotation = []

    Threads.@threads for i in eachindex(simulated_small_trees)

        small_seqnams, small_nucs, small_tre, small_alphavec, small_betavec = simulated_small_trees[i]
        big_seqnams, big_nucs, big_tre, big_alphavec, big_betavec = simulated_big_trees[i]

        # Perform the init_2_grid calculation only once 

        f_grid_small = alphabetagrid(small_seqnams, small_nucs, newick(small_tre))
        f_grid_big = alphabetagrid(big_seqnams, big_nucs, newick(big_tre))


        dispatch = CodonMolecularEvolution.FUBARweightedpos()

        smooth_small_result_df, smooth_small_tuple = smoothFUBAR(dispatch, f_grid_small, string(i) * "_small_smooth"; K=K, HMC_samples=iters, plots=false)
        smooth_big_result_df, smooth_big_tuple = smoothFUBAR(dispatch, f_grid_big, string(i) * "_big_smooth"; K=K, HMC_samples=iters, plots=false)



        dirichlet_small_result_df, dirichlet_small_θ = FUBAR(f_grid_small, string(i) * "_small_dirichlet", plots=false)
        dirichlet_big_result_df, dirichlet_big_θ = FUBAR(f_grid_big, string(i) * "_big_dirichlet", plots=false)


        push!(small_alignment_wide, smooth_small_tuple.global_posterior_probability)
        push!(big_alignment_wide, smooth_big_tuple.global_posterior_probability)

        write(string(i) * "_small_mixiness.txt", string(smooth_small_tuple.mixing))
        write(string(i) * "_big_mixiness.txt", string(smooth_big_tuple.mixing))

        push!(series_annotation, "ω="*string(round(maximum(small_betavec ./ small_alphavec),digits=2)))


    end

    

    scatter(small_alignment_wide, big_alignment_wide, color="blue", legend=:outertopright, markerstrokewidth=0, alpha=0.75, label = "Predicted probabilites", series_annotation = text.(series_annotation, :left, 5), size = (1000,1000))
    title!("Scatter plot of predictions on alignments of varying sizes")
    xlabel!("Prediction on "*string(nsites)*" sites")
    ylabel!("Prediction on "*string(nsites_big)*" sites")
    plot!([0, 1], [0, 1], line=(1, :red), label="y=x", markerstrokewidth=0)
    savefig("scatter_alignment_wide.svg")

    result_frame = DataFrame(small_alignment_wide = small_alignment_wide, big_alignment_wide = big_alignment_wide)
    CSV.write("alignment_wide_test_small_big.csv",result_frame)

end



function calculate_p(samples, K)

    ψ_samples = [samples[i].z.θ[K+1] for i in 1:length(samples)]

    ψ_positive_given_data = sum(ψ_samples .> 0) / length(ψ_samples)

    return ψ_positive_given_data

end

function calculate_mixiness(samples, K)

    ψ_samples = [samples[i].z.θ[K+1] for i in 1:length(samples)]

    switches = 0

    for i in 1:(length(ψ_samples)-1)

        if ψ_samples[i] * ψ_samples[i+1] < 0
            switches = switches + 0.5
        end

    end

    return switches / length(ψ_samples)

end

if parse(Bool, argument_dictionary["only_simulate"])
    if argument_dictionary["simulation_type"] == "independent"
        generate_independent_trees(argument_dictionary)
    else
        generate_dependent_trees(argument_dictionary)
    end
    exit()
end
if argument_dictionary["simulation_type"] == "independent"
    run_independent_tree_simulation(argument_dictionary)
elseif argument_dictionary["simulation_type"] == "dependent"
    run_dependent_tree_simulation(argument_dictionary)
end