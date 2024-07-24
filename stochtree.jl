using CodonMolecularEvolution, MolecularEvolution, FASTX, Distributions, EvalMetrics, Plots, DataFrames, CSV

# include("CodonMolecularEvolution.jl/src/CodonMolecularEvolution.jl")

#=
#Generate a tree, stochastically
negative_stochastic_tree = CodonMolecularEvolution.standard_tree_sim(500)
#Generate alpha,beta values:
negative_alphavec, negative_betavec = exp.(randn(1000) .* 0.1),rand(Beta(0.5,0.5),1000)
#Simulate the data
negative_nucs, negative_seqnams, negative_tre = CodonMolecularEvolution.sim_alphabeta_seqs(negative_alphavec, negative_betavec, 
negative_stochastic_tree, outpath = "negative_simdata_FUBAR",
    CodonMolecularEvolution.demo_nucmat, CodonMolecularEvolution.demo_f3x4) #Note the nucmat and f3x4, not from real data

#Generate a tree, stochastically
positive_stochastic_tree = CodonMolecularEvolution.standard_tree_sim(500)
#Generate alpha,beta values:
positive_alphavec, positive_betavec = exp.(randn(1000) .* 0.1),exp.(randn(1000) .* 0.1)
#Simulate the data
positive_nucs, positive_seqnams, positive_tre = CodonMolecularEvolution.sim_alphabeta_seqs(positive_alphavec, positive_betavec, 
positive_stochastic_tree, outpath = "positive_simdata_FUBAR",
    CodonMolecularEvolution.demo_nucmat, CodonMolecularEvolution.demo_f3x4) #Note the nucmat and f3x4, not from real data


#Run FUBAR on the data returned from the sim:
negative_result_df, negative_for_retabulation, negative_samples = CodonMolecularEvolution.restricted_smoothFUBAR(negative_seqnams, negative_nucs, newick(negative_tre), "negative_simtest_FUBAR");
positive_result_df, positive_for_retabulation, positive_samples = CodonMolecularEvolution.restricted_smoothFUBAR(positive_seqnams, positive_nucs, newick(positive_tre), "positive_simtest_FUBAR");
=#

argument_dictionary = Dict(
  "directory" => ".",
  "positive_sites" => "1",
  "ntrees" => "30",
  "ntaxa" => "50",
  "nsites" => "500",
  "K" => "150",
  "iters" => "2500",
  "tree_directory" => "",
  "only_simulate" => "false")

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

simulated_negative_trees = []
simulated_positive_trees = []

positive_sites = parse(Int64,argument_dictionary["positive_sites"])

ntrees = parse(Int64,argument_dictionary["ntrees"])
ntaxa = parse(Int64, argument_dictionary["ntaxa"])
nsites = parse(Int64, argument_dictionary["nsites"])
K = parse(Int64, argument_dictionary["K"])
iters = parse(Int64, argument_dictionary["iters"])

main_directory = argument_dictionary["directory"]




function simulate_single_positive_site(nsites, σα, rate, multiplier,pos_sites)

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
    while (sum(positive_betavec .> positive_alphavec) < pos_sites) | length(positive_alphavec) < nsites
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


if argument_dictionary["tree_directory"] == ""

    Threads.@threads for i in 1:ntrees

        local_negative_tree = CodonMolecularEvolution.standard_tree_sim(ntaxa)
   
        negative_alphavec, negative_betavec, positive_alphavec, positive_betavec = simulate_single_positive_site(nsites, 0.2, 0.3, 2,positive_sites)

    

        negative_nucs, negative_seqnams, negative_tre = CodonMolecularEvolution.sim_alphabeta_seqs(negative_alphavec, negative_betavec, 
            local_negative_tree, outpath = string(i)*"_negative_simdata_FUBAR",
        CodonMolecularEvolution.demo_nucmat, CodonMolecularEvolution.demo_f3x4) #Note the nucmat and f3x4, not from real data

        push!(simulated_negative_trees, ( negative_seqnams,negative_nucs, negative_tre, negative_alphavec, negative_betavec))


        local_positive_tree = CodonMolecularEvolution.standard_tree_sim(ntaxa)
    
        single_strong_index = rand(1:length(positive_alphavec)) 

        positive_betavec[single_strong_index] = 2 * positive_alphavec[single_strong_index]

        positive_nucs, positive_seqnams, positive_tre = CodonMolecularEvolution.sim_alphabeta_seqs(positive_alphavec, positive_betavec, 
            local_positive_tree, outpath = string(i)*"_positive_simdata_FUBAR",
        CodonMolecularEvolution.demo_nucmat, CodonMolecularEvolution.demo_f3x4) #Note the nucmat and f3x4, not from real data

        push!(simulated_positive_trees, (positive_seqnams, positive_nucs, positive_tre, positive_alphavec, positive_betavec))

        negative_dataframe = DataFrame([negative_alphavec negative_betavec],["alphavec","betavec"])
        positive_dataframe = DataFrame([positive_alphavec positive_betavec],["alphavec","betavec"])

        CSV.write(string(i)*"_negative_vector_simdata.csv",negative_dataframe)
        CSV.write(string(i)*"_positive_vector_simdata.csv",positive_dataframe)

    end

else

    trees = readdir(argument_dictionary["tree_directory"])

    simulation_count = Int64(sum(occursin.(".tre",trees)) / 2)

    

    for i in 1:simulation_count 

        main_directory = argument_dictionary["tree_directory"]*"/"*string(i)

        positive_seqnams, positive_nucs = read_fasta(main_directory*"_positive_simdata_FUBAR.fasta")
        positive_tre = read_newick_tree(main_directory*"_positive_simdata_FUBAR.tre")
        negative_seqnams, negative_nucs = read_fasta(main_directory*"_negative_simdata_FUBAR.fasta")
        negative_tre = read_newick_tree(main_directory*"_negative_simdata_FUBAR.tre")

        positive_dataframe = CSV.read(main_directory*"_positive_vector_simdata.csv", DataFrame)
        negative_dataframe = CSV.read(main_directory*"_negative_vector_simdata.csv", DataFrame)

        negative_alphavec, negative_betavec = (negative_dataframe.alphavec, negative_dataframe.betavec)
        positive_alphavec, positive_betavec = (positive_dataframe.alphavec, positive_dataframe.betavec)

        push!(simulated_positive_trees, (positive_seqnams, positive_nucs, positive_tre, positive_alphavec, positive_betavec))
        push!(simulated_negative_trees, ( negative_seqnams,negative_nucs, negative_tre, negative_alphavec, negative_betavec))

        
    end

end


if parse(Bool, argument_dictionary["only_simulate"])

    exit()

end




function calculate_p(samples, K)

    ψ_samples = [samples[i].z.θ[K + 1] for i in 1:length(samples)]

    ψ_positive_given_data = sum(ψ_samples .> 0) / length(ψ_samples)

    return ψ_positive_given_data

end

function calculate_mixiness(samples, K)

    ψ_samples = [samples[i].z.θ[K + 1] for i in 1:length(samples)]

    switches = 0

    for i in 1:(length(ψ_samples)-1)

        if ψ_samples[i] * ψ_samples[i+1] < 0
            switches = switches + 0.5
        end

    end

    return switches / length(ψ_samples)

end

targets = zeros(2 * ntrees)
scores = zeros(2 * ntrees)

println("Running with: "*string(Threads.nthreads())*" threads.")
println(length(simulated_negative_trees))

result_tuples = []

plotting_tuples = []
dirichlet_plotting_tuples = []
Threads.@threads for i in eachindex(simulated_negative_trees)
    try 
        negative_seqnams, negative_nucs, negative_tre, negative_alphavec, negative_betavec = simulated_negative_trees[i]
        positive_seqnams, positive_nucs, positive_tre, positive_alphavec, positive_betavec = simulated_positive_trees[i]

    # Perform the init_2_grid calculation only once 
    
        f_grid_negative = CodonMolecularEvolution.pre_init2grid(negative_seqnams, negative_nucs, newick(negative_tre),string(i)*"_negative")
        f_grid_positive = CodonMolecularEvolution.pre_init2grid(positive_seqnams, positive_nucs, newick(positive_tre),string(i)*"_positive")


        dispatch = CodonMolecularEvolution.FUBARweightedpos()

        smooth_negative_result_df, smooth_negative_tuple, smooth_negative_samples, smooth_negative_ℓπ = CodonMolecularEvolution.smoothFUBAR_precomputed_f(dispatch,f_grid_negative,string(i)*"_negative_smooth"; K = K, HMC_samples = iters)
        smooth_positive_result_df, smooth_positive_tuple, smooth_positive_samples, smooth_positive_ℓπ = CodonMolecularEvolution.smoothFUBAR_precomputed_f(dispatch, f_grid_positive, string(i)*"_positive_smooth"; K = K, HMC_samples = iters)

        dirichlet_negative_result_df, dirichlet_negative_tuple = CodonMolecularEvolution.FUBAR_precomputed_f(f_grid_negative, string(i)*"_negative_dirichlet")
        dirichlet_positive_result_df, dirichlet_positive_tuple = CodonMolecularEvolution.FUBAR_precomputed_f(f_grid_positive, string(i)*"_positive_dirichlet")

        scores[2*i - 1] = calculate_p(smooth_negative_samples, K)
        scores[2*i] = calculate_p(smooth_positive_samples, K)
        targets[2*i] = 1

        push!(plotting_tuples, (smooth_negative_ℓπ, smooth_negative_samples, smooth_negative_tuple[1], f_grid_negative, string(i)*"_negative_smooth", 200)) # Last is an ugly constant for n_adapts
        push!(plotting_tuples, (smooth_positive_ℓπ, smooth_positive_samples, smooth_positive_tuple[1], f_grid_positive, string(i)*"_positive_smooth", 200))

        push!(dirichlet_plotting_tuples, (dirichlet_negative_tuple[1],f_grid_negative, string(i)*"_negative_dirichlet"))
        push!(dirichlet_plotting_tuples, (dirichlet_positive_tuple[1],f_grid_positive, string(i)*"_positive_dirichlet"))

        push!(result_tuples, (smooth_negative_result_df, smooth_positive_result_df, dirichlet_negative_result_df, dirichlet_positive_result_df, negative_alphavec, negative_betavec, positive_alphavec, positive_betavec))

        write(string(i)*"_negative_mixiness.txt",string(calculate_mixiness(smooth_negative_samples, K)))
        write(string(i)*"_positive_mixiness.txt",string(calculate_mixiness(smooth_positive_samples, K)))
    catch e

    end
end

#=
for i in eachindex(plotting_tuples)
    ℓπ, samples, θ, f_grid, analysis_name, n_adapts  = plotting_tuples[i]

    CodonMolecularEvolution.core_plots(ℓπ, samples, θ, f_grid, analysis_name, n_adapts)
    CodonMolecularEvolution.other_plots(ℓπ, samples, f_grid, analysis_name, n_adapts)
end

for i in eachindex(dirichlet_plotting_tuples)
    θ, f, analysis_name = dirichlet_plotting_tuples[i]

    CodonMolecularEvolution.FUBAR_plot_from_θ(θ, f, analysis_name)
end
=#
#= 
Threads.@threads for i in 1:length(simulated_negative_trees)

    negative_seqnams, negative_nucs, negative_tre, negative_alphavec, negative_betavec = simulated_negative_trees[i]
    positive_seqnams, positive_nucs, positive_tre, positive_alphavec, positive_betavec = simulated_positive_trees[i]

    smooth_negative_result_df, negative_for_retabulation, negative_samples,negative_ℓπ = CodonMolecularEvolution.test_positive_selection_smoothFUBAR(negative_seqnams, negative_nucs, newick(negative_tre), string(i)*"_negative_simtest_FUBAR", HMC_samples = 500);
    
    dirichlet_negative_result_df = CodonMolecularEvolution.FUBAR(negative_seqnams, negative_nucs, newick(negative_tre),string(i)*"_negative_simtest_dirichlet")[1]

    
    scores[2*i - 1] = calculate_p(negative_samples, 50)



    smooth_positive_result_df, positive_for_retabulation, positive_samples, positive_ℓπ = CodonMolecularEvolution.test_positive_selection_smoothFUBAR(positive_seqnams, positive_nucs, newick(positive_tre), string(i)*"_positive_simtest_FUBAR", HMC_samples = 500);

    dirichlet_positive_result_df = CodonMolecularEvolution.FUBAR(positive_seqnams, positive_nucs, newick(positive_tre),string(i)*"_positive_simtest_dirichlet")[1]

    push!(result_tuples, (smooth_negative_result_df, smooth_positive_result_df, 
    dirichlet_negative_result_df, dirichlet_positive_result_df, negative_alphavec, negative_betavec, positive_alphavec, positive_betavec))

    scores[2*i] = calculate_p(positive_samples, 50)
    targets[2*i] = 1
    push!(plotting_tuples, (negative_for_retabulation[2], negative_samples, negative_seqnams, negative_nucs, negative_tre, negative_ℓπ ,string(i)*"_negative_simtest_FUBAR"))
    push!(plotting_tuples, (positive_for_retabulation[2], positive_samples, positive_seqnams, positive_nucs, positive_tre, positive_ℓπ ,string(i)*"_positive_simtest_FUBAR"))

end
=#
# STARTED 20 TREES at 12:45 6 cores

# prplot(targets, scores)
# savefig("prplot.svg")
# rocplot(targets, scores)
# savefig("rocplot.svg")

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

rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

rff_outperformance = sum((positive_sites_smooth .> 0.5) .* (positive_sites_dirichlet .< 0.5)) / length(positive_sites_smooth)
dirichlet_outperformance = sum((positive_sites_smooth .< 0.5) .* (positive_sites_dirichlet .> 0.5)) / length(positive_sites_smooth)
scatter(positive_sites_smooth,positive_sites_dirichlet,label="Positive site posteriors",size=(1000,1000),color="red",legend = :outertopright, markerstrokewidth = 0, alpha = 0.75)

# plot!(rectangle(0.5,0.5,0.5,0.0), opacity = 0.3, label = "Region of RFF outperformance, share="*string(round(rff_outperformance*100,sigdigits=2))*"%",color="green")
# plot!(rectangle(0.5,0.5,0,0.5), opacity = 0.3, label = "Region of Dirichlet outperformance, share ="*string(round(dirichlet_outperformance*100,sigdigits=2))*"%",color="red")



rff_outperformance = sum((negative_sites_smooth .< 0.5) .* (negative_sites_dirichlet .> 0.5)) / length(negative_sites_smooth)
dirichlet_outperformance = sum((negative_sites_smooth .> 0.5) .* (negative_sites_dirichlet .< 0.5)) / length(negative_sites_smooth)
scatter!(negative_sites_smooth,negative_sites_dirichlet,label="Negative site posteriors",size=(1000,1000),color="blue", legend = :outertopright, markerstrokewidth = 0, alpha = 0.75)
title!("Scatter plot of Dirichlet prediction vs RFF prediction")
xlabel!("Smooth prediction")
ylabel!("Dirichlet prediction")
# plot!(rectangle(0.5,0.5,0.0,0.5), opacity = 0.3, label = "Region of RFF outperformance, share="*string(round(rff_outperformance*100,sigdigits=2))*"%",color="green")
# plot!(rectangle(0.5,0.5,0.5,0.0), opacity = 0.3, label = "Region of Dirichlet outperformance, share ="*string(round(dirichlet_outperformance*100,sigdigits=2))*"%",color="red")
plot!([0,1], [0,1], line=(1, :red), label = "y=x",markerstrokewidth=0)
savefig("scatter.svg")

smooth_scores = [positive_sites_smooth; negative_sites_smooth]
smooth_targets = Int64.([exp.(zeros(length(positive_sites_smooth))); zeros(length(negative_sites_smooth))])
dirichlet_scores = [positive_sites_dirichlet; negative_sites_dirichlet]
dirichlet_targets = Int64.([exp.(zeros(length(positive_sites_dirichlet))); zeros(length(negative_sites_dirichlet))])

rocplot(smooth_targets, smooth_scores, label = "RFF", fillalpha = 0)
rocplot!(dirichlet_targets, dirichlet_scores, label = "Dirichlet", fillalpha = 0)
savefig("roc.svg")
prplot(smooth_targets, smooth_scores, label = "RFF", fillalpha = 0)
prplot!(dirichlet_targets, dirichlet_scores, label = "Dirichlet", fillalpha = 0)
savefig("prc.svg")


histogram(positive_sites_smooth, color = "red", bins = 0:0.05:1, alpha = 0.6, xlabel = "P(β > α | data)", label = "Positive class")
histogram!(negative_sites_smooth, color = "blue", bins = 0:0.05:1, alpha = 0.6, label = "Negative class")
title!("RFF classification histogram")
savefig("rff_histogram.svg")

histogram(positive_sites_dirichlet, color = "red", bins = 0:0.05:1, alpha = 0.6, xlabel = "P(β > α | data)", label = "Positive class")
histogram!(negative_sites_dirichlet, color = "blue", bins = 0:0.05:1, alpha = 0.6, label = "Negative class")
title!("Dirichlet classification histogram")
savefig("dirichlet_histogram.svg")

#=
function compare_result_dfs(df_a, df_b, beta_vec, alpha_vec)
    
    log_prediction_value_a = 0
    log_prediction_value_b = 0

    for i in 1:length(df_a.positive_posterior)

        if beta_vec[i] > alpha_vec[i]
            log_prediction_value_a += log(df_a.positive_posterior[i])
            log_prediction_value_b += log(df_b.positive_posterior[i])
        else
            log_prediction_value_a += log(1 - df_a.positive_posterior[i])
            log_prediction_value_b += log(1 - df_b.positive_posterior[i])
        end
    end

    return log_prediction_value_a, log_prediction_value_b

end

negative_smooth_pvs = []
negative_dirichlet_pvs = []
positive_smooth_pvs = []
positive_dirichlet_pvs = []

for i in eachindex(result_tuples)

    smooth_negative_result_df, smooth_positive_result_df, dirichlet_negative_result_df, dirichlet_positive_result_df, negative_alphavec, negative_betavec, positive_alphavec, positive_betavec = result_tuples[i]

    negative_smooth_pv, negative_dirichlet_pv = compare_result_dfs(smooth_negative_result_df, dirichlet_negative_result_df, negative_alphavec, negative_betavec)
    positive_smooth_pv, positive_dirichlet_pv = compare_result_dfs(smooth_positive_result_df, dirichlet_positive_result_df, positive_alphavec, positive_betavec)

    push!(negative_smooth_pvs,negative_smooth_pv)
    push!(negative_dirichlet_pvs,negative_dirichlet_pv)
    push!(positive_smooth_pvs, positive_smooth_pv)
    push!(positive_dirichlet_pvs, positive_dirichlet_pv)


end

prediction_value_df = DataFrame(negative_smooth = negative_smooth_pvs, negative_dirichlet = negative_dirichlet_pvs, positive_smooth = positive_smooth_pvs, positive_dirichlet = positive_dirichlet_pvs)

CSV.write("prediction_values.csv",prediction_value_df)
=#
