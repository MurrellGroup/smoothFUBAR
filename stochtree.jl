using MolecularEvolution, FASTX, Distributions, EvalMetrics, Plots, DataFrames, CSV

include("CodonMolecularEvolution.jl/src/CodonMolecularEvolution.jl")

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
    while sum(positive_betavec .> positive_alphavec) < pos_sites || length(positive_alphavec) < nsites
        alpha = exp(randn() * σα)
        beta = rand(Exponential(rate))
        if length(positive_alphavec)  < nsites
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

    negative_seqnams, negative_nucs, negative_tre, negative_alphavec, negative_betavec = simulated_negative_trees[i]
    positive_seqnams, positive_nucs, positive_tre, positive_alphavec, positive_betavec = simulated_positive_trees[i]

    # Perform the init_2_grid calculation only once 

    f_grid_negative = CodonMolecularEvolution.alphabetagrid(negative_seqnams, negative_nucs, newick(negative_tre))
    f_grid_positive = CodonMolecularEvolution.alphabetagrid(positive_seqnams, positive_nucs, newick(positive_tre))


    dispatch = CodonMolecularEvolution.FUBARweightedpos()

    smooth_negative_result_df, smooth_negative_θ, negative_samples = CodonMolecularEvolution.smoothFUBAR(dispatch,f_grid_negative,string(i)*"_negative_smooth"; K = K, HMC_samples = iters)
    smooth_positive_result_df, smooth_positive_θ, positive_samples = CodonMolecularEvolution.smoothFUBAR(dispatch, f_grid_positive, string(i)*"_positive_smooth"; K = K, HMC_samples = iters)

    dirichlet_negative_result_df, dirichlet_negative_θ  = CodonMolecularEvolution.FUBAR(f_grid_negative, string(i)*"_negative_dirichlet")
    dirichlet_positive_result_df, dirichlet_positive_θ = CodonMolecularEvolution.FUBAR(f_grid_positive, string(i)*"_positive_dirichlet")
    


    scores[2*i - 1] = calculate_p(negative_samples, K)
    scores[2*i] = calculate_p(positive_samples, K)
    targets[2*i] = 1

    push!(result_tuples, (smooth_negative_result_df, smooth_positive_result_df, dirichlet_negative_result_df, dirichlet_positive_result_df, negative_alphavec, negative_betavec, positive_alphavec, positive_betavec))

    write(string(i)*"_negative_mixiness.txt",string(calculate_mixiness(negative_samples, K)))
    write(string(i)*"_positive_mixiness.txt",string(calculate_mixiness(positive_samples, K)))

end


prplot(targets, scores)
savefig("prplot.svg")
rocplot(targets, scores)
savefig("rocplot.svg")

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