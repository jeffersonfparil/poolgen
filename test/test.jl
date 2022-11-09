githubci = parse(Bool, ARGS[1])

### Load libraries
using Distributed
if githubci
    using Pkg
    Pkg.add(url="https://github.com/jeffersonfparil/poolgen.git")
    using poolgen ### Load poolgen first so we can compile now and no precompilation for each process
    Distributed.addprocs(length(Sys.cpu_info())-1)
    @everywhere using poolgen ### Load poolgen for each process
else
    include("/data-weedomics-2/poolgen/src/poolgen.jl")  ### Load poolgen first so we can compile now and no precompilation for each process
    Distributed.addprocs(length(Sys.cpu_info())-1)
    @everywhere include("/data-weedomics-2/poolgen/src/poolgen.jl") ### Load poolgen for each process
end

println("##########################################")
println("Pileup to syncx conversion")
println("##########################################")
@time poolgen.pileup2syncx("test/test_1.pileup")

@time poolgen.pileup2syncx("test/test_2.pileup")

println("##########################################")
println("Filtering")
println("##########################################")
maximum_missing_fraction = 0.50
alpha1 = 0.05
maf = 0.001
alpha2 = 0.50
minimum_coverage = 10
@time poolgen.filter("test/test_1.pileup",
                     "pileup",
                     maximum_missing_fraction=maximum_missing_fraction,
                     alpha1=alpha1,
                     maf=maf,
                     alpha2=alpha2,
                     minimum_coverage=minimum_coverage,
                     out="")

@time poolgen.filter("test/test_1.pileup",
                     "syncx",
                     maximum_missing_fraction=maximum_missing_fraction,
                     alpha1=alpha1,
                     maf=maf,
                     alpha2=alpha2,
                     minimum_coverage=minimum_coverage,
                     out="")

# mv(string("test/test_1-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, ".syncx"),
#    string("test/test_1-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, "-FROM_PILEUP.syncx"))
@time poolgen.filter("test/test_1.syncx",
                     alpha1=alpha1,
                     maf=maf,
                     alpha2=alpha2,
                     minimum_coverage=minimum_coverage,
                     out="")

maximum_missing_fraction = 0.90
alpha1 = 0.01
maf = 0.001
alpha2 = 0.50
minimum_coverage = 1
@time poolgen.filter("test/test_2.pileup",
                     "pileup",
                     maximum_missing_fraction=maximum_missing_fraction,
                     alpha1=alpha1,
                     maf=maf,
                     alpha2=alpha2,
                     minimum_coverage=minimum_coverage,
                     out="")

@time poolgen.filter("test/test_2.pileup",
                     "syncx",
                     maximum_missing_fraction=maximum_missing_fraction,
                     alpha1=alpha1,
                     maf=maf,
                     alpha2=alpha2,
                     minimum_coverage=minimum_coverage,
                     out="")

# mv(string("test/test_2-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, ".syncx"),
#    string("test/test_2-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, "-FROM_PILEUP.syncx"))
@time poolgen.filter("test/test_2.syncx",
                     maximum_missing_fraction=maximum_missing_fraction,
                     alpha1=alpha1,
                     maf=maf,
                     alpha2=alpha2,
                     minimum_coverage=minimum_coverage,
                     out="")

println("##########################################")
println("Imputation")
println("##########################################")
window_size=20
model=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2]
distance=true
@time poolgen.impute("test/test_1.pileup",
                    window_size=window_size,
                    model=model,
                    distance=distance,
                    out="")

@time poolgen.impute("test/test_2.pileup",
                    window_size=window_size,
                    model=model,
                    distance=distance,
                    out="")

println("##########################################")
println("Simulation")
println("##########################################")
n = 5                   ### number of founders
m = 1000                ### number of loci
l = 135_000_000         ### total genome length
k = 5                   ### number of chromosomes
ϵ = Int(1e+15)          ### some arbitrarily large number to signify the distance at which LD is nil
a = 2                   ### number of alleles per locus
vec_chr_lengths = [0]   ### chromosome lengths
vec_chr_names = [""]    ### chromosome names 
dist_noLD = 500_000     ### distance at which LD is nil (related to ϵ)
o = 1_000               ### total number of simulated individuals
t = 10                  ### number of random mating constant population size generation to simulate
nQTL = 10               ### number of QTL to simulate
heritability = 0.5      ### narrow(broad)-sense heritability as only additive effects are simulated
LD_chr = ""             ### chromosome to calculate LD decay from
LD_n_pairs = 10_000     ### number of randomly sampled pairs of loci to calculate LD
plot_LD = false         ### simulate# plot simulated LD decay
npools = 50
@time map, bim, ped, fam, syncx, csv = poolgen.simulate(n=n,
                                                        m=m,
                                                        l=l,
                                                        k=k,
                                                        ϵ=ϵ,
                                                        a=a,
                                                        vec_chr_lengths=vec_chr_lengths,
                                                        vec_chr_names=vec_chr_names,
                                                        dist_noLD=dist_noLD,
                                                        o=o,
                                                        t=t,
                                                        nQTL=nQTL,
                                                        heritability=heritability,
                                                        npools=npools,
                                                        LD_chr=LD_chr,
                                                        LD_n_pairs=LD_n_pairs,
                                                        plot_LD=plot_LD)

# println("##########################################")
# println("GP")
# println("##########################################")
vec_models = ["OLS", "ELASTIC", "LMM"]
vec_MM_models = ["GBLUP", "RRBLUP"]
maf = 0.001
FE_method = ["CANONICAL", "N<<P"][2]
alpha = 1.0
covariate = ["", "XTX", "COR"][2]
MM_method = ["ML", "REML"][1]
inner_optimizer = ["LBFGS", "BFGS", "SimulatedAnnealing", "GradientDescent", "NelderMead"][1]
optim_trace = false
mat_models = hcat(vcat(vec_models, vec_models[end]), vcat(repeat([vec_MM_models[1]], length(vec_models)), vec_MM_models[end]))
# for i in 1:size(mat_models, 1)
#     model = mat_models[i, 1]
#     MM_model = mat_models[i, 2]
#     println(model)
#     model == "LMM" ? println(MM_model) : nothing
#     @time out = poolgen.genomic_prediction(syncx=syncx, 
#                                            phenotype=csv, 
#                                            model=model,
#                                            maf=maf,
#                                            FE_method=FE_method,
#                                            alpha=alpha,
#                                            covariate=covariate,
#                                            MM_model=MM_model,
#                                            MM_method=MM_method,
#                                            inner_optimizer=inner_optimizer,
#                                            optim_trace=optim_trace)
# end

println("##########################################")
println("GP cross-validation")
println("##########################################")
nfold = 2
nrep = 3
save_plots = false
save_predictions = true
save_summary_plot = false
for i in 1:size(mat_models, 1)
    model = mat_models[i, 1]
    MM_model = mat_models[i, 2]
    println(model)
    model == "LMM" ? println(MM_model) : nothing
    @time out = poolgen.genomic_prediction_CV(nfold=nfold,
                                              nrep=nrep,
                                              model=model,
                                              syncx=syncx,
                                              maf=maf,
                                              phenotype=csv,
                                              FE_method=FE_method,
                                              alpha=alpha,
                                              GBLUP_K=covariate,
                                              MM_model=MM_model,
                                              MM_method=MM_method,
                                              inner_optimizer=inner_optimizer,
                                              optim_trace=optim_trace,
                                              save_plots=save_plots,
                                              save_predictions=save_predictions,
                                              save_summary_plot=save_summary_plot)
end

# println("##########################################")
# println("GP cross-validation - WeedOmics data")
# println("##########################################")
# using Distributed
# include("/data-weedomics-2/poolgen/src/poolgen.jl")
# Distributed.addprocs(length(Sys.cpu_info())-1)
# @everywhere include("/data-weedomics-2/poolgen/src/poolgen.jl")
# syncx = "/data-weedomics-2/2.d_40_populations_model_validation/Not_imputed.syncx"
# # # syncx = "/data-weedomics-2/2.d_40_populations_model_validation/Imputed_geno-mean.syncx"
# # phenotype = "/data-weedomics-2/2.d_40_populations_model_validation/Lolium_phenotype_data-READY.csv"
# phenotype = "/data-weedomics-2/2.d_40_populations_model_validation/phenotype_data.csv"
# delimiter=","
# header=true
# id_col=1
# # phenotype_col=13
# # # Y = poolgen.user_functions.functions.LOAD(phenotype, ",", true, 1, [11])
# # # for p in Y.iid
# # #     println(p)
# # # end
# # # idx = (match.(Regex("ACC"), Y.iid) .!= nothing) .& (match.(Regex("Pool"), Y.iid) .== nothing)
# # # y = Y.phe
# # # y[.!idx, 1] .= missing
# # # Φ = poolgen.user_functions.functions.Phenotype(Y.iid, Y.tid, y)
# # # poolgen.user_functions.functions.SAVE(Φ, replace(phenotype, "Lolium_phenotype_data-READY.csv" => "Lolium_glyphosate_resistance_populations_only.csv"), ",", ["id", "Glyphosate_resistance"])
# # # phenotype = "/data-weedomics-2/2.d_40_populations_model_validation/Lolium_glyphosate_resistance_populations_only.csv"
# # # phenotype_col=2

# # nfold = 10
# nrep = 1
# save_plots = false
# save_predictions = true
# save_summary_plot = true
# maf = 1/(42*4)
# # filter_genotype = true
# filter_genotype = false
# # transform_phenotype = true
# transform_phenotype = false
# standardise = false
# FE_method = ["CANONICAL", "N<<P"][2]
# GBLUP_K = ["", "XTX", "COR"][2]
# MM_method = ["ML", "REML"][1]
# inner_optimizer = ["LBFGS", "BFGS", "SimulatedAnnealing", "GradientDescent", "NelderMead"][1]
# optim_trace = false
# # mat_models = hcat(["OLS", "GLMNET", "GLMNET", "MM"],
# #                   ["", "", "", "GBLUP"],
# #                   [0, 0.0, 1.0, 0])
# mat_models = hcat(["OLS", "GLMNET", "GLMNET"],
#                   ["", "", ""],
#                   [0, 0.0, 1.0])
# vec_phenotype_cols = vcat(collect(10:18), collect(20:22), collect(28:33))
# Y = poolgen.user_functions.functions.LOAD(phenotype, ",", true, 1, vec_phenotype_cols)
# for j in 1:length(vec_phenotype_cols)
#     # j = 1
#     pheno = Y.tid[j]
#     phenotype_col = vec_phenotype_cols[j]
#     nfold = sum(.!ismissing.(Y.phe[:,j]))
#     for i in 1:size(mat_models, 1)
#         # i = 1
#         model = mat_models[i, 1]
#         MM_model = mat_models[i, 2]
#         alpha = mat_models[i, 3]
#         if model == "MM"
#             out = string(pheno, "-", model, "_", MM_model, ".tsv")
#         elseif model == "GLMNET"
#             out = string(pheno, "-", model, "-alpha_", round(alpha, digits=4), ".tsv")
#         else
#             out = string(pheno, "-", model, ".tsv")
#         end
#         println("===============================================")
#         println(out)
#         println("===============================================")
#         @time out = poolgen.genomic_prediction_CV(nfold=nfold,
#                                                 nrep=nrep,
#                                                 model=model,
#                                                 syncx=syncx,
#                                                 maf=maf,
#                                                 phenotype=phenotype,
#                                                 delimiter=delimiter,
#                                                 header=header,
#                                                 id_col=id_col,
#                                                 phenotype_col=phenotype_col,
#                                                 filter_genotype=filter_genotype,
#                                                 transform_phenotype=transform_phenotype,
#                                                 standardise=standardise,
#                                                 FE_method=FE_method,
#                                                 alpha=alpha,
#                                                 GBLUP_K=GBLUP_K,
#                                                 MM_model=MM_model,
#                                                 MM_method=MM_method,
#                                                 inner_optimizer=inner_optimizer,
#                                                 optim_trace=optim_trace,
#                                                 save_plots=save_plots,
#                                                 save_predictions=save_predictions,
#                                                 save_summary_plot=save_summary_plot,
#                                                 out=out)
#     end
# end
