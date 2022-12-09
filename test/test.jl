githubci = parse(Bool, ARGS[1])

### Load libraries
using Distributed
if githubci
    using Pkg
    Pkg.add(url="https://github.com/jeffersonfparil/poolgen.git")
    using poolgen ### Load poolgen first so we can compile now and no precompilation for each process
    Distributed.addprocs(length(Sys.cpu_info()))
    @everywhere using poolgen ### Load poolgen for each process
else
    include("src/poolgen.jl")  ### Load poolgen first so we can compile now and no precompilation for each process
    Distributed.addprocs(length(Sys.cpu_info()))
    @everywhere include("src/poolgen.jl") ### Load poolgen for each process
end

poolgen.convert("test/test.sync")
poolgen.pileup2syncx("test/test.pileup")

poolgen.filter("test/test.pileup", "pileup", maximum_missing_fraction=0.20, alpha1=0.01, maf=0.02, alpha2=0.55, minimum_coverage=7)
poolgen.filter("test/test.pileup", "syncx", maximum_missing_fraction=0.20, alpha1=0.01, maf=0.02, alpha2=0.55, minimum_coverage=7)
poolgen.filter("test/test.syncx", maximum_missing_fraction=0.20, alpha1=0.01, maf=0.02, alpha2=0.55, minimum_coverage=7)

for model in ["Mean", "OLS", "RR", "LASSO", "GLMNET"]
    for distance in [false, true]
        poolgen.impute("test/test.pileup", window_size=20, model=model, distance=distance)
    end
end

poolgen.simulate(n=5, m=10_000, l=135_000_000, k=5, ϵ=Int(1e+15), a=2, vec_chr_lengths=[0], vec_chr_names=[""], dist_noLD=500_000, o=100, t=10, nQTL=10, heritability=0.5, npools=5, LD_chr="", LD_n_pairs=10_000, plot_LD=true, out_geno="test/test-SIMULATED-GENOTYPES", out_pheno="test/test-SIMULATED-PHENOTYPES")

poolgen.gwalpha(syncx="test/test.syncx", py_phenotype="test/test.py", maf=0.001, penalty=false, out="")

for model in ["OLS", "GLMNET", "MM"]
    for GBLUP_K in ["XTX", "COR"]
        for MM_model in ["GBLUP", "RRBLUP"]
            for MM_method in ["ML", "REML"]
                poolgen.genomic_prediction(syncx="test/test.syncx",
                                        phenotype="test/test.csv",
                                        model=model,
                                        filter_genotype=true,
                                        transform_phenotype=true,
                                        standardise=false,
                                        maf=0.01,
                                        delimiter=",",
                                        header=true,
                                        id_col=1,
                                        phenotype_col=2,
                                        missing_strings=["NA", "NAN", "NaN", "missing", ""],
                                        FE_method=["CANONICAL", "N<<P"][2],
                                        alpha=1.0,
                                        GBLUP_K=GBLUP_K,
                                        MM_model=MM_model,
                                        MM_method=MM_method,
                                        MM_inner_optimizer=["GradientDescent", "LBFGS", "BFGS", "SimulatedAnnealing","NelderMead"][1],
                                        MM_optim_trace=true,
                                        out="")
            end
        end
    end
end

for model in ["OLS", "GLMNET", "MM"]
    for GBLUP_K in ["XTX", "COR"]
        for MM_model in ["GBLUP", "RRBLUP"]
            for MM_method in ["ML", "REML"]
                poolgen.genomic_prediction_CV(nrep=10,
                                              nfold=10,
                                              syncx="test/test_Lr.syncx",
                                              phenotype="test/test_Lr.csv",
                                              model=model,
                                              maf=0.001,
                                              delimiter=",",
                                              header=true,
                                              id_col=1,
                                              phenotype_col=2,
                                              missing_strings=["NA", "NAN", "NaN", "missing", ""],
                                              filter_genotype=true,
                                              transform_phenotype=true,
                                              standardise=false,
                                              FE_method=["CANONICAL", "N<<P"][2],
                                              alpha=1.0,
                                              MM_model=MM_model,
                                              MM_method=MM_method,
                                              MM_inner_optimizer=["GradientDescent", "LBFGS", "BFGS", "SimulatedAnnealing", "NelderMead"][1],
                                              MM_optim_trace=false,
                                              GBLUP_K=GBLUP_K,
                                              save_plots=false,
                                              save_predictions=false,
                                              save_summary_plot=false,
                                              out="")
            end
        end
    end
end
