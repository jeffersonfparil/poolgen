githubci = parse(Bool, ARGS[1])

### Load libraries
using Distributed
Distributed.addprocs(length(Sys.cpu_info())-1)
if githubci
    using Pkg
    Pkg.add(url="https://github.com/jeffersonfparil/poolgen.git")
    @everywhere using poolgen
else
    @everywhere include("/data-weedomics-2/poolgen/src/poolgen.jl")
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
println("GP")
println("##########################################")
nfold = 2
nrep = 3
model = ["OLS", "ELASTIC-NET", "LMM"][1]
syncx = "test/Simulated-16663168544.syncx"
maf = 0.001
phenotype = "test/Simulated-16663168544.csv"
delimiter = ","
header = true
id_col = 1
phenotype_col = 2
missing_strings = ["NA", "NAN", "NaN", "missing", ""]
FE_method = ["CANONICAL", "N<<P"][2]
alpha = 1.0
covariate = ["", "XTX", "COR"][2]
MM_model = ["GBLUP", "RRBLUP"][1]
MM_method = ["ML", "REML"][1]
inner_optimizer=["LBFGS", "BFGS", "SimulatedAnnealing", "GradientDescent", "NelderMead"][1]
optim_trace = false
out = ""
@time poolgen.genomic_prediction_CV(nfold=nfold, nrep=nrep, model=model, syncx=syncx, maf=maf, phenotype=phenotype, delimiter=delimiter, header=header, id_col=id_col, phenotype_col=phenotype_col, missing_strings=missing_strings, FE_method=FE_method, alpha=alpha, covariate=covariate, MM_model=MM_model, MM_method=MM_method, inner_optimizer=inner_optimizer, optim_trace=optim_trace, out=out)
