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
println("Test GP workflow")
println("##########################################")
# @time poolgen.pileup2syncx("test/test_2.pileup",
@time poolgen.pileup2syncx("test/test_1.pileup",
                           out="test/test_GP_workflow-01_RAW.syncx")
@time poolgen.filter("test/test_GP_workflow-01_RAW.syncx",
                     maximum_missing_fraction=0.50,
                     alpha1=0.01,
                     maf=0.0,
                     alpha2=0.50,
                     minimum_coverage=1,
                     out="test/test_GP_workflow-02_FILTERED.syncx")
@time poolgen.impute("test/test_GP_workflow-02_FILTERED.syncx",
                    window_size=100,
                    model="OLS",
                    distance=true,
                    out="test/test_GP_workflow-03_IMPUTED.syncx")


# 
n = 5                 ### number of founders
m = 10_000            ### number of loci
l = 135_000_000       ### total genome length
k = 5                 ### number of chromosomes
ϵ = Int(1e+15)        ### some arbitrarily large number to signify the distance at which LD is nil
a = 2                 ### number of alleles per locus
vec_chr_lengths = [0] ### chromosome lengths
vec_chr_names = [""]  ### chromosome names 
dist_noLD = 500_000   ### distance at which LD is nil (related to ϵ)
o = 1_000               ### total number of simulated individuals
t = 10                ### number of random mating constant population size generation to simulate
nQTL = 10             ### number of QTL to simulate
heritability = 0.5    ### narrow(broad)-sense heritability as only additive effects are simulated
LD_chr = ""           ### chromosome to calculate LD decay from
LD_n_pairs = 10_000   ### number of randomly sampled pairs of loci to calculate LD
plot_LD = false       ### simulate# plot simulated LD decay
@time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
npools = 50
@time X, y = POOL(_X, _y, npools)
nfold = 10
nrep = 3
@time syncx, phenotype = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X, y)
maf = 0.001
delimiter = ","
header = true
id_col = 1
phenotype_col = 2
missing_strings = ["NA", "NAN", "NaN", "missing", ""]
FE_method = ["CANONICAL", "N<<P"][2]
out = ""
cv_tsv = CV_OLS_MULTIVAR(nfold, nrep, syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, FE_method, out)
###########################
nfold = 10
nrep = 3
syncx = "Simulated-16663168544.syncx"
maf = 0.0001
phenotype = "Simulated-16663168544.csv"
delimiter = ","
header = true
id_col = 1
phenotype_col = 2
missing_strings = ["NA", "NAN", "NaN", "missing", ""]
# FE_method = ["CANONICAL", "N<<P"][2]
# model  = poolgen.user_functions.functions.OLS_MULTIVAR; params = ["N<<P"]
# model  = poolgen.user_functions.functions.ELA_MULTIVAR; params = [1.0]
# #########################
using Optim
_covariate = ["", "XTX", "COR"][2]
_model = ["GBLUP", "RRBLUP"][1]
_method = ["ML", "REML"][1]
_FE_method = ["CANONICAL", "N<<P"][2]
_inner_optimizer=[LBFGS(), BFGS(), SimulatedAnnealing(), GradientDescent(), NelderMead()][1]
_optim_trace = false
model = poolgen.user_functions.functions.LMM_MULTIVAR; params = [_covariate, _model, _method, _FE_method, _inner_optimizer, _optim_trace]
out = ""
# poolgen.user_functions.functions.CV_MULTIVAR(nfold, nrep, syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, FE_method, out)
poolgen.user_functions.functions.CV_MULTIVAR(nfold, nrep, syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, model, params, out)
##########################
