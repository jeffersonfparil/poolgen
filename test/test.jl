githubci = parse(Bool, ARGS[1])

### Load libraries
using Distributed
Distributed.addprocs(length(Sys.cpu_info())-1)
if githubci
    using Pkg
    Pkg.add(url="https://github.com/jeffersonfparil/poolgen.git")
    @everywhere using poolgen
else
    @everywhere include("/home/jeffersonfparil/Documents/poolgen/src/poolgen.jl")
end

println("##########################################")
println("Pileup to syncx conversion")
println("##########################################")
@time poolgen.pileup2syncx("test/test_1.pileup",
               out="")

@time poolgen.pileup2syncx("test/test_2.pileup",
               out="")

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

mv(string("test/test_1-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, ".syncx"),
   string("test/test_1-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, "-FROM_PILEUP.syncx"))
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

mv(string("test/test_2-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, ".syncx"),
   string("test/test_2-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, "-FROM_PILEUP.syncx"))
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
