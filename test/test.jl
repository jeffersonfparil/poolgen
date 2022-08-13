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
alpha1 = 0.05
maf = 0.001
alpha2 = 0.50
cov = 10
@time poolgen.filter("test/test_1.pileup",
                     "pileup",
                     alpha1=alpha1,
                     maf=maf,
                     alpha2=alpha2,
                     cov=cov,
                     out="")

@time poolgen.filter("test/test_1.pileup",
                     "syncx",
                     alpha1=alpha1,
                     maf=maf,
                     alpha2=alpha2,
                     cov=cov,
                     out="")

mv(string("test/test_1-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, ".syncx"),
   string("test/test_1-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, "-FROM_PILEUP.syncx"))
@time poolgen.filter("test/test_1.syncx",
                     alpha1=alpha1,
                     maf=maf,
                     alpha2=alpha2,
                     cov=cov,
                     out="")

alpha1 = 0.05
maf = 0.0001
alpha2 = 0.50
cov = 1
@time poolgen.filter("test/test_2.pileup",
                     "pileup",
                     alpha1=alpha1,
                     maf=maf,
                     alpha2=alpha2,
                     cov=cov,
                     out="")

@time poolgen.filter("test/test_2.pileup",
                     "syncx",
                     alpha1=alpha1,
                     maf=maf,
                     alpha2=alpha2,
                     cov=cov,
                     out="")

mv(string("test/test_2-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, ".syncx"),
   string("test/test_2-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, "-FROM_PILEUP.syncx"))
@time poolgen.filter("test/test_2.syncx",
                     alpha1=alpha1,
                     maf=maf,
                     alpha2=alpha2,
                     cov=cov,
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

# println("##########################################")
# println("Test GP workflow")
# println("##########################################")

# using Distributed
# Distributed.addprocs(length(Sys.cpu_info())-1)
# @everywhere include("/home/jeffersonfparil/Documents/poolgen/src/poolgen.jl")

# @time poolgen.pileup2syncx("test/test_2.pileup",
#                            out="test/test_GP_workflow-01_RAW.syncx")
# @time poolgen.impute("test/test_GP_workflow-01_RAW.syncx",
#                     window_size=100,
#                     model="OLS",
#                     distance=true,
#                     out="test/test_GP_workflow-02_IMPUTED.syncx")
# @time poolgen.filter("test/test_GP_workflow-02_IMPUTED.syncx",
#                      alpha1=0.01,
#                      maf=0.0001,
#                      alpha2=0.50,
#                      cov=1,
#                      out="test/test_GP_workflow-02_FILTERED.syncx")
