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
               alpha1=alpha1,
               maf=maf,
               alpha2=alpha2,
               cov=cov,
               outype="pileup",
               out="")

@time poolgen.filter("test/test_2.pileup",
               alpha1=alpha1,
               maf=maf,
               alpha2=alpha2,
               cov=cov,
               outype="pileup",
               out="")

println("##########################################")
println("Imputation")
println("##########################################")
window_size=20
model=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2]
distance=true
syncx_imputed=""
@time poolgen.impute("test/test_1.pileup",
                    window_size=window_size,
                    model=model,
                    distance=distance,
                    syncx_imputed=syncx_imputed)

@time poolgen.impute("test/test_2.pileup",
                    window_size=window_size,
                    model=model,
                    distance=distance,
                    syncx_imputed=syncx_imputed)
