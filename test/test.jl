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

poolgen.impute("test/test.pileup", window_size=20, model="Mean", distance=false)
poolgen.impute("test/test.pileup", window_size=20, model="OLS", distance=false)
poolgen.impute("test/test.pileup", window_size=20, model="OLS", distance=true)
poolgen.impute("test/test.pileup", window_size=20, model="RR", distance=false)
poolgen.impute("test/test.pileup", window_size=20, model="RR", distance=true)
poolgen.impute("test/test.pileup", window_size=20, model="LASSO", distance=false)
poolgen.impute("test/test.pileup", window_size=20, model="LASSO", distance=true)
poolgen.impute("test/test.pileup", window_size=20, model="GLMNET", distance=false)
poolgen.impute("test/test.pileup", window_size=20, model="GLMNET", distance=true)


# poolgen.simulate
# poolgen.gwalpha
# poolgen.genomic_prediction
# poolgen.genomic_prediction_CV