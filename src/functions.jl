### Naming convention:
### (1) variable names: snake_case
### (2) structure names: CamelCase
### (3) function names: SCREAMING
### (4) function names user-exposed: snake_case

module functions

using Distributed
using Random
using GLMNet
using MultivariateStats
using ProgressMeter
using LinearAlgebra; BLAS.set_num_threads(1) ### Set the number of threads used for linear algebra to 1 to allow for parallel execution functions with matrix computations
using Distributions
using Optim
using StatsBase
using Dates
using Plots

include("structs.jl")
using .structs: PileupLine, SyncxLine, LocusAlleleCounts, Window, PhenotypeLine, Phenotype, MinimisationError
include("functions_io.jl")
using ProgressMeter, Distributions
include("functions_filterTransform.jl")
using Plots
include("functions_simulate.jl")
include("functions_impute.jl")
include("functions_linearModel.jl")
include("functions_gp.jl")
include("functions_gwas.jl")

end
