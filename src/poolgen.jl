module poolgen

include("functions.jl")
using .functions: convert
using .functions: pileup2syncx
using .functions: filter
using .functions: impute
using .functions: simulate
using .functions: gwalpha
using .functions: genomic_prediction
using .functions: genomic_prediction_CV

end
