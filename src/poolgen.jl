module poolgen

include("user_functions.jl")
using .user_functions: convert
using .user_functions: pileup2syncx
using .user_functions: filter
using .user_functions: impute
using .user_functions: simulate
using .user_functions: gwalpha
# using .user_functions: genomic_prediction
using .user_functions: genomic_prediction_CV

end
