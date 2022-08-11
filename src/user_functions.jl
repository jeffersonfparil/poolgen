### Naming convention:
### (1) variable names: snake_case
### (2) structure names: CamelCase
### (3) function names: SCREAMING
### (4) function names user-exposed: snake_case

module user_functions

using Distributed
using ProgressMeter

include("structs.jl")
include("functions.jl")
using .structs: PileupLine, LocusAlleleCounts, Window
using .functions: SPLIT, MERGE, PILEUP2SYNCX, FILTER, IMPUTE, SAVE

function pileup2syncx(pileup::String; out::String="")
    # using Distributed
    # Distributed.addprocs(length(Sys.cpu_info())-1)
    # @everywhere using ProgressMeter
    # @everywhere include("structs.jl")
    # @everywhere include("functions.jl")
    # @everywhere using .structs: PileupLine, LocusAlleleCounts, Window
    # @everywhere using .functions: SPLIT, MERGE, PILEUP2SYNCX
    # pileup = "/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup"    
    # # pileup = "/home/jeffersonfparil/Documents/poolgen/test/test_2.pileup"    
    # out = ""
    ### Define output file if not specified
    if out == ""
        out = string(join(split(pileup, ".")[1:(end-1)], "."), ".syncx")
    end
    ### Find file positions for parallel processing
    threads = length(Distributed.workers())
    positions_init, positions_term = SPLIT(threads, pileup)
    ### Parse to convert from pileup to syncx format
    @time filenames_out = @sync @showprogress @distributed (append!) for i in 1:length(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        digits = length(string(length(positions_init)))
        id = lpad(i, (digits-length(string(i))), "0")
        tmp = string(pileup, "-FILTERED-", id, ".tmp")
        filename = PILEUP2SYNCX(pileup, init, term, tmp)
        [filename]
    end
    ### Sort the output files from parallel processing and merge
    MERGE(filenames_out, out)
    return(out)
end

function filter(pileup::String; alpha1::Float64=0.05, maf::Float64=0.01, alpha2::Float64=0.50, cov::Int64=5, outype::String=["pileup", "syncx"][1], out::String="")::String
    # using Distributed
    # Distributed.addprocs(length(Sys.cpu_info())-1)
    # @everywhere using ProgressMeter
    # @everywhere include("structs.jl")
    # @everywhere include("functions.jl")
    # @everywhere using .structs: PileupLine, LocusAlleleCounts, Window
    # @everywhere using .functions: SPLIT, FILTER
    # pileup = "/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup"    
    # # pileup = "/home/jeffersonfparil/Documents/poolgen/test/test_2.pileup"    
    # alpha1 = 0.05
    # maf = 0.001
    # alpha2 = 0.50
    # cov = 1
    # outype = ["pileup", "syncx"][1]
    # out = ""
    ### Define output file if not specified
    if out == ""
        if outype == "pileup"
            out = string(join(split(pileup, ".")[1:(end-1)], "."), "-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, ".pileup")
        elseif outype == "syncx"
            out = string(join(split(pileup, ".")[1:(end-1)], "."), "-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", cov, ".syncx")
        end
    end
    ### Find file positions for parallel processing
    threads = length(Distributed.workers())
    positions_init, positions_term = SPLIT(threads, pileup)
    ### Filter loci by minimum allele frequency (alpha1 and maf) and minimum coverage (alpha2 and cov)
    @time filenames_out = @sync @showprogress @distributed (append!) for i in 1:length(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        digits = length(string(length(positions_init)))
        id = lpad(i, (digits-length(string(i))), "0")
        tmp = string(pileup, "-FILTERED-", id, ".tmp")
        filename = FILTER(pileup, init, term, alpha1, maf, alpha2, cov, outype, tmp)
        [filename]
    end
    ### Sort the output files from parallel processing and merge
    MERGE(filenames_out, out)
    return(out)
end

function impute(pileup_with_missing::String; window_size::Int=100, model::String=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance::Bool=true, syncx_imputed::String="")::String
    # using Distributed
    # Distributed.addprocs(length(Sys.cpu_info())-1)
    # @everywhere using ProgressMeter
    # @everywhere include("structs.jl")
    # @everywhere include("functions.jl")
    # @everywhere using .functions: SPLIT, MERGE, IMPUTE
    # pileup_with_missing = "/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup"
    # # pileup_with_missing = "/home/jeffersonfparil/Documents/poolgen/test/test_2.pileup"
    # window_size = 20
    # model = "LASSO"
    # distance = true
    # syncx_imputed = "/home/jeffersonfparil/Documents/poolgen/test/test-IMPUTED.syncx"
    ### Define output file if not specified
    if syncx_imputed == ""
        syncx_imputed = string(join(split(pileup_with_missing, '.')[1:(end-1)], '.'), "-IMPUTED-window_", window_size, "-model_", model, "-distance_", distance, ".syncx")
    end
    ### Define the full path to the input and output files since calling functions within @distributed loop will revert back to the root directory from where julia was executed from
    if dirname(pileup_with_missing) == ""
        pileup_with_missing = string(pwd(), "/", pileup_with_missing)
    end
    if dirname(syncx_imputed) == ""
        syncx_imputed = string(pwd(), "/", syncx_imputed)
    end
    ### Find file positions for parallel processing
    threads = length(Distributed.workers())
    positions_init, positions_term = SPLIT(threads, pileup_with_missing, window_size)
    ### Impute
    @time filenames_out = @sync @showprogress @distributed (append!) for i in 1:length(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        digits = length(string(length(positions_init)))
        id = lpad(i, (digits-length(string(i))), "0")
        tmp = string(pileup_with_missing, "-IMPUTED-", id, ".syncx.tmp")
        filename = IMPUTE(pileup_with_missing, init, term, window_size, model, distance, tmp)
        [filename]
    end
    ### Sort the chunks so we maintain the one-to-one correspondence between input and output loci arrangement
    MERGE(filenames_out, window_size, syncx_imputed)
    return(syncx_imputed)
end

end