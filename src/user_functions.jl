### Naming convention:
### (1) variable names: snake_case
### (2) structure names: CamelCase
### (3) function names: SCREAMING
### (4) function names user-exposed: snake_case

module user_functions

using Distributed
using ProgressMeter

include("functions.jl")
using .functions: PileupLine, SyncxLine, LocusAlleleCounts, Window
using .functions: PARSE, SPLIT, MERGE, PILEUP2SYNCX, LOAD, FILTER, IMPUTE, SAVE

function pileup2syncx(pileup::String; out::String="")::String
    # using Distributed
    # Distributed.addprocs(length(Sys.cpu_info())-1)
    # @everywhere using ProgressMeter
    # @everywhere include("functions.jl")
    # @everywhere using .functions: PileupLine, LocusAlleleCounts, Window
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
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        digits = length(string(length(positions_init)))
        id = lpad(i, (digits-length(string(i))), "0")
        tmp = string(pileup, "-PILEUP2SYNCX-", id, ".tmp")
        filename = PILEUP2SYNCX(pileup, init, term, tmp)
        [filename]
    end
    ### Sort the output files from parallel processing and merge
    MERGE(filenames_out, out)
    return(out)
end

function filter(pileup::String, outype::String; maximum_missing_fraction::Float64=0.10, alpha1::Float64=0.05, maf::Float64=0.01, alpha2::Float64=0.50, minimum_coverage::Int64=5, out::String="")::String
    # using Distributed
    # Distributed.addprocs(length(Sys.cpu_info())-1)
    # @everywhere using ProgressMeter
    # @everywhere include("functions.jl")
    # @everywhere using .functions: PileupLine, LocusAlleleCounts, Window
    # @everywhere using .functions: SPLIT, MERGE, FILTER
    # pileup = "/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup"    
    # # pileup = "/home/jeffersonfparil/Documents/poolgen/test/test_2.pileup"    
    # outype = ["pileup", "syncx"][1]
    # maximum_missing_fraction = 0.10
    # alpha1 = 0.05
    # maf = 0.001
    # alpha2 = 0.50
    # minimum_coverage = 1
    # out = ""
    ### Define output file if not specified
    if out == ""
        if outype == "pileup"
            out = string(join(split(pileup, ".")[1:(end-1)], "."), "-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, ".pileup")
        elseif outype == "syncx"
            out = string(join(split(pileup, ".")[1:(end-1)], "."), "-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, ".syncx")
        end
    end
    ### Find file positions for parallel processing
    threads = length(Distributed.workers())
    positions_init, positions_term = SPLIT(threads, pileup)
    ### Filter loci by minimum allele frequency (alpha1 and maf) and minimum coverage (alpha2 and minimum_coverage)
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        digits = length(string(length(positions_init)))
        id = lpad(i, (digits-length(string(i))), "0")
        tmp = string(pileup, "-FILTERED-", id, ".tmp")
        filename = FILTER(pileup, outype, init, term, maximum_missing_fraction, alpha1, maf, alpha2, minimum_coverage, tmp)
        [filename]
    end
    ### Sort the output files from parallel processing and merge
    MERGE(filenames_out, out)
    return(out)
end

function filter(syncx::String; maximum_missing_fraction::Float64=0.10, alpha1::Float64=0.05, maf::Float64=0.01, alpha2::Float64=0.50, minimum_coverage::Int64=5, out::String="")::String
    # using Distributed
    # Distributed.addprocs(length(Sys.cpu_info())-1)
    # @everywhere using ProgressMeter
    # @everywhere include("structs.jl")
    # @everywhere include("functions.jl")
    # @everywhere include("user_functions.jl")
    # @everywhere using .structs: PileupLine, LocusAlleleCounts, Window
    # @everywhere using .functions: SPLIT, MERGE, FILTER
    # @everywhere using .user_functions: pileup2syncx
    # # syncx = pileup2syncx("/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup")
    # syncx = pileup2syncx("/home/jeffersonfparil/Documents/poolgen/test/test_2.pileup")
    # maximum_missing_fraction = 0.90
    # alpha1 = 0.05
    # maf = 0.001
    # alpha2 = 0.50
    # minimum_coverage = 1
    # out = ""
    ### Define output file if not specified
    if out == ""
        out = string(join(split(syncx, ".")[1:(end-1)], "."), "-FILTERED-missing_", maximum_missing_fraction, "-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, ".syncx")
    end
    ### Find file positions for parallel processing
    threads = length(Distributed.workers())
    positions_init, positions_term = SPLIT(threads, syncx)
    ### Filter loci by minimum allele frequency (alpha1 and maf) and minimum coverage (alpha2 and minimum_coverage)
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        digits = length(string(length(positions_init)))
        id = lpad(i, (digits-length(string(i))), "0")
        tmp = string(syncx, "-FILTERED-", id, ".tmp")
        filename = FILTER(syncx, init, term, maximum_missing_fraction, alpha1, maf, alpha2, minimum_coverage, tmp)
        [filename]
    end
    ### Sort the output files from parallel processing and merge
    MERGE(filenames_out, out)
    return(out)
end

function impute(filename::String; window_size::Int=100, model::String=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance::Bool=true, distance_nPCA::Int=3, out::String="")::String
    # using Distributed
    # Distributed.addprocs(length(Sys.cpu_info())-1)
    # @everywhere using ProgressMeter
    # @everywhere include("functions.jl")
    # @everywhere include("user_functions.jl")
    # @everywhere using .functions: PileupLine, SyncxLine, LocusAlleleCounts, Window
    # @everywhere using .functions: PARSE, SPLIT, MERGE, PILEUP2SYNCX, IMPUTE
    # @everywhere using .user_functions: pileup2syncx
    # filename = "/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup"
    # # filename = "/home/jeffersonfparil/Documents/poolgen/test/test_2.pileup"
    # # filename = "/home/jeffersonfparil/Documents/poolgen/test/test_1.syncx"
    # # filename = "/home/jeffersonfparil/Documents/poolgen/test/test_2.syncx"
    # window_size = 20
    # # window_size = 100
    # model = "LASSO"
    # distance = true
    # out = ""
    ### Define output file if not specified
    if out == ""
        out = string(join(split(filename, '.')[1:(end-1)], '.'), "-IMPUTED-window_", window_size, "-model_", model, "-distance_", distance, ".syncx")
    end
    ### Define the full path to the input and output files since calling functions within @distributed loop will revert back to the root directory from where julia was executed from
    if dirname(filename) == ""
        filename = string(pwd(), "/", filename)
    end
    if dirname(out) == ""
        out = string(pwd(), "/", out)
    end
    ### Convert to syncx if the input file is in pileup format: output syncx - filename of syncx file
    syncx = try
        file = open(filename, "r")
        PARSE(SyncxLine(1, readline(file)))
        close(file)
        filename
    catch
        pileup2syncx(filename)
    end
    ### Find file positions for parallel processing
    threads = length(Distributed.workers())
    positions_init, positions_term = SPLIT(threads, syncx, window_size)
    idx = collect(1:length(positions_term))[positions_term .== maximum(positions_term)][1]
    positions_init = positions_init[1:idx]
    positions_term = positions_term[1:idx]
    ### Impute
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        digits = length(string(length(positions_init)))
        id = lpad(i, (digits-length(string(i))), "0")
        tmp = string(syncx, "-IMPUTED-", id, ".syncx.tmp")
        filename = IMPUTE(syncx, init, term, window_size, model, distance, distance_nPCA, tmp)
        [filename]
    end
    ### Sort the chunks so we maintain the one-to-one correspondence between input and output loci arrangement
    if length(filenames_out) > 1
        MERGE(filenames_out, window_size, out)
    else
        MERGE(filenames_out, out)
    end
    return(out)
end

end