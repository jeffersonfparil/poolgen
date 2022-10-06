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
using .functions: PARSE, SAVE, SPLIT, MERGE, CONVERT, PILEUP2SYNCX, FILTER, IMPUTE, GWALPHA
using .functions: GWAS_PLOT_MANHATTAN, ESTIMATE_LOD
using .functions: SIMULATE, POOL, EXPORT_SIMULATED_DATA

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
###########
### I/O ###
###########
function convert(syncx_or_sync::String; out::String="")::String
    # using Distributed
    # Distributed.addprocs(length(Sys.cpu_info())-1)
    # @everywhere using ProgressMeter
    # @everywhere include("functions.jl")
    # @everywhere using .functions: PileupLine, LocusAlleleCounts, Window
    # @everywhere using .functions: SPLIT, MERGE, CONVERT
    # syncx_or_sync = "/home/jeffersonfparil/Documents/poolgen/test/test_3.syncx"
    # # syncx_or_sync = "/home/jeffersonfparil/Documents/poolgen/test/test_3-A.syncx"
    # out = ""
    ### Define output file if not specified
    if out == ""
        extension = [split(syncx_or_sync, ".")[end][end] == 'x' ? ".sync" : ".syncx"][1]
        out = string(join(split(syncx_or_sync, ".")[1:(end-1)], "."), extension)
    end
    ### Find file positions for parallel processing
    threads = length(Distributed.workers())
    positions_init, positions_term = SPLIT(threads, syncx_or_sync)
    digit = length(string(length(positions_init)))
    ### Convert into syncx or sync format
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        id = lpad(i, (digit+1)-length(string(i)), "0")
        tmp = string(syncx_or_sync, "-CONVERT2SYNCXOR2SYNC-", id, ".tmp")
        filename = CONVERT(syncx_or_sync, init, term, tmp)
        [filename]
    end
    ### Sort the output files from parallel processing and merge
    MERGE(filenames_out, out)
    return(out)
end

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
    digit = length(string(length(positions_init)))
    ### Parse to convert from pileup to syncx format
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        id = lpad(i, (digit+1)-length(string(i)), "0")
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
            out = string(join(split(pileup, ".")[1:(end-1)], "."), "-FILTERED-missing_", maximum_missing_fraction, "-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, ".pileup")
        elseif outype == "syncx"
            out = string(join(split(pileup, ".")[1:(end-1)], "."), "-FILTERED-missing_", maximum_missing_fraction, "-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, ".syncx")
        end
    end
    ### Find file positions for parallel processing
    threads = length(Distributed.workers())
    positions_init, positions_term = SPLIT(threads, pileup)
    digit = length(string(length(positions_init)))
    ### Filter loci by minimum allele frequency (alpha1 and maf) and minimum coverage (alpha2 and minimum_coverage)
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        id = lpad(i, (digit+1)-length(string(i)), "0")
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
    digit = length(string(length(positions_init)))
    ### Filter loci by minimum allele frequency (alpha1 and maf) and minimum coverage (alpha2 and minimum_coverage)
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        id = lpad(i, (digit+1)-length(string(i)), "0")
        tmp = string(syncx, "-FILTERED-", id, ".tmp")
        filename = FILTER(syncx, init, term, maximum_missing_fraction, alpha1, maf, alpha2, minimum_coverage, tmp)
        [filename]
    end
    ### Sort the output files from parallel processing and merge
    MERGE(filenames_out, out)
    return(out)
end

function impute(filename::String; window_size::Int=100, model::String=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance::Bool=true, out::String="")::String
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
    digit = length(string(length(positions_init)))
    ### Impute
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        id = lpad(i, (digit+1)-length(string(i)), "0")
        tmp = string(syncx, "-IMPUTED-", id, ".syncx.tmp")
        filename = IMPUTE(syncx, init, term, window_size, model, distance, tmp)
        [filename]
    end
    ### Merge the chunks
    if length(filenames_out) > 1
        MERGE(filenames_out, window_size, out)
    else
        MERGE(filenames_out, out)
    end
    return(out)
end

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################
### SIMULATION ###
##################
function simulate(;n::Int64, m::Int64, l::Int64, k::Int64, ϵ::Int64=Int(1e+15), a::Int64=2, vec_chr_lengths::Vector{Int64}=Int64.([0]), vec_chr_names::Vector{String}=[""], dist_noLD::Int64=10_000, o::Int64=1_000, t::Int64=10, nQTL::Int64=10, heritability::Float64=0.5, npools::Int64=5, plot_LD::Bool=true)::Tuple{String, String, String, String, String, String}
    # n = 5                 ### number of founders
    # m = 10_000            ### number of loci
    # l = 135_000_000       ### total genome length
    # k = 5                 ### number of chromosomes
    # ϵ = Int(1e+15)        ### some arbitrarily large number to signify the distance at which LD is nil
    # a = 2                 ### number of alleles per locus
    # vec_chr_lengths = [0] ### chromosome lengths
    # vec_chr_names = [""]  ### chromosome names 
    # dist_noLD = 500_000   ### distance at which LD is nil (related to ϵ)
    # o = 100               ### total number of simulated individuals
    # t = 10                ### number of random mating constant population size generation to simulate
    # nQTL = 10             ### number of QTL to simulate
    # heritability = 0.5    ### narrow(broad)-sense heritability as only additive effects are simulated
    # LD_chr = ""           ### chromosome to calculate LD decay from
    # LD_n_pairs = 10_000   ### number of randomly sampled pairs of loci to calculate LD
    # plot_LD = true        ### simulated LD decay
    @time vec_chr, vec_pos, X, y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, plot_LD)
    G, p = POOL(X, y, npools)
    map, bim, ped, fam = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X, y)
    syncx, csv = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, G, p)
    return(map, bim, ped, fam, syncx, csv)
    # ################################################################################
    # ################################################################################
    # ################################################################################
    # ### GWAS test
    # using Statistics, Distributions, Plots
    # QTL_idx = collect(1:length(b))[b .> 0.0]
    # vec_chr_QTL = vec_chr[QTL_idx]
    # vec_pos_QTL = vec_pos[QTL_idx]
    # q = length(QTL_idx)
    # # n, m = size(X)
    # # b_hat = G \ (p .- mean(p))
    # # p1 = Plots.scatter(abs.(b_hat), legend=false, xlab="Position", ylab="|Estimated effect|", title="Pool-GWAS", markerstrokewidth=0.001, markeralpha=0.1)
    # # for i in 1:q
    # #     qtl = QTL_idx[i]
    # #     Plots.plot!(p1, [qtl, qtl], [0, 1], seriestype=:straightline)
    # #     Plots.annotate!(p1, qtl, 0, Plots.text(string("β", i, "\n", round(b[qtl])), 0, 7, :bottom))
    # # end
    # # b_hat2 = X \ (y .- mean(y))
    # # p2 = Plots.scatter(abs.(b_hat2), legend=false, xlab="Position", ylab="|Estimated effect|", title="Indi-GWAS", markerstrokewidth=0.001, markeralpha=0.1)
    # # for i in 1:q
    # #     qtl = QTL_idx[i]
    # #     Plots.plot!(p2, [qtl, qtl], [0, 1], seriestype=:straightline)
    # #     Plots.annotate!(p2, qtl, 0, Plots.text(string("β", i, "\n", round(b[qtl])), 0, 7, :bottom))
    # # end
    # # p3 = Plots.plot(p1, p2, layout=(2, 1)) ### Pool-GWAS performs better than Ind-GWAS!!!!

    # ### Pause julia
    # ### Ctrl+z
    # ### test in plink and gemma
    # # wget 'https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip'
    # # wget 'https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz'
    # # unzip plink_linux_x86_64_20220402.zip; rm plink_linux_x86_64_20220402.zip
    # # gunzip gemma-0.98.5-linux-static-AMD64.gz; mv gemma-0.98.5-linux-static-AMD64 gemma
    # # chmod +x plink
    # # chmod +x gemma
    # # name=$(ls | grep ".fam$" | head -n1 | cut -d '.' -f1)
    # # time ./plink --file "$name" --make-bed --out "$name"
    # # time ./gemma -bfile "$name" -lm 1 -o "$name"

    # ### Resume julia
    # ### fg
    # using CSV, DataFrames
    # cd("output")
    # fname = readdir()[match.(Regex(".assoc.txt"), readdir()) .!= nothing][1]
    # file = open(fname, "r")
    # data = CSV.read(file, DataFrames.DataFrame)
    # close(file)

    # vec_chr_gemma = string.(data.chr)
    # vec_pos_gemma = data.ps
    # vec_lod_gemma = -log10.(data.p_wald)

    # ### NEED TO IMPROVE GWAS_PLOT_MANHATTAN TO MODEL COORDINATES PER CHROMOSOME IN A MORE INTUITIVE AND EASILY ACCESSIBLE WAY!!!

    # p1, chr_names, pos_init, pos_term,  = GWAS_PLOT_MANHATTAN(vec_chr_gemma, vec_pos_gemma, vec_lod_gemma)
    # GWAS_PLOT_SIMULATED_QTL!(p1, vec_chr_QTL, vec_pos_QTL, b, chr_names, pos_init, pos_term)


    # vec_lod_pool_OLS = ESTIMATE_LOD(b_hat)


end

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#############################
### QUANTITATIVE GENETICS ###
#############################
function gwalpha(;syncx::String, py_phenotype::String, maf::Float64=0.001, penalty::Bool=true, out::String="")::String
    # using Distributed
    # Distributed.addprocs(length(Sys.cpu_info())-1)
    # @everywhere using ProgressMeter
    # @everywhere include("functions.jl")
    # @everywhere include("user_functions.jl")
    # @everywhere using .functions: PileupLine, SyncxLine, LocusAlleleCounts, Window
    # @everywhere using .functions: SPLIT, MERGE, GWALPHA
    # syncx = "/home/jeffersonfparil/Documents/poolgen/test/test_3.syncx"
    # py_phenotype = "/home/jeffersonfparil/Documents/poolgen/test/test_3-pheno-any-filename.py"
    # maf = 0.001
    # penalty = false
    # out = ""
    ### Define output file if not specified
    if out == ""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-GWAlpha-maf_", round(maf, digits=5), ".tsv")
    end
    ### Define the full path to the input and output files since calling functions within @distributed loop will revert back to the root directory from where julia was executed from
    if dirname(syncx) == ""
        syncx = string(pwd(), "/", syncx)
    end
    if dirname(out) == ""
        out = string(pwd(), "/", out)
    end
    ### Find file positions for parallel processing
    threads = length(Distributed.workers())
    positions_init, positions_term = SPLIT(threads, syncx)
    digit = length(string(length(positions_init)))
    ### Impute
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        id = lpad(i, (digit+1)-length(string(i)), "0")
        tmp = string(syncx, "-GWAlpha-", id, ".tsv.tmp")
        filename = GWALPHA(syncx, py_phenotype, init, term, maf, penalty, tmp)
        [filename]
    end
    ### Merge the chunks
    MERGE(filenames_out, out)
    return(out)
end

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
###########################
### POPULATION GENETICS ###
###########################



end