### DATA FILTERING AND TRANSFORMATION FUNCTIONS

####### TEST ########
# include("structs.jl")
# using .structs: PileupLine, SyncxLine, LocusAlleleCounts, Window, PhenotypeLine, Phenotype
# include("functions_io.jl")
# using ProgressMeter, Distributions
#####################

### FILTERING FUNCTIONS
function FILTER(line::Union{PileupLine, SyncxLine}, maximum_missing_fraction::Float64, alpha1::Float64=0.05, maf::Float64=0.01, alpha2::Float64=0.50, minimum_coverage::Int64=5)::Tuple{Bool, Window}
    ####### TEST ########
    # filename = "../test/test.pileup"
    # file = open(filename, "r")
    # line = PileupLine(1, readline(file)); close(file)
    # # filename = "../test/test.syncx"
    # # file = open(filename, "r")
    # # line = SyncxLine(1, readline(file)); close(file)
    # close(file)
    # maximum_missing_fraction = 0.90
    # alpha1 = 0.05
    # maf = 0.01
    # alpha2 = 0.50
    # minimum_coverage = 2
    #####################
    window = PARSE([PARSE(line)])
    ### copy the allele counts and set missing to zero
    X = copy(window.cou)
    X[ismissing.(X)] .= 0
    ### calculate the coverage per pool and count the number of zero coverga pools and the maximum number of zero coverage pools
    coverages = sum(X, dims=1)
    zero_cov_count = sum(coverages[1,:] .== 0)
    maximum_zero_cov_count = Int(round(length(coverages[1,:])*maximum_missing_fraction))
    ### Filter-out zero coverage pools for the calculation of minimum allele frequency and minimum coverage thresholds
    idx = coverages[1,:] .> 0
    X = X[:, idx]
    ### Return false if there are only at most a single non-zero coverage pool
    if prod(size(X)) > 7
        coverages = coverages[:, idx]
        frequencies = X ./ coverages
        frequencies[isnan.(frequencies)] .= 0.0
        _, p = size(frequencies)
        ### allele frequencies (alpha1)
        idx = collect(1:7)[(sum(frequencies, dims=2) .== maximum(sum(frequencies, dims=2)))[:,1]][1]
        frq = sort(frequencies[idx,:])
        emaxbf = frq[(end-Int(ceil(alpha1*p)))] ### expected maximum biallelic frequency at alpha1
        eminbf = 1 - emaxbf                     ### expected minimum biallelic frequency at alpha1
        ### coverage (alpha2)
        cov = sort(coverages[1,:])[Int(ceil(alpha2*p))]
        ### output
        out = (zero_cov_count <= maximum_zero_cov_count) & (eminbf >= maf) & (emaxbf < (1-maf)) & (eminbf < 1) & (emaxbf > 0) & (cov >= minimum_coverage)
    else
        out = false
    end
    return(out, window)
end

function FILTER(pileup::String, outype::String, init::Int, term::Int, maximum_missing_fraction::Float64=0.10, alpha1::Float64=0.05, maf::Float64=0.01, alpha2::Float64=0.50, minimum_coverage::Int64=5, out::String="")::String
    ####### TEST ########
    # pileup = "../test/test.pileup"
    # outype = ["pileup", "syncx"][1]
    # file = open(pileup, "r")
    # seekend(file)
    # final_position = position(file)
    # close(file)
    # threads = 4
    # vec_positions = Int.(round.(collect(0:(final_position/threads):final_position)))
    # init = vec_positions[2]
    # term = vec_positions[3]
    # maximum_missing_fraction = 0.50
    # alpha1 = 0.05
    # maf = 0.01
    # alpha2 = 0.50
    # minimum_coverage = 2
    # out = ""
    #####################
    ### Output file
    if out == ""
        if outype == "pileup"
            out = string(join(split(pileup, ".")[1:(end-1)], "."), "-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, ".pileup")
        elseif outype == "syncx"
            out = string(join(split(pileup, ".")[1:(end-1)], "."), "-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, ".syncx")
        end
    end
    ### Open pileup file
    file = open(pileup, "r")
    ### Set position to the next line if the current position is a truncated line
    if (init>0)
        seek(file, init-1)
        line = readline(file)
        if line != ""
            init = position(file)
        end
    end
    seek(file, init)
    ### Filter
    while position(file) < term
        line = PileupLine(1, readline(file));
        test, window = FILTER(line, maximum_missing_fraction, alpha1, maf, alpha2, minimum_coverage)
        if test
            if outype=="pileup"
                SAVE(line, out)
            elseif outype=="syncx"
                SAVE(window, out)
            end
        end
    end
    close(file)
    return(out)
end

function FILTER(syncx::String, init::Int, term::Int, maximum_missing_fraction::Float64=0.10, alpha1::Float64=0.05, maf::Float64=0.01, alpha2::Float64=0.50, minimum_coverage::Int64=5, out::String="")::String
    ####### TEST ########
    # syncx = "../test/test.syncx"
    # file = open(syncx, "r")
    # seekend(file)
    # final_position = position(file)
    # close(file)
    # threads = 4
    # vec_positions = Int.(round.(collect(0:(final_position/threads):final_position)))
    # init = vec_positions[2]
    # term = vec_positions[3]
    # maximum_missing_fraction = 0.50
    # alpha1 = 0.05
    # maf = 0.01
    # alpha2 = 0.50
    # minimum_coverage = 2
    # out = ""
    #####################
    ### Output file
    if out == ""
       out = string(join(split(syncx, ".")[1:(end-1)], "."), "-FILTERED-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, ".syncx")
    end
    ### Set position to the next line if the current position is a truncated line
    file = open(syncx, "r")
    if (init>0)
        seek(file, init-1)
        line = readline(file)
        if line != ""
            init = position(file)
        end
    end
    seek(file, init)
    ### Filter
    while position(file) < term
        line = SyncxLine(1, readline(file));
        test, window = FILTER(line, maximum_missing_fraction, alpha1, maf, alpha2, minimum_coverage)
        if test
            SAVE(window, out)
        end
    end
    close(file)
    return(out)
end

function FILTER(χ::Window, maf::Float64=0.001, δ::Float64=1e-10, remove_insertions::Bool=true, remove_minor_alleles::Bool=false, remove_correlated_alleles::Bool=false, θ::Float64=1.00, centre::Bool=false)::Tuple{Matrix{Float64}, Vector{String}, Vector{Int64}, Vector{String}}
    ####### TEST ########
    # χ = LOAD("../test/test.syncx", false)
    # maf = 0.001
    # δ = 1e-10
    # remove_insertions = true
    # remove_minor_alleles = false
    # remove_correlated_alleles = true
    # θ = 0.99
    # centre = true
    #####################
    X = Float64.(χ.cou')
    vec_idx = []
    ### If we are outputting frequencies we need to remove insertion counts since they're also counted as the bases the are.
    if remove_insertions
        for i in 1:7:size(X,2)
            # i = 1
            g = X[:, i:(i+6)]
            g[:, 5] .= 0
            X[:, i:(i+6)] = g ./ sum(g, dims=2)
        end
    end
    ### Keep only p-1 alleles per locus where p is the number of polymorphic alleles in a locus
    if remove_minor_alleles
        S = reshape(sum(X, dims=1), 7, Int(size(X,2)/7))
        n, m = size(S)
        for j in 1:m
            # j = 1
            x = S[:, j]
            idx_nonzero = x .!= 0
            x = x[idx_nonzero]
            idx_sort = sortperm(x)
            x = x[idx_sort]
            idx_minus_one_allele = collect(2:length(x))

            idx = collect(1:7)[idx_nonzero][idx_sort][idx_minus_one_allele]
            idx = (j-1)*7 .+ idx

            append!(vec_idx, idx)
        end
        X = X[:, vec_idx]
    else
        vec_idx = collect(1:size(X,2))
    end
    ### Filter by minimum allele frequency
    if (maf > 0) && (maf < 1)
        # println("Filtering by minor allele frequency.")
        idx = (minimum(X, dims=1)[1,:] .>= maf) .& (maximum(X, dims=1)[1,:] .<= 1-maf)
        X = X[:, idx]
        vec_idx = vec_idx[idx] ### update index
    end
    ### Remove non-polymorphic loci
    # println("Removing non-polymorphc loci.")
    idx = (Distributions.var(X, dims=1) .> δ)[1,:]
    X = X[:, idx]
    vec_idx = vec_idx[idx] ### update index
    p = length(vec_idx)
    ### Remove fixed alleles, and collinear alleles
    ### Iteratively calculate the correlations between alleles across loci which should be more efficient as we are testing and breaking if the correlation is greater than or equal to θ
    if remove_correlated_alleles
        idx = zeros(Bool, p)
        idx[end] = true
        for i in 1:(p-1)
            test_θ = true
            for j in (i+1):p
                # i = 1; j = 10
                if abs(cor(X[:, i], X[:, j])) <= θ
                    test_θ *= true
                else
                    test_θ *= false
                    break
                end
            end
            if test_θ
                idx[i] = true
            end
        end
        vec_idx = vec_idx[idx] ### update index
        X = X[:, idx]
    end
    ### centre on 0, i.e. freq=0 becomes -0.5 and freq=1 becomes +0.5
    if centre
        X = X .- 0.5
    end
    ### Output
    vec_chr = repeat(χ.chr, inner=7)[vec_idx]
    vec_pos = repeat(χ.pos, inner=7)[vec_idx]
    vec_ale = repeat(["A", "T", "C", "G", "INS", "DEL", "N"], outer=Int(size(χ.cou,1)/7))[vec_idx]
    return(X, vec_chr, vec_pos, vec_ale)
end

### TRANSFORMATION FUNCTIONS
function PILEUP2SYNCX(pileup::String, init::Int, term::Int, out::String="")::String
    ####### TEST ########
    # pileup = "../test/test.pileup"
    # file = open(pileup, "r")
    # seekend(file)
    # final_position = position(file)
    # close(file)
    # threads = 4
    # vec_positions = Int.(round.(collect(0:(final_position/threads):final_position)))
    # init = vec_positions[2]
    # term = vec_positions[3]
    # out = ""
    #####################
    ### Output file
    if out == ""
        out = string(join(split(pileup, ".")[1:(end-1)], "."), "-CONVERTED.syncx")
    end
    ### Open pileup file
    file = open(pileup, "r")
    ### Set position to the next line if the current position is a truncated line
    if (init>0)
        seek(file, init-1)
        line = readline(file)
        if line != ""
            init = position(file)
        end
    end
    seek(file, init)
    while position(file) < term
        SAVE(PARSE([PARSE(PileupLine(1, readline(file)))]), out)
    end
    close(file)
    return(out)
end

function TRANSFORM(y::Vector{T}, absolute_threshold::Float64=0.1, maxiter::Int=10)::Tuple{Vector{T}, Vector{String}, Vector{T}} where T <: Number
    ####### TEST ########
    # using UnicodePlots
    # y = rand(100)
    # UnicodePlots.histogram(y)
    # absolute_threshold = 0.1
    # maxiter = 10
    #####################
    vec_fun = []
    vec_max = []
    iter = 0
    while ((skewness(y) > absolute_threshold) || (skewness(y) < -absolute_threshold)) && (iter < maxiter)
        iter += 1
        min_abs = abs(minimum(vcat(y, 0.0)))
        if skewness(y) > 1.00
            push!(vec_fun, "++")
            push!(vec_max, min_abs)
            y = log.(y .+ vec_max[end])
        elseif skewness(y) > 0.00
            push!(vec_fun, "+")
            push!(vec_max, min_abs)
            y = sqrt.(y .+ vec_max[end])
        elseif skewness(y) < -1.00
            push!(vec_fun, "--")
            push!(vec_max, maximum(y .+ min_abs))
            y = log.(vec_max[end] .- y)
        elseif skewness(y) < 0.00
            push!(vec_fun, "-")
            push!(vec_max, maximum(y .+ min_abs))
            y = sqrt.(vec_max[end] .- y)
        else
            break
        end
    end
    return(y, vec_fun, vec_max)
end

function UNTRANSFORM(y::Vector{T}, vec_fun::Vector{String}, vec_max::Vector{T})::Vector{T} where T <: Number
    ####### TEST ########
    # using UnicodePlots
    # y_untransformed = rand(100)
    # y, vec_fun, vec_max = TRANSFORM(y)
    # UnicodePlots.histogram(y_untransformed)
    # UnicodePlots.histogram(y)
    #####################
    vec_fun = reverse(vec_fun)
    vec_max = reverse(vec_max)
    for i in 1:length(vec_fun)
        if vec_fun[i] == "++"
            y = exp.(y) .- 1
        elseif vec_fun[i] == "+"
            y = y.^2
        elseif vec_fun[i] == "--"
            y = .-(exp.(y) .- vec_max[i])
        elseif vec_fun[i] == "-"
            y = .-((y.^2) .- vec_max[i])
        end
    end
    return(y)
end
