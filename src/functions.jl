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

include("structs.jl")
using .structs: PileupLine, SyncxLine, LocusAlleleCounts, Window, PhenotypeLine, Phenotype

### DATA PARSING AND EXTRACTION
function PARSE(line::PileupLine, minimum_quality=20)::LocusAlleleCounts
    lin = split(line.lin, '\t')
    chr = lin[1]
    pos = parse(Int, lin[2])
    ref = lin[3][1]
    dep = parse.(Int, lin[4:3:end])
    sta = lin[5:3:end]
    qua = lin[6:3:end]
    p = length(dep)
    A=[]; T=[]; C=[]; G=[]; I=[]; D=[]; N=[]
    for i in 1:p
        vec_sta_per_pool = split(replace(replace(sta[i], Regex("\\^.") => ""), "\$" => ""), "") ### priot to splitting - remove alignment start and mapping quality strings, as well as alignment end marker
        vec_qua_per_pool = split(qua[i], "")
        j = k = 0
        push!(A, 0); push!(T, 0); push!(C, 0); push!(G, 0); push!(I, 0); push!(D, 0); push!(N, 0)
        while j < length(vec_sta_per_pool)
            j += 1
            k += 1
            # println(string("i=", i))
            # println(string("j=", j))
            # println(string("k=", k))
            
            ### allele read qualities
            q = try
                    if vec_qua_per_pool[k][1] != '*'
                        Int(vec_qua_per_pool[k][1]) - 33
                    else 
                        0
                    end
                catch
                    0
                end

            ### allele states
            str_state = uppercase(vec_sta_per_pool[j])
            if (str_state == "+") | (str_state == "-")
                ### remove insertion and deletion sequences
                if str_state == "+"
                    # push!(s, "I")
                    q>=minimum_quality ? I[i] += 1 : nothing
                else
                    # push!(s, 'D')
                    D[i] += 1 ### no need to test for quality because the allele is deleted
                    k -= 1
                end
                l = 1
                count = parse(Int, vec_sta_per_pool[j+l])
                while count != "error"
                    l += 1
                    count = try
                        parse(Int, vec_sta_per_pool[j+l])
                    catch
                        "error"
                    end
                end
                j = j + l + parse(Int, string(vec_sta_per_pool[(j+1):(j+(l-1))]...))
            elseif ((str_state==".") | 
                    (str_state==",") | 
                    (str_state ∈ ["A", "T", "C", "G"])
                   ) & (q>=minimum_quality)
                if str_state ∈ ["A", "T", "C", "G"]
                    a = uppercase(str_state)[1]
                else
                    a = uppercase(ref)[1]
                end
                if a == 'A'
                    A[i] += 1
                elseif a == 'T'
                    T[i] += 1
                elseif a == 'C'
                    C[i] += 1
                elseif a == 'G'
                    G[i] += 1
                else
                    N[i] += 1
                end
            else
                # push!(s, 'N')
                N[i] += 1
            end
            
        end

    end
    return(LocusAlleleCounts(chr, pos, ref, dep, A, T, C, G, I, D, N))
end

function PARSE(line::SyncxLine)::LocusAlleleCounts
    # file = open("/home/jeffersonfparil/Documents/poolgen/test/test_2.syncx", "r")
    # line = SyncxLine(1, readline(file))
    lin = split(line.lin, "\t")
    chr = lin[1]
    pos = parse(Int, lin[2])
    vco = map(x -> x=="missing" ? missing : parse(Int, x), vcat(split.(lin[3:end], ":")...))
    cou = (convert(Matrix{Any}, reshape(vco, 7, Int(length(vco)/7))))
    cou[ismissing.(cou)] .= 0
    ref = 'N'
    dep = sum(cou, dims=1)[1,:]
    A = cou[1,:]
    T = cou[2,:]
    C = cou[3,:]
    G = cou[4,:]
    I = cou[5,:]
    D = cou[6,:]
    N = cou[7,:]
    return(LocusAlleleCounts(chr, pos, ref, dep, A, T, C, G, I, D, N))
end

function PARSE(window::Vector{LocusAlleleCounts})::Window
    n = length(window)
    p = length(window[1].dep)
    chr = []
    pos = []
    ref = []
    cou = Array{Any,2}(missing, (n*7, p))
    for i in 1:n
        push!(chr, window[i].chr)
        push!(pos, window[i].pos)
        push!(ref, window[i].ref)
        idx = window[i].dep .> 0
        cou[((i-1)*7)+1, idx] = window[i].A[idx]
        cou[((i-1)*7)+2, idx] = window[i].T[idx]
        cou[((i-1)*7)+3, idx] = window[i].C[idx]
        cou[((i-1)*7)+4, idx] = window[i].G[idx]
        cou[((i-1)*7)+5, idx] = window[i].I[idx]
        cou[((i-1)*7)+6, idx] = window[i].D[idx]
        cou[((i-1)*7)+7, idx] = window[i].N[idx]
    end
    imp = zeros((n*7, p))
    return(Window(chr, pos, ref, cou, imp))
end

function PARSE(line::PhenotypeLine, trait_names::Vector{String}=[""], missing_strings::Vector{String}=["NA", ""])::Phenotype
    # phenotype = "/home/jeffersonfparil/Documents/poolgen/test/test.csv"
    # file = open(phenotype, "r")
    # header = PhenotypeLine(1, readline(file), ",", 1, [2,3,4,5,10,11])
    # line = PhenotypeLine(1, readline(file), ",", 1, [2,3,4,5,10,11])
    # close(file)
    # missing_strings = ["NA", ""]
    lin = split(line.lin, line.dlm)
    y = map(x -> sum(x .== missing_strings)>0 ? missing : parse(Float64, x), lin[line.trc])
    ids = [string(lin[line.idc])]
    if trait_names==[""]
        trait_names = repeat([""], length(y))
    end
    Y = convert(Array{Any}, reshape(y, 1, length(y)))
    return(Phenotype(ids, trait_names, Y))
end

function PARSE(lines::Vector{Phenotype}, rename_traits::Vector{String}=[""])::Phenotype
    # phenotype = "/home/jeffersonfparil/Documents/poolgen/test/test.csv"
    # file = open(phenotype, "r")
    # header = PhenotypeLine(1, readline(file), ",", 1, [2,3,4,5,10,11])
    # line1 = PARSE(PhenotypeLine(1, readline(file), ",", 1, [2,3,4,5,10,11]))
    # line2 = PARSE(PhenotypeLine(1, readline(file), ",", 1, [2,3,4,5,10,11]))
    # line3 = PARSE(PhenotypeLine(1, readline(file), ",", 1, [2,3,4,5,10,11]))
    # lines = [line1, line2, line3]
    # close(file)
    # rename_traits = ["a", "b", "c", "d", "e", "f"]
    ids = lines[1].iid
    trait_names = lines[1].tid
    Y = lines[1].phe
    for i in 2:length(lines)
        # i = 1
        if trait_names == lines[i].tid
            append!(ids, lines[i].iid)
            Y = vcat(Y, lines[i].phe)
        else
            continue
        end
    end
    if rename_traits != [""]
        trait_names = rename_traits
    end
    @assert (length(ids), length(trait_names)) == size(Y)
    return(Phenotype(ids, trait_names, Y))
end

function EXTRACT(window::Window, locus::Int)::Window
    Window([window.chr[locus]],
           [window.pos[locus]],
           [window.ref[locus]],
           window.cou[(7*(locus-1))+1:(locus*7), :],
           window.imp[(7*(locus-1))+1:(locus*7), :])
end

function EXTRACT(window::Window, loci::UnitRange{Int})::Window
    Window(window.chr[loci],
           window.pos[loci],
           window.ref[loci],
           window.cou[(7*(loci.start-1))+1:(loci.stop*7), :],
           window.imp[(7*(loci.start-1))+1:(loci.stop*7), :])
end

function SLIDE!(window::Window, locus::LocusAlleleCounts)::Window
    new = PARSE([locus])
    window.chr[1:(end-1)] = window.chr[2:end]
    window.pos[1:(end-1)] = window.pos[2:end]
    window.ref[1:(end-1)] = window.ref[2:end]
    window.cou[1:(end-7), :] = window.cou[8:end, :]
    window.imp[1:(end-7), :] = window.imp[8:end, :]
    window.chr[end] = new.chr[1]
    window.pos[end] = new.pos[1]
    window.ref[end] = new.ref[1]
    window.cou[(end-6):end, :] = new.cou[1:7, :]
    window.imp[(end-6):end, :] = zeros(size(new.cou))
    return(window)
end

### I/O
function SPLIT(threads::Int, pileup::String)::Tuple{Vector{Int64}, Vector{Int64}}
    # threads = 4
    # pileup = "/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup"
    ####
    file = open(pileup, "r")
    seekend(file)
    temp_positions = Int.(round.(collect(0:(position(file)/threads):position(file))))
    close(file)

    positions_init = [temp_positions[1]]
    positions_term = []

    file = open(pileup, "r")

    for i in 2:(length(temp_positions)-1)
        ### Navigate to the next line if the current line is truncated from the starting side
        seek(file, temp_positions[i]-1)
        line = readline(file)
        if line == ""
            seek(file, temp_positions[i])
        end
        push!(positions_init, position(file))
        push!(positions_term, position(file))
    end
    push!(positions_term, temp_positions[end])
    close(file)
    return(Int.(positions_init), Int.(positions_term))
end

function SPLIT(threads::Int, pileup::String, window_size::Int)::Tuple{Vector{Int64}, Vector{Int64}}
    # threads = 4
    # pileup = "/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup"
    # window_size = 20
    ####
    file = open(pileup, "r")
    seekend(file)
    temp_positions = Int.(round.(collect(0:(position(file)/threads):position(file))))
    close(file)

    positions_init = [temp_positions[1]]
    positions_term = []

    file = open(pileup, "r")

    for i in 2:(length(temp_positions)-1)
        ### Navigate to the next line if the current line is truncated from the starting side
        seek(file, temp_positions[i]-1)
        line = readline(file)
        if line == ""
            seek(file, temp_positions[i])
        end
        push!(positions_init, position(file))
        for j in 1:window_size
            readline(file)
        end
        push!(positions_term, position(file))
    end
    push!(positions_term, temp_positions[end])
    close(file)
    return(Int.(positions_init), Int.(positions_term))
    ### TEST
    # file = open(pileup, "r")
    # n_lines = []
    # chr_pos = []
    # for i in eachindex(positions_init)
    #     # i = 1
    #     seek(file, positions_init[i])
    #     n = 0
    #     while position(file) < positions_term[i]
    #         pos = parse(Int, split(readline(file), "\t")[2])
    #         n+=1
    #         push!(chr_pos, pos)
    #     end
    #     push!(n_lines, n)
    # end
    # close(file)
    # hcat(positions_init, positions_term, n_lines, cumsum(n_lines))

    # for i in 1:(length(cumsum(n_lines))-1)
    #     # i = 1
    #     pos = cumsum(n_lines)[i]
    #     tail = chr_pos[(pos-(window_size-1)):pos]
    #     head = chr_pos[(pos+1):(pos+window_size)]
    #     if tail == head
    #         println(string("The tail of chunk ", i, " matches the head of chunk ", i+1, "!"))
    #     end
    # end
end

function MERGE(filenames_out::Vector{String}, out::String)::String
    ### Sort the output files from parallel processing
    sort!(filenames_out)
    ### Merge
    file_out = open(out, "w")
    for i in eachindex(filenames_out)
        if isfile(filenames_out[i])
            file_in = open(filenames_out[i], "r")
            while !eof(file_in)
                write(file_out, string(readline(file_in), "\n"))
            end
            close(file_in)
            ### clean up
            rm(filenames_out[i])
        end
    end
    close(file_out)
    return(out)
end

function MERGE(filenames_out::Vector{String}, window_size::Int64, out::String)::String
    ### Sort the chunks so we maintain the one-to-one correspondence between input and output loci arrangement
    sort!(filenames_out)
    ### Trim-off overhanging windows and merge
    file_out = open(out, "w")
    for i in eachindex(filenames_out)
        if isfile(filenames_out[i])
            if i < length(filenames_out)
                ### trim trailing window from the first and intervening chunks
                lines = 0
                file_in = open(filenames_out[i], "r")
                while !eof(file_in)
                    lines += 1
                    readline(file_in);
                end
                close(file_in)
                max_line = lines - window_size
            else
                ### do not trim the last chunk
                max_line = Inf
            end
            file_in = open(filenames_out[i], "r")
            j = 0
            while (!eof(file_in)) & (j < max_line)
                j += 1
                write(file_out, string(readline(file_in), "\n"))
            end
            close(file_in)
            ### clean up
            rm(filenames_out[i])    ### syncx chunks
        end
    end
    close(file_out)
    return(out)
end

function SAVE(line::PileupLine, filename::String)
    file = open(filename, "a")
    write(file, string(line.lin, '\n'))
    close(file)
end

function SAVE(window::Window, filename::String)
    OUT = hcat(window.chr, window.pos)
    n, m = size(window.cou)
    for i in 1:m
        # i = 1
        counts = reshape(window.cou[:,i], (7, Int(n/7)))
        OUT = hcat(OUT, [join(x,':') for x in eachcol(counts)])
    end
    out = join([join(x,'\t') for x in eachrow(OUT)], '\n')
    file = open(filename, "a")
    write(file, string(out, '\n'))
    close(file)
end

function PILEUP2SYNCX(pileup::String, init::Int, term::Int, out::String="")::String
    # pileup = "/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup"
    # # pileup = "/home/jeffersonfparil/Documents/poolgen/test/test_2.pileup"
    # file = open(pileup, "r")
    # seekend(file)
    # final_position = position(file)
    # close(file)
    # threads = 4
    # vec_positions = Int.(round.(collect(0:(final_position/threads):final_position)))
    # init = vec_positions[2]
    # term = vec_positions[3]
    # out = ""
    ### Output file
    if out == ""
        out = string(join(split(pileup, ".")[1:(end-1)], "."), ".syncx")
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

function LOAD(phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_cols::Vector{Int}=[0], missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""])::Phenotype
    # phenotype = "/home/jeffersonfparil/Documents/poolgen/test/test.csv"
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_cols = [0]
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    file = open(phenotype, "r")
    if header
        header_line = convert(Vector{String}, split(readline(file), delimiter))
    end
    lines = []
    while !eof(file)
        l = readline(file)
        m = length(split(l, delimiter))
        if phenotype_cols==[0]
            phenotype_cols = collect(1:m)[collect(1:m) .!= id_col]
        end
        push!(lines, PARSE(PhenotypeLine(1, l, delimiter, id_col, phenotype_cols), [""], missing_strings))
    end
    close(file)
    return(PARSE(convert(Vector{Phenotype}, lines), header_line[phenotype_cols]))
end

function LOAD(syncx::String, count::Bool)::Window
    # include("user_functions.jl"); using .user_functions: pileup2syncx
    # syncx = pileup2syncx("/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup")
    # syncx = pileup2syncx("/home/jeffersonfparil/Documents/poolgen/test/test_2.pileup")
    # count = false
    file = open(syncx, "r")
    loci = []
    while !eof(file)
        push!(loci, PARSE(SyncxLine(1, readline(file))))
    end
    close(file)
    GENOTYPE = PARSE(convert(Vector{LocusAlleleCounts}, loci))
    if !count
        for i in 1:7:size(GENOTYPE.cou,1)
            # i = 1
            g = GENOTYPE.cou[i:(i+6), :]
            GENOTYPE.cou[i:(i+6), :] = g ./ sum(g, dims=1)
        end
    end
    return(GENOTYPE)
end

### FILTER
function FILTER(line, maximum_missing_fraction::Float64, alpha1::Float64=0.05, maf::Float64=0.01, alpha2::Float64=0.50, minimum_coverage::Int64=5)::Tuple{Bool, Window}
    # filename = "/home/jeffersonfparil/Documents/poolgen/test/test_2.pileup"
    # # filename = "/home/jeffersonfparil/Documents/poolgen/test/test_2.syncx"
    # file = open(filename, "r")
    # line = PileupLine(1, readline(file)); close(file)
    # # line = SyncxLine(1, readline(file)); close(file)
    # close(file)
    # maximum_missing_fraction = 0.90
    # alpha1 = 0.05
    # maf = 0.01
    # alpha2 = 0.50
    # minimum_coverage = 2
    ## parse input line: pileup or syncx
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
    ### Return false if there are no non-zero coverage pools
    if prod(size(X)) > 0
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
    # pileup = "/home/jeffersonfparil/Documents/poolgen/test/test_2.pileup"
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
    # include("user_functions.jl"); using .user_functions: pileup2syncx
    # # syncx = pileup2syncx("/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup")
    # syncx = pileup2syncx("/home/jeffersonfparil/Documents/poolgen/test/test_2.pileup")
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

function FILTER(phenotype::Phenotype, maximum_missing_fraction::Float64, alpha1::Float64, alpha2::Float64)
    PHENOTYPE = LOAD("/home/jeffersonfparil/Documents/poolgen/test/test_1.csv", ",")
end

### IMPUTE
function IMPUTE!(window::Window, model::String=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance::Bool=true)::Window
    n, p = size(window.cou)
    ### Find the indices of pools with missing data.
    ### These will be used independently and iteratively as our response variables
    idx_pools = sum(ismissing.(window.cou), dims=1)[1,:] .> 0
    ### If we have at least one pool with no missing data, then we proceed with imputation
    if sum(.!idx_pools) >= 1
        ### Explanatory variables
        X = Int.(window.cou[:, .!idx_pools])

         ### Distance covariate (only add if the window is within a single chromosome)
         if distance & (length(unique(window.chr))==1)
            m = length(window.pos)
            D = zeros(Int, m, m)
            for i in 1:m
                for j in 1:m
                    D[i,j] = abs(window.pos[i] - window.pos[j])
                end
            end
            Z = MultivariateStats.projection(MultivariateStats.fit(PCA, repeat(D, inner=(7,1)); maxoutdim=3)) ### using the first 3 PCs by default
            X = hcat(X, Z)
        end

        for j in collect(1:p)[idx_pools]
            # j = collect(1:p)[idx_pools][1]
            y = window.cou[:, j]
            idx_loci = ismissing.(y)
            y_train = Int.(y[.!idx_loci])
            X_train = X[.!idx_loci, :]
            nf, pf = size(X_train)
           
            ### Train models
            if model == "Mean"
                β = append!([0.0], repeat([1/pf], pf))
            elseif model == "OLS"
                β = try
                    hcat(ones(nf), X_train) \ y_train
                catch
                    try
                        LinearAlgebra.pinv(hcat(ones(nf), X_train)'*hcat(ones(nf), X_train)) * (hcat(ones(nf), X_train)'*y_train)
                    catch
                        missing
                    end
                end
            elseif (model == "RR") | (model == "LASSO") | (model == "GLMNET")
                model=="RR" ? alpha1=0 : model=="LASSO" ? alpha1=1 : alpha1=0.5
                β = try
                    try
                        GLMNet.coef(GLMNet.glmnetcv(hcat(ones(nf), X_train), y_train, alpha1=alpha1, tol=1e-7)) # equivalend to mod.path.betas[:, argmin(mod)]
                    catch
                        GLMNet.coef(GLMNet.glmnetcv(hcat(ones(nf), X_train), y_train, alpha1=alpha1, tol=1e-3)) # equivalend to mod.path.betas[:, argmin(mod)]
                    end
                catch
                    β = missing
                end
            end

            ### Impute
            if !ismissing(β)
                X_valid = X[idx_loci, :]
                y_imputed = Int.(round.(hcat(ones(sum(idx_loci)), X_valid) * β))
                y_imputed[y_imputed .< 0] .= 0 ### collapse negative counts to zero (negative imputations are only a problem on OLS)
                # y_imputed .-= minimum(y_imputed)
                ### use the average imputed value
                window.imp[idx_loci, j] .+= 1
                y_imputed_mean = append!([], ((window.cou[idx_loci, j] .* window.imp[idx_loci, j]) .+ y_imputed) ./ (window.imp[idx_loci, j] .+ 1))
                y_imputed_mean[ismissing.(y_imputed_mean)] = y_imputed
                window.cou[idx_loci, j] = Int.(round.(y_imputed_mean))
            end
        end
    else
        ### If don't have a single pool with no missing data, then we return the input window without imputing
        nothing
    end
    return(window)
end

function IMPUTE(syncx::String, init::Int, term::Int, window_size::Int=100, model::String=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance::Bool=true, out::String="")::String
    # syncx = "/home/jeffersonfparil/Documents/poolgen/test/test_1.syncx"
    # # syncx = "/home/jeffersonfparil/Documents/poolgen/test/test_2.syncx"
    # file = open(syncx, "r")
    # seekend(file)
    # threads = 4
    # vec_positions = Int.(round.(collect(0:(position(file)/threads):position(file))))
    # close(file)
    # init = vec_positions[2]
    # term = vec_positions[3]
    # window_size = 20
    # model = "LASSO"
    # distance = true
    # out = ""
    ### Output syncx
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-IMPUTED.syncx")
    end

    ### Open syncx file
    file = open(syncx, "r")
    ### Set position to the next line if the current position is a truncated line
    if (init>0)
        seek(file, init-1)
        line = readline(file)
        if line != ""
            init = position(file)
        end
    end
    seek(file, init)
    ### Impute
    i = init
    j = 0
    window = []
    while (i < term) & (!eof(file))
        if window == []
            while (j < window_size) & (i < term) & (!eof(file))
                j += 1
                locus = PARSE(SyncxLine(j, readline(file)))
                i = position(file)
                push!(window, locus)
            end
            window = PARSE(Array{LocusAlleleCounts}(window))
            IMPUTE!(window, model, distance)
            SAVE(EXTRACT(window, 1), out)
        end
        if !eof(file)
            j += 1
            line = SyncxLine(j, readline(file))
            i = position(file)
            locus = try
                        PARSE(line)
                    catch
                        break
                    end
            SLIDE!(window, locus)
            SAVE(EXTRACT(window, 1), out)
            IMPUTE!(window, model, distance)
        end
    end
    if !eof(file)
        SAVE(EXTRACT(window, 2:window_size), out)
    else
        SAVE(window, out)
    end
    close(file)
    return(out)
end

# ### MODEL
# include("/home/jeffersonfparil/Documents/poolgen/src/user_functions.jl")
# using .user_functions: pileup2syncx
# GENOTYPE = LOAD(pileup2syncx("/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup"), false)
# PHENOTYPE = LOAD("/home/jeffersonfparil/Documents/poolgen/test/test_1.csv", ",")

# X = Float64.(GENOTYPE.cou')
# Y = Float64.(PHENOTYPE.phe[:,1])
# Y = 

# B = X \ Y[:,1]


### MISC: ### SPARSITY SIMULATION AND CROSS-VALIDATION
function SIMULATESPARSITY(filename, read_length::Int=100, missing_loci_fraction::Float64=0.50, missing_pools_fraction::Float64=0.25, pileup_simulated_missing::String="")
    ###########################################################
    ### TEST
    # filename = "/data-weedomics-1/test_human.pileup"
    # read_length = 150
    # missing_loci_fraction = 0.50
    # missing_pools_fraction = 0.25
    # pileup_simulated_missing = ""
    ###########################################################
    println("Simulating missing data.")
    if pileup_simulated_missing==""
        pileup_simulated_missing = string(join(split(filename, '.')[1:(end-1)], '.'), "-SIMULATED_MISSING.pileup")
    end
    println("Counting lines.")
    @show loci_count = countlines(filename)
    println("Counting pools.")
    file_temp = open(filename, "r")
    i = -1
    if i < 0
        i += 1
        line = readline(file_temp)
    end
    close(file_temp)
    @show pool_count = Int((length(split(line, '\t')) - 3) / 3)
    println("Randomly sampling loci chunks to set to missing.")
    chunk_count = Int(ceil(loci_count / read_length))
    println("Counting the number of loci and pool which will be set to missing.")
    @show missing_loci_count = Int(round(chunk_count*missing_loci_fraction))
    @show missing_pool_count = Int(round(pool_count*missing_pools_fraction))
    ### Proceed if we will be simulating at least 1 missing datapoint
    if missing_loci_count*missing_pool_count > 0
        println("Randomly sample chunks of loci which will be set to missing.")
        random_chunk = sort(Random.randperm(chunk_count)[1:missing_loci_count])
        random_chunk = ((random_chunk .- 1) .* read_length) .+ 1
        println("Open input and output files, and initialise the iterators")
        FILE = open(filename, "r")
        file_out = open(pileup_simulated_missing, "w")
        i = 0
        j = 1
        println("Simulate missing data.")
        pb = ProgressMeter.Progress(loci_count, i)
        while !eof(FILE)
            i += 1; ProgressMeter.next!(pb)
            line = readline(FILE)
            ### while-looping to properly deal with the last read line
            while i == random_chunk[j]
                ### if we reach the last missing chunk, then stop incrementing
                length(random_chunk) == j ? j : j += 1
                ### parse the tab-delimited line
                vec_line = split(line, '\t')
                ### extract the scaffold or chromosome name
                scaffold_or_chromosome = vec_line[1]
                ### randomly choose pools to get missing data
                idx_pool_rand_missing = Random.randperm(pool_count)[1:missing_pool_count]
                idx_pool_rand_missing = (((idx_pool_rand_missing .- 1) .* 3) .+ 1) .+ 3
                position_ini = parse(Int, vec_line[2])
                position_end = position_ini + (read_length - 1) ### less initial position twice since we're counting the initial position as part of the read length and we've already written it before the forst iteration of the while-loop
                while (scaffold_or_chromosome == vec_line[1]) & (parse(Int, vec_line[2]) <= position_end) & !eof(FILE)
                    ### Set to missing each of the randomly sampled pools in the current locus
                    for k in idx_pool_rand_missing
                        vec_line[k:k+2] = ["0", "*", "*"]
                    end
                    ### Write-out the line with simulated missing data
                    line = join(vec_line, '\t')
                    write(file_out, string(line, '\n'))
                    i += 1; ProgressMeter.next!(pb)
                    line = readline(FILE)
                    vec_line = split(line, '\t')
                end
            end
            ### Write-out the line without missing data
            write(file_out, string(line, '\n'))
        end
        close(FILE)
        close(file_out)
        println("##############################################################")
        println("Missing data simulation finished! Please find the output file:")
        println(pileup_simulated_missing)
        println("##############################################################")
        return(pileup_simulated_missing)
    else
        println("Did not simulate any missing data. Because the fraction of missing loci and pools parameter was too low.")
        return(1)
    end
end

function CROSSVALIDATE(syncx_without_missing, syncx_with_missing, out, csv_out="")
    # syncx_without_missing = "test.syncx"
    # syncx_with_missing = "test-SIMULATED_MISSING.syncx"
    # out = "test-SIMULATED_MISSING-IMPUTED.syncx"
    # csv_out=""
    ### NOTE: we should have the same exact locus corresponding per row across these three files
    file = open(syncx_without_missing, "r")
    p = collect(eachindex(split(replace(readline(file), '\t'=>':'), ':'))) ### Number of columns delimited by tabs per pool and chromosome coordinate and delimited by colon between alleles
    close(file)
    ### output file
    if csv_out == ""
        csv_out = string("Imputation_cross_validation_output-", time(), ".csv")
    end
    ### extract expected and imputed allele counts
    file_without_missing = open(syncx_without_missing, "r")
    file_with_missing    = open(syncx_with_missing, "r")
    file_imputed         = open(out, "r")
    missing_counter = 0
    imputed_counter = 0
    while !eof(file_without_missing)
        c = split(replace(readline(file_without_missing), '\t'=>':'), ':')
        m = split(replace(readline(file_with_missing), '\t'=>':'), ':')
        i = split(replace(readline(file_imputed), '\t'=>':'), ':')

        if ((c[1] == m[1] == i[1]) & (c[2] == m[2] == i[2])) == false
            println("Error!")
            println(string(syncx_without_missing, ", ", syncx_with_missing, ", and ", out, " are sorted differently!"))
            println("Please sort the loci in the same order. Exiting now.")
            exit()
        end
        missings = (m .== "missing")
        unimputed = (i .== "missing")
        idx = (missings) .& (.!unimputed)
        missing_counter += sum(missings)
        if sum(idx) > 0
            imputed_counter += sum(missings)
            ### extract allele counts
            expected = parse.(Int, c[idx])
            imputed = parse.(Int, i[idx])
            ### calculate allele frequencies
            n = length(expected)
            _expected = Int.(reshape(expected, (7, Int(n/7))))
            _imputed = Int.(reshape(imputed, (7, Int(n/7))))
            expected_freq = vec(_expected ./ sum(_expected, dims=1))
            imputed_freq = vec(_imputed ./ sum(_imputed, dims=1))
            ### convert undefined quotients (0/0 = NaN) to zero
            expected_freq[isnan.(expected_freq)] .= 0.0
            imputed_freq[isnan.(imputed_freq)] .= 0.0
            ### save imputed locus data
            file_out = open(csv_out, "a")
            write(file_out, string(join([join(x, ',') for x in eachrow(hcat(expected, imputed, expected_freq, imputed_freq))], '\n'),'\n') )
            close(file_out)
        end
    end
    close(file_without_missing)
    close(file_with_missing)
    close(file_imputed)
    ### write out proportion of imputed missing data at the last line of the output file
    file_out = open(csv_out, "a")
    write(file_out, string(join([imputed_counter/missing_counter, "", "", ""], ','), "\n"))
    close(file_out)
    ### output
    return(csv_out)
end


end
