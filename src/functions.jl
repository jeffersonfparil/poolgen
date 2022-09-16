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
using LinearAlgebra
using Distributions
using Plots

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
    lin = map(x -> x=="0:0:0:0:0:0:0" ? "missing:missing:missing:missing:missing:missing:missing" : string(x), split(line.lin, "\t"))
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
        ### Set missing loci from missing:missing:missing:missing:missing:missing:missing into 0:0:0:0:0:0:1
        for j in 1:Int(n/7)
            if sum(ismissing.(counts[:, j])) > 0
                counts[:, j] = [0, 0, 0, 0, 0, 0, 0]
            end
        end
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
    # phenotype = "/home/jeffersonfparil/Documents/poolgen/test/test_3-pheno-any-filename.csv"
    # delimiter = ","
    # header = true
    # id_col = 3 ### pool IDs
    # phenotype_cols = [8, 9] ### colorimetric resistance metric, and survival rate
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
        push!(lines, PARSE(PhenotypeLine(1, l, delimiter, id_col, phenotype_cols), header_line[phenotype_cols], missing_strings))
    end
    close(file)
    return(PARSE(convert(Vector{Phenotype}, lines), header_line[phenotype_cols]))
end

function LOAD(syncx::String, count::Bool)::Window
    # include("user_functions.jl"); using .user_functions: pileup2syncx
    # syncx = pileup2syncx("/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup")
    # syncx = pileup2syncx("/home/jeffersonfparil/Documents/poolgen/test/test_2.pileup")
    # syncx = "/home/jeffersonfparil/Documents/poolgen/test/test_3.syncx"
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
            g[5, :] .= 0 ### And since we are outputting frequencies we need to remove insertion counts since they're also counted as the bases the are.
            GENOTYPE.cou[i:(i+6), :] = g ./ sum(g, dims=1)
        end
    end
    return(GENOTYPE)
end

function LOAD(py_phenotype::String)::Tuple{String, Float64, Float64, Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    # py_phenotype = "/home/jeffersonfparil/Documents/poolgen/test/test_3-pheno-any-filename.py"
    file = open(py_phenotype, "r")
    phenotype_name = replace(replace(split(readline(file), "=")[end], "\"" => ""), ";" => "")
    sigma = parse(Float64, replace(replace(split(readline(file), "=")[end], "\"" => ""), ";" => ""))
    min = parse(Float64, replace(replace(split(readline(file), "=")[end], "\"" => ""), ";" => ""))
    max = parse(Float64, replace(replace(split(readline(file), "=")[end], "\"" => ""), ";" => ""))
    perc = parse.(Float64, split(split(split(readline(file), "[")[2], "]")[1], ","))
    q = parse.(Float64, split(split(split(readline(file), "[")[2], "]")[1], ","))
	# n = length(q) + 1
    bins = append!([x for x in perc], 1) - append!(zeros(1), perc)
    close(file)
    return(phenotype_name, sigma, min, max, perc, q, bins)
end

function CONVERT(syncx_or_sync::String, init::Int, term::Int, out::String="")::String
    # syncx_or_sync = "/home/jeffersonfparil/Documents/poolgen/test/test_3.syncx"
    # # syncx_or_sync = "/home/jeffersonfparil/Documents/poolgen/test/test_3-test.sync"
    # file = open(syncx_or_sync, "r")
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
        extension = [split(syncx_or_sync, ".")[end][end] == 'x' ? ".sync" : ".syncx"][1]
        out = string(join(split(syncx_or_sync, ".")[1:(end-1)], "."), extension)
    end
    ### Set position to the next line if the current position is a truncated line
    file = open(syncx_or_sync, "r")
    if (init>0)
        seek(file, init-1)
        line = readline(file)
        if line != ""
            init = position(file)
        end
    end
    seek(file, init)
    file_out = open(out, "a")

    p = length(split(split(readline(file), "\t")[end], ":"))
    ### Filter
    if p == 7
        ### SYNCX to SYNC
        seek(file, init)
        idx = [1, 2, 3, 4, 6, 7] ### Remove column 5, i.e. INSERTION COUNT COLUMN
        while position(file) < term
            line = split(readline(file), "\t")
            n = length(line) - 2
            X = parse.(Int, reshape(vcat(split.(line[3:end], ":")...), 7, n))'
            X = X[:, idx]
            for i in 1:n
                line[i+2] = join(X[i,:], ":")
            end
            line = vcat(line[1:2], "N", line[3:end])
            line = string(join(line, "\t"), "\n")
            write(file_out, line)
        end
    else
        ### SYNC to SYNCX
        seek(file, init)
        while position(file) < term
            line = split(readline(file), "\t")
            n = length(line) - 3
            X = parse.(Int, reshape(vcat(split.(line[4:end], ":")...), 6, n))'
            X = hcat(X[:,1:4], zeros(Int64, n), X[:,5:end])
            for i in 1:n
                line[i+3] = join(X[i,:], ":")
            end
            line = vcat(line[1:2], line[4:end]) ### remove reference allele column
            line = string(join(line, "\t"), "\n")
            write(file_out, line)
        end
    end
    close(file)
    close(file_out)
    return(out)
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

function FILTER(GENOTYPE::Window, epsilon::Float64=1e-10)::Tuple{Matrix{Float64}, Vector{Int64}}
    # GENOTYPE = LOAD("../test/test_3.syncx", false)
    X = Float64.(GENOTYPE.cou')
    ### Keep only p-1 alleles per locus where p is the number of polymorphic alleles in a locus
    _X = reshape(sum(X, dims=1), 7, Int(size(X,2)/7))
    vec_idx = []
    for j in axes(_X, 2)
        x = _X[:, j]
        idx_nonzeros = collect(1:7)[x .!= 0]
        idx = (j-1)*7 .+ idx_nonzeros[x[idx_nonzeros] .!= minimum(x[idx_nonzeros])]
        append!(vec_idx, idx)
    end
    X = X[:, vec_idx]
    ### Remove non-polymorphic loci
    idx = (Distributions.var(X, dims=1) .> epsilon)[1,:]
    X = X[:, idx]
    ### update index
    vec_idx = vec_idx[idx]
    ### centre on 0, i.e. freq=0 becomes -0.5 and freq=1 becomes +0.5
    X = (X .- 0.5)
    return(X, vec_idx)
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
            Z = MultivariateStats.projection(MultivariateStats.fit(PCA, repeat(D, inner=(7,1)); maxoutdim=1))
            # X = hcat(X, Z)
            X = hcat(X, X[:,1] .* Z[:,1])
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
            else
                println("Sorry ", model, " model is not implemented.")
                return(1)
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
    # # syncx = "/home/jeffersonfparil/Documents/poolgen/test/test_4.syncx"
    # file = open(syncx, "r")
    # seekend(file)
    # threads = 4
    # vec_positions = Int.(round.(collect(0:(position(file)/threads):position(file))))
    # close(file)
    # init = vec_positions[2] # init = 0
    # term = vec_positions[3] # file = open(syncx, "r"); seekend(file); term = position(file);  close(file)
    # window_size = 100
    # model = "OLS"
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

### Remove highly collinear alleles
function REMOVE_COLLINEAR_ALLELES_AND_LOCI(X::Matrix{Float64}, vec_idx::Vector{Float64}; cor_threshold::Float64=0.90)::Tuple{Matrix{Float64}, Vector{Float64}}
    # cor_threshold = 0.90
    idx = []
    @showprogress for j1 in axes(X, 2)
        vec_j2 = collect((j1+1):size(X, 2))
        for j2 in vec_j2
            #   j1 = 1
            #   j2 = 3
            if abs(cor(X[:,j1], X[:,j2])) > cor_threshold
                append!(idx, j2)
            end
        end
    end
    # idx_bk = copy(idx)
    idx = unique(sort(idx))
    X = X[:, idx]
    ### AGAIN! Update index
    vec_idx = vec_idx[idx]
    ### Output
    return(X, vec_idx)
end

function NLL_BETA(beta::Array{Float64,1}, data_A::Array{Float64,1}, data_B::Array{Float64,1})
	-sum(
		 log.(10,
		 	 (Distributions.cdf.(Distributions.Beta(beta[1], beta[2]), data_A)) .-
		 	 (Distributions.cdf.(Distributions.Beta(beta[1], beta[2]), append!(zeros(1), data_A[1:(length(data_A)-1)])))
		     )
		 )     -sum(
		 log.(10,
		 	 (Distributions.cdf.(Distributions.Beta(beta[3], beta[4]), data_B)) .-
		 	 (Distributions.cdf.(Distributions.Beta(beta[3], beta[4]), append!(zeros(1), data_B[1:(length(data_B)-1)])))
		     )
		 )
end

function GWALPHA(syncx::String, py_phenotype::String, init::Int64, term::Int64, out::String="")
    # syncx = "/home/jeffersonfparil/Documents/poolgen/test/test_3.syncx"
    # py_phenotype = "/home/jeffersonfparil/Documents/poolgen/test/test_3-pheno-any-filename.py"
    # file = open(syncx, "r")
    # seekend(file)
    # threads = 4
    # vec_positions = Int.(round.(collect(0:(position(file)/threads):position(file))))
    # close(file)
    # init = vec_positions[2] # init = 0
    # term = vec_positions[3] # file = open(syncx, "r"); seekend(file); term = position(file);  close(file)
    # out = ""
    ### Output syncx
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-GWAlpha.tsv")
    end
    ### Extract phenotype information
    phenotype_name, sigma, min, max, perc, q, bins = LOAD(py_phenotype)
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
    ### GWAlpha iteratively
    i = init
    while (i < term) & (!eof(file))
        i = position(file)
        locus = PARSE([PARSE(SyncxLine(j, readline(file)))])
        freqs = locus.cou ./ sum(locus.cou, dims=2)
        freqs[isnan.(freqs)] .= 0.0
	    allele_freqs = sum(freqs' .* bins, dims=1)
    end

end


function GWAlpha_ML_iterator(COUNTS::SharedArray{Int64,3}, snp::Int64, BINS::Array{Float64,1}, MAF::Float64, MIN::Float64, MAX::Float64, SD::Float64)
	OUT_ALPHA = convert(Array{Float64}, zeros(6))
	OUT_allele = zeros(6)
	OUT_snp = zeros(6)
	OUT_1 = zeros(6)
	OUT_pA = convert(Array{Float64}, zeros(6))
	# #parse allele counts from the sync file
	# COUNTS = zeros(Int64, NPOOLS, 6)
	# for i in 1:NPOOLS
	# 	COUNTS[i,:] = [parse(Int64, x) for x in split.(SYNC[snp, 4:(NPOOLS+3)], [':'])[i]]
	# end
	#convert to frequencies per pool
	counts = COUNTS[:,:,snp]
	FREQS = counts ./ ( sum(counts, dims=2) .+ 1e-10 ) #added 1e-10 to the denominator to avoid NAs in pools with no allele counts (zero depth; which should actually have been filtered out after mpileup using awk)
	allele_freqs = sum(FREQS .* BINS, dims=1)
	#iterate across alleles while filtering by MAF
	if (sum(counts) != 0.0)
		if (minimum(allele_freqs[allele_freqs .!= 0.0]) >= MAF) & (maximum(allele_freqs) < (1.0 - MAF)) #locus filtering by mean MAF
			for allele in 1:6
				if (allele_freqs[allele] > 0.0) & (maximum(FREQS[:,allele]) < 0.999999)  #allele filtering remove alleles with no counts and that the number of pools with allele frequency close to one should not occur even once!
				# if (sum(FREQS[:,allele] .== 0.0) < NPOOLS) & (sum(FREQS[:,allele] .> 0.999999) < 1) #filter-out alleles with at least 1 pool fixed for that allele because it causes a failure in the optimization
					freqA = FREQS[:, allele]
					pA = sum(freqA .* BINS)
					pB = 1 - pA
					BINA = (freqA .* BINS) ./ pA
					BINB = ( (1 .- freqA) .* BINS ) ./ (1-pA)
					percA = cumsum(BINA)
					percB = cumsum(BINB)

					### optimize (minimize) -log-likelihood of these major allele frequencies modelled as a beta distribution
					# using Nelder-Mead optimization or Box minimisation (try-catch if one or the other fails with preference to Nelder-Mead)
					lower_limits = [1e-20, 1e-20, 1e-20, 1e-20]
					upper_limits = [1.0, 1.0, 1.0, 1.0]
					initial_values = [0.1, 0.1, 0.1, 0.1]
					BETA = try
						Optim.optimize(beta->NLL_BETA(beta, percA, percB), initial_values, NelderMead())
					catch
						try
							Optim.optimize(beta->NLL_BETA(beta, percA, percB), lower_limits, upper_limits, initial_values)
						catch ### lower limits of 1e-20 to 1e-6 causes beta dist parameter values to shrink to zero somehow - so we're setting lower limits to 1e-5 instead
							lower_limits = [1e-5, 1e-5, 1e-5, 1e-5]
							Optim.optimize(beta->NLL_BETA(beta, percA, percB), lower_limits, upper_limits, initial_values)
						end
					end
					MU_A = MIN + ((MAX-MIN)*BETA.minimizer[1]/(BETA.minimizer[1]+BETA.minimizer[2]))
					MU_B = MIN + ((MAX-MIN)*BETA.minimizer[3]/(BETA.minimizer[3]+BETA.minimizer[4]))

					### compute alpha
					W_PENAL = 2*sqrt(pA*pB)
					ALPHA = W_PENAL*(MU_A - MU_B) / SD
					OUT_ALPHA[allele] =  ALPHA
					OUT_allele[allele] =  allele
					OUT_snp[allele] =  snp
					OUT_1[allele] =  1
					OUT_pA[allele] =  pA
				else
					## for explicit zero effects for null (low to none) frequency alleles
					OUT_ALPHA[allele] =  0.0
					OUT_allele[allele] =  allele
					OUT_snp[allele] =  snp
					OUT_1[allele] =  0
					OUT_pA[allele] =  allele_freqs[allele]
				end
			end
		end
	end
	return([OUT_ALPHA, OUT_allele, OUT_snp, OUT_1, OUT_pA])
end

function GWAlpha_ML(;filename_sync::String, filename_phen_py::String, MAF::Float64)
	### load the sync and phenotype files
	SYNC = DelimitedFiles.readdlm(filename_sync, '\t')
	phen = DelimitedFiles.readdlm(filename_phen_py, '\t')

	### gather phenotype specifications
	NPOOLS = length(split(phen[5], ['=', ',', '[', ']', ';'])) - 3 #less the first leading and trailing elements
	if length(split(phen[1], ['=', '\"'])) < 3
		global NAME = split(phen[1], ['=', '\''])[3]
	else
		global NAME = split(phen[1], ['=', '\"'])[3]
	end
	SD = parse.(Float64, split(phen[2], ['=',';'])[2])
	MIN = parse.(Float64, split(phen[3], ['=',';'])[2])
	MAX = parse.(Float64, split(phen[4], ['=',';'])[2])
	PERC = parse.(Float64, split(phen[5], ['=', ',', '[', ']', ';'])[3:(NPOOLS+1)])
	QUAN = parse.(Float64, split(phen[6], ['=', ',', '[', ']', ';'])[3:(NPOOLS+1)])
	BINS = append!([x for x in PERC], 1) - append!(zeros(1), PERC)

	### gather genotype (allele frequency) specifications
	NSNP = size(SYNC)[1]
	n_pools_sync = size(SYNC)[2] - 3
	if NPOOLS != n_pools_sync
		println("The number of pools with phenotype data does not match the number of pools with allele frequency data!")
		println("Please check you input files :-)")
		println("Remove leading and intervening whitespaces in the phenotype file.")
		exit()
	else
		n_pools_sync = nothing #clear out contents of this redundant n_pools variable
	end

	### creating an allele counts SharedArray from SYNC to minimize memory usage
	### input shared array
	COUNTS = SharedArrays.SharedArray{Int64,3}(NPOOLS, 6, size(SYNC)[1])
	progress_bar = ProgressMeter.Progress(NPOOLS, dt=1, desc="Converting the sync file into a SharedArray of allele counts: ",  barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
	for i in 1:NPOOLS
		COUNTS[i,:,:] = parse.(Int64, hcat(split.(SYNC[:, 4:(NPOOLS+3)], [':'])[:, i]...))
		ProgressMeter.update!(progress_bar, i)
	end
	### output shared arrays
	alpha_out = SharedArrays.SharedArray{Float64,1}(NSNP*6)
	allele_id = SharedArrays.SharedArray{Int64,1}(NSNP*6)
	locus_id = SharedArrays.SharedArray{Int64,1}(NSNP*6)
	locus_w_eff = SharedArrays.SharedArray{Int64,1}(NSNP*6)
	allele_freq = SharedArrays.SharedArray{Float64,1}(NSNP*6)

	if length(Distributed.procs()) > 1
		println("Performing GWAlpha_ML in parallel...")
	 	@time x = @sync @distributed for snp in 1:NSNP
			# println(snp)
			OUT = GWAlpha_ML_iterator(COUNTS, snp, BINS, MAF, MIN, MAX, SD)
			idx = collect(((6*snp)-5):(6*snp))
			alpha_out[idx] = OUT[1]
			allele_id[idx] = OUT[2]
			locus_id[idx] = OUT[3]
			locus_w_eff[idx] = OUT[4]
			allele_freq[idx] = OUT[5]
		end
	else
		for snp in 1:NSNP
			# println(snp)
			OUT = GWAlpha_ML_iterator(COUNTS, snp, BINS, MAF, MIN, MAX, SD)
			idx = collect(((6*snp)-5):(6*snp))
			alpha_out[idx] = OUT[1]
			allele_id[idx] = OUT[2]
			locus_id[idx] = OUT[3]
			locus_w_eff[idx] = OUT[4]
			allele_freq[idx] = OUT[5]
		end
	end
	ALPHA_OUT = alpha_out[locus_w_eff .== 1]
	ALLELE_ID_INT = allele_id[locus_w_eff .== 1]
	LOCUS_ID = locus_id[locus_w_eff .== 1]
	LOCUS_W_EFF = locus_w_eff[locus_w_eff .== 1]
	ALLELE_FREQ = allele_freq[locus_w_eff .== 1]

	### estimate heuristic p-values
	P_VALUES, LOD = significance_testing_module.estimate_pval_lod(convert(Array{Float64,1}, ALPHA_OUT))
	### output
	ALLELE_ID = repeat(["N"], inner=length(ALLELE_ID_INT))
	for i in 1:length(ALLELE_ID) #convert int allele ID into corresponding A, T, C, G, N, DEL
		if ALLELE_ID_INT[i] == 1; ALLELE_ID[i] = "A"
		elseif ALLELE_ID_INT[i] == 2; ALLELE_ID[i] = "T"
		elseif ALLELE_ID_INT[i] == 3; ALLELE_ID[i] = "C"
		elseif ALLELE_ID_INT[i] == 4; ALLELE_ID[i] = "G"
		elseif ALLELE_ID_INT[i] == 5; ALLELE_ID[i] = "N"
		elseif ALLELE_ID_INT[i] == 6; ALLELE_ID[i] = "DEL"
		end
	end
	OUT = (CHROM=convert(Array{Any,1},SYNC[LOCUS_ID,1]), POS=convert(Array{Int64,1},SYNC[LOCUS_ID,2]), ALLELE=convert(Array{Any,1},ALLELE_ID), FREQ=convert(Array{Float64,1},ALLELE_FREQ), ALPHA=convert(Array{Float64},ALPHA_OUT), PVALUES=convert(Array{Float64},P_VALUES), LOD=convert(Array{Float64},LOD))
	return(OUT)
end













function BEST_FITTING_DISTRIBUTION(vec_b::Vector{Float64})
    DIST_NAMES =   [Distributions.Bernoulli, Distributions.Beta, Distributions.Binomial, Distributions.Categorical,
                    Distributions.DiscreteUniform, Distributions.Exponential, Distributions.Normal, Distributions.Gamma,
                    Distributions.Geometric, Distributions.Laplace, Distributions.Pareto, Distributions.Poisson,
                    Distributions.InverseGaussian, Distributions.Uniform]
    DIST_INSTANCES = [try Distributions.fit_mle(D, vec_b); catch nothing; end for D in DIST_NAMES]
    NEG_LOGLIK = [try -sum(Distributions.logpdf.(D, vec_b)); catch nothing; end for D in DIST_INSTANCES]
    DISTRIBUTIONS_DF = hcat((DIST_NAMES[NEG_LOGLIK .!= nothing],
                            DIST_INSTANCES[NEG_LOGLIK .!= nothing],
                            NEG_LOGLIK[NEG_LOGLIK .!= nothing])...)
    D = try
        (DISTRIBUTIONS_DF[argmin(DISTRIBUTIONS_DF[:,3]), 2], DISTRIBUTIONS_DF[argmin(DISTRIBUTIONS_DF[:,3]), 1])
    catch
        (nothing, "Failure to fit into any of the distributions tested.")
    end
    return(D)
end

function ESTIMATE_LOD(vec_b::Vector{Float64})::Vector{Float64}
    D, D_name = BEST_FITTING_DISTRIBUTION(vec_b)
    println(string("Distribution used: ", D_name))
    if (D == nothing) | (std(vec_b) == 0)
        PVAL = repeat([1.0], inner=length(vec_b))
        LOD = PVAL .- 1.0
    else
        PVAL = [Distributions.cdf(D, x) <= Distributions.ccdf(D, x) ? 2*Distributions.cdf(D, x) : 2*Distributions.ccdf(D, x) for x in vec_b]
        LOD = -log.(10, PVAL)
    end
    return(LOD)
end

function GWAS_SIMPLREG(y::Vector{Float64}, X::Matrix{Float64})::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    n, m = size(X)
    k = 2 ### 1st df for intercept, and the second df for the allele
    vec_b = []
    vec_t = []
    vec_pval = []
    for i in 1:m
        # i = 1
        x = hcat(ones(n), X[:, i])
        b = try
                inv(x' * x) * x' * y
            catch
                pinv(x' * x) * x' * y
            end
        ε = y - (x * b)
        V_ε = (ε' * ε) / (n-k)
        se_ε = sqrt(V_ε)
        COV_b = V_ε * inv(x' * x)
        se_b = sqrt.([COV_b[1,1], COV_b[k,k]])
        t = b ./ se_b
        pval = Distributions.ccdf(Distributions.Chisq(k-1), t.^2)
        # pval = Distributions.cdf(Distributions.Normal(), -abs.(t)) + Distributions.ccdf(Distributions.Normal(), +abs.(t))
        append!(vec_b,    b[k])
        append!(vec_t,    t[k])
        append!(vec_pval, pval[k])
    end
    vec_lod = -log10.(vec_pval)
    # ### TEST PLOTS
    # LOD           = -log10.(vec_pval)
    # LOD_expected  = -log10.(cdf(Distributions.Uniform(0, 1), (collect(1:m) .- 0.5) ./ m))
    # LOD_threshold = -log10(0.05/m)
    # p1 = Plots.scatter(vec_b, legend=false)
    # p2 = Plots.scatter(LOD)
    # p3 = Plots.scatter(sort(LOD_expected), sort(LOD))
    # Plots.plot!([0,1], [0,1], seriestype=:straightline)
    # p4 = Plots.scatter(sort(LOD_expected)[1:(end-10)], sort(LOD)[1:(end-10)])
    # Plots.plot!([0,1], [0,1], seriestype=:straightline)
    # Plots.plot(p1, p2, p3, p4, layout=(2,2))
    return(vec_b, vec_t, vec_pval, vec_lod)
end

function GWAS_PLOT_MANHATTAN(vec_chr::Vector{String}, vec_pos::Vector{Int64}, vec_lod::Vector{Float64}, title::String="")::Plots.Plot{Plots.GRBackend}
    ### add lengths of all the chromosome to find xlim!!!!
    names_chr = sort(unique(vec_chr)) 
    l = 0
    for chr in names_chr
        # chr = sort(unique(vec_chr))[1]
        idx = vec_chr .== chr
        pos = vec_pos[idx]
        l += maximum(pos)
    end

    p1 = Plots.scatter([0], [0], xlims=[0, l], ylims=[minimum(vec_lod), maximum(vec_lod)], legend=false, markersize=0, markerstrokewidth=0, title=title)
    x0 = 0
    for chr in names_chr
        # chr = sort(unique(vec_chr))[1]
        idx = vec_chr .== chr
        pos = vec_pos[idx]
        x = x0 .+ pos
        y = vec_lod[idx]
        Plots.scatter!(x, y, legend=false,
                    markerstrokewidth=0.001,
                    markeralpha=0.4)
        x0 = maximum(x)
    end
    m = length(vec_lod)
    LOD_threshold = -log10(0.05/m)
    Plots.plot!([0,1], [LOD_threshold,LOD_threshold], seriestype=:straightline, legend=false)
    return(p1)
end

function GWAS_PLOT_QQ(vec_lod::Vector{Float64}, title::String="")::Plots.Plot{Plots.GRBackend}
    m = length(vec_lod)
    LOD_expected  = -log10.(cdf(Distributions.Uniform(0, 1), (collect(1:m) .- 0.5) ./ m))
    p2 = Plots.scatter(sort(LOD_expected), sort(vec_lod), markerstrokewidth=0.001, markeralpha=0.4, legend=false)
    Plots.plot!([0,1], [0,1], seriestype=:straightline, legend=false, title=title)
    return(p2)
end

function GWAS_PLOT(vec_chr::Vector{String}, vec_pos::Vector{Int64}, vec_lod::Vector{Float64})::Plots.Plot{Plots.GRBackend}
    p1 = GWAS_PLOT_MANHATTAN(vec_chr, vec_pos, vec_lod)
    p2 = GWAS_PLOT_QQ(vec_lod)
    p3 = Plots.plot(p1, p2, layout=(2,1))
    return(p3)
end

function GWAS_MULTIREG(y::Vector{Float64}, X::Matrix{Float64})::Vector{Float64}
    X = hcat(ones(size(X, 1)), Float64.(X))
    # b = X \ Y; ### Slow
    b = X' * inv(X * X') * y; ### Wayyyy faster! ~12 times faster!
    return(b[2:end])
end

### Proposed new method
function MODEL_SIMULATE_SIMPLREG(y::Vector{Float64}, n_degrees::Int64=3, n_sim::Int64=100, round_digits::Int64=4)::Vector{Float64}
    # n_degrees = 2
    # n_sim = 100
    n = length(y)
    k = 1 + n_degrees
    P = ones(n)
    for i in 1:n_degrees
        P = hcat(P, collect(1:n).^i)
    end
    b = inv(P' * P) * P' * y
    # e = y .- (P*b)
    # Ve = (e' * e) ./ (n-k)
    # Vb = Ve .* inv(P' * P)
    P_sim = ones(n_sim)
    for i in 1:n_degrees
        P_sim = hcat(P_sim, collect(range(1, n, length=100)).^i)
    end
    y_sim = round.(P_sim * b, digits=round_digits)
    return(y_sim)
end

function GWAS_SIMULREG(y::Vector{Float64}, X::Matrix{Float64})::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    _, m = size(X)
    y_sim = MODEL_SIMULATE_SIMPLREG(y)
    n = length(y_sim)
    k = 2
    vec_b = []
    vec_t = []
    vec_pval = []
    @time for j in 1:m
        # j = 37417
        # j = 103647
        x = X[:, j]
        x_sim = MODEL_SIMULATE_SIMPLREG(x)
        if var(x_sim) != 0.0
            X_sim = hcat(ones(length(x_sim)), x_sim)
            b = inv(X_sim' * X_sim) * X_sim' * y_sim
            e = y_sim .- (X_sim*b)
            Ve = (e' * e) ./ (n-k)
            Vb = Ve .* inv(X_sim' * X_sim)
            sb = sqrt.([Vb[1,1], Vb[k,1k]])
            t = b ./ sb
            pval = Distributions.ccdf(Distributions.Chisq(k-1), t.^2)
        else
            b = [0.0, 0.0]
            t = [0.0, 0.0]
            pval = [1.0, 1.0]
        end
        # p1 = Plots.plot(y, legend=false)
        # p2 = Plots.plot(x, legend=false)
        # p3 = Plots.plot(y_sim, legend=false)
        # p4 = Plots.plot(x_sim, legend=false)
        # Plots.plot(p1, p2, p3, p4, layout=(2,2))
        append!(vec_b, b[k])
        append!(vec_t, t[k])
        append!(vec_pval, pval[k])
    end
    vec_lod = -log10.(vec_pval)
    return(vec_b, vec_t, vec_pval, vec_lod)
end


end
