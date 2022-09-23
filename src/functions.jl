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
using Optim
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

function LOAD(phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_cols::Vector{Int}=[8,9], missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""])::Phenotype
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

function LOAD_OUT(tsv_gwalpha::String, include_all_sites::Bool=false)::Tuple{Vector{String}, Vector{Int64}, Vector{String}, Vector{Float64}, Vector{Float64}}
    # tsv_gwalpha = "/home/jeffersonfparil/Documents/poolgen/test/test_3-GWAlpha-maf_0.001.tsv"
    file = open(tsv_gwalpha, "r")
    vec_chr = []
    vec_pos = []
    vec_allele = []
    vec_freq = []
    vec_alpha = []
    while !eof(file)
        line = split(readline(file), "\t")
        if (include_all_sites == false) & (line[4] != "0")
            push!(vec_chr, line[1])
            push!(vec_pos, parse(Int, line[2]))
            push!(vec_allele, line[3])
            push!(vec_freq, parse(Float64, line[5]))
            push!(vec_alpha, parse(Float64, line[6]))
        end
    end
    close(file)
    return(vec_chr, vec_pos, vec_allele, vec_freq, vec_alpha)
end

function LOAD_PY(py_phenotype::String)::Tuple{String, Float64, Float64, Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}}
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

### SIMULATION (using Int64 for the genotype encoding so we can include multi-allelic loci in the future)
function BUILD_FOUNDING_GENOMES(n::Int64, m::Int64, l::Int64, k::Int64, ϵ::Int64=Int(1e+15), a::Int64=2, vec_chr_lengths=[], vec_chr_names=[])::Tuple{Vector{String}, Vector{Int64}, Vector{Int64}, Array{Int64, 3}}
    # n = 2                 ### number of founders
    # m = 10_000            ### number of loci
    # l = 2_300_000         ### total genome size
    # k = 7                 ### number of chromosomes
    # ϵ = Int(1e+15)        ### an arbitrarily large Int64 number to indicate no LD, i.e. the distance between the termini of 2 adjacent chromosomes is infinitely large LD-wise but we don;t want to use Inf as it is not Int64
    # a = 2                 ### number of alleles per locus (biallelic or a=2 by default)
    # vec_chr_lengths = []  ### optional vector of chromosome lengths
    # vec_chr_names = []    ### optional vector of chromosome names
    ### Define lengths of each chromosome
    if vec_chr_lengths == []
        vec_chr_lengths = repeat([Int(floor(l / k))], k)
        if sum(vec_chr_lengths)  < l
            vec_chr_lengths[end] = vec_chr_lengths[end] + (l - sum(vec_chr_lengths))
        end
    end
    ### Initiate the chromosome names output vector
    vec_chr = []
    ### Generate founders with two sets of chromosomes (randomly dsitributed 1's and 0's)
    B = Distributions.Binomial(a-1, 0.5)
    F = rand(B, 2, n, m)
    ### Sample SNP (or loci) locations from a uniform distrbution
    U = Distributions.Uniform(1, l)
    vec_pos = [0]
    while length(vec_pos) < m
        vec_pos = sort(unique(Int.(ceil.(rand(U, 2*m))))[1:m])
    end
    ### Calculate the preliminary distances between loci (we will correct for the distances between terminal loci of adjacent chromosomes below)
    vec_dist = append!([0], diff(vec_pos))
    chr_counter = 1
    for i in 1:m
        # i = 1000
        if vec_pos[i] > cumsum(vec_chr_lengths)[chr_counter]
             ### set the distance between the terminal loci of two adjacent chromosomes
             ### to some arbitrarily large number
            vec_dist[i] = ϵ
            ### Move to the next chromosome if we have not reached the final chromosome yet
            if chr_counter < length(vec_chr_lengths)
                chr_counter += 1
            end
        end
        ### Subtract the cummulative chromome length from the commulative position,
        ### so that the position in each chromosome starts at 1
        if chr_counter > 1
            vec_pos[i] -= cumsum(vec_chr_lengths)[chr_counter-1]
        end
        ### Populate the chromosome names output vector
        if vec_chr_names == []
            append!(vec_chr, [string(chr_counter)])
        else
            append!(vec_chr, [vec_chr_names[chr_counter]])
        end
    end
    ### Test plots where we should see k peaks each corresponding to a chromosome terminal
    # p1 = Plots.plot(vec_dist, title="distances", legend=false)
    # p2 = Plots.plot(vec_pos, title="positions", legend=false)
    # Plots.plot(p1, p2, layout=(1,2))
    return(vec_chr, vec_pos, vec_dist, F)
end

function MEIOSIS(X::Matrix{Int64}, distances::Vector{Int64}, dist_noLD::Int64=10_000)::Vector{Int64}
    # m = 1_000
    # X = sample([0, 1], (2, m))
    # ### Sample locations from a uniform distrbution    
    # U = Distributions.Uniform(1, 100*m)
    # positions = [0]
    # while length(positions) < m
    #     positions = sort(unique(Int.(ceil.(rand(U, 2*m))))[1:m])
    # end
    # distances = append!([0], diff(positions))
    # dist_noLD = 10_000
    ### Find the number of biallelic loci (X: 2 x m)
    _, m = size(X)
    ### Random sampling of first locus in  the chromosome
    vec_idx = zeros(Bool, m)
    vec_idx[1] = sample([true, false])
    ### Iterate across the remaining loci and determine if each locus belong to one or the other homologous chromosome
    for i in 2:m
        # i = 2
        p_LD = 1 - minimum([(1/2), (1/2) * (distances[i]/dist_noLD)]) ### probability of linkage (adjusted from 0.5 as a factor of the pairwise distance between loci)
        linked = rand(Distributions.Binomial(1, p_LD)) == 1
        if linked
            vec_idx[i] =  vec_idx[i-1]
        else
            vec_idx[i] = !vec_idx[i-1]
        end
    end
    Plots.plot(vec_idx) ### should show blocks of 0's and 1's
    ### Define the resulting allele
    vec_gamete = zeros(Int64, m)
    vec_gamete[  vec_idx] = X[1,   vec_idx] ### homologous chromosome 1
    vec_gamete[.!vec_idx] = X[2, .!vec_idx] ### homologous chromosome 2
    return(vec_gamete)
end

function SIMULATE(G::Array{Int64, 3}, vec_dist::Vector{Int64}, dist_noLD::Int64=10_000, o::Int64=10_000, t::Int64=1)::Array{Int64, 3}
    # n = 5                   ### number of founders
    # m = 10_000              ### number of loci
    # l = 20_000       ### total genome length
    # k = 1                   ### number of chromosomes
    # ϵ = Int(1e+15)          ### some arbitrarily large number to signify the distance at which LD is nil
    # a = 2                   ### number of alleles per locus
    # vec_chr_lengths = []    ### chromosome lengths
    # vec_chr_names = []      ### chromosome names 
    # vec_chr, vec_pos, vec_dist, G = BUILD_FOUNDING_GENOMES(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names)
    # dist_noLD = 10_000     ### distance at which LD is nil (related to ϵ)
    # o = 1_000               ### total number of simulated individuals
    # t = 10
    ### Simulate t genetations of random mating starting from the founders with a constant population size of o per generation
    _, n, m = size(G)
    P = zeros(Int64, 2, o, m)
    @showprogress for _ in 1:t
        for i in 1:o
            # i = 1
            g = G[:, sample(collect(1:n)), :]
            P[1, i, :] = MEIOSIS(g, vec_dist, dist_noLD) ### homologous chromosome 1
            P[2, i, :] = MEIOSIS(g, vec_dist, dist_noLD) ### homologous chromosome 2
        end
        G = copy(P)
        _, n, m  = size(G)
    end
    return(P)
end

function LD(P::Array{Int64, 3}, n_pairs::Int64=10_000)
    # n = 4                  ### number of founders
    # m = 10_000              ### number of loci
    # l = 135_000_000         ### total genome length
    # k = 5                   ### number of chromosomes
    # ϵ = Int(1e+15)          ### some arbitrarily large number to signify the distance at which LD is nil
    # a = 2                   ### number of alleles per locus
    # vec_chr_length = []    ### chromosome lengths
    # vec_chr_names = []      ### chromosome names
    # dist_noLD = 250_000     ### distance at which LD is nil (related to ϵ)
    # o = 1_000               ### total number of simulated individuals
    # t = 100                  ### number of o-sized random mating generations to simulate
    # vec_chr, vec_pos, vec_dist, F = BUILD_FOUNDING_GENOMES(n, m, l, k, ϵ, a, vec_chr_length, vec_chr_names) 
    # P = SIMULATE(F, vec_dist, dist_noLD, 10, 1)
    # P = SIMULATE(P, vec_dist, dist_noLD, 100, 1)
    # P = SIMULATE(P, vec_dist, dist_noLD, 1_000, 1)
    # P = SIMULATE(P, vec_dist, dist_noLD, 10_000, 1)
    # n_pairs = 2_000
    ### Calculate LD on the first chromosome
    _, n, m = size(P)
    mat_pairs = zeros(2, 1)
    while size(mat_pairs, 2) < n_pairs
        mat_pairs = Int.(ceil.(rand(2, 2*n_pairs) .* m))
        idx = abs.(mat_pairs[1,:] .- mat_pairs[2,:]) .> 0.00
        mat_pairs = mat_pairs[:, idx]
        mat_pairs = mat_pairs[:, 1:n_pairs]
    end
    vec_r2 = []
    vec_dist = []
    @showprogress for i in 1:n_pairs
        # i = 19
        idx = mat_pairs[:, i]
        P_2_adjacent_loci = P[:, :, idx]
        locus_1 = P_2_adjacent_loci[:, :, 1]
        alleles_1 = sort(unique(locus_1))
        P_alleles_1 = [sum(a .== locus_1) for a in alleles_1] ./ prod(size(locus_1))

        locus_2 = P_2_adjacent_loci[:, :, 2]
        alleles_2 = sort(unique(locus_2))
        P_alleles_2 = [sum(a .== locus_2) for a in alleles_2] ./ prod(size(locus_2))
        
        vec_2_loci_genotypes = [join(P_2_adjacent_loci[1, i, :]) for i in 1:o]
        append!(vec_2_loci_genotypes, [join(P_2_adjacent_loci[2, i, :]) for i in 1:o])
        genotypes = sort(unique(vec_2_loci_genotypes))
        P_2_loci_genotypes = [sum(g .== vec_2_loci_genotypes) for g in genotypes] ./ length(vec_2_loci_genotypes)

        vec_D = []
        for g in genotypes
            # g = genotypes[4]
            _alleles_ = parse.(Int64, String.(split(g, "")))
            
            P_AxB = P_alleles_1[alleles_1 .== _alleles_[1]][end] * P_alleles_2[alleles_2 .== _alleles_[2]][end]
            P_AB  = P_2_loci_genotypes[genotypes .== g][end]
            append!(vec_D, abs(P_AB - P_AxB))
        end
        D = mean(vec_D)
        append!(vec_r2, D^2 / (prod(P_alleles_1)*prod(P_alleles_2)))
        append!(vec_dist, abs(diff(vec_pos[idx])[end]))
    end
    Plots.scatter(vec_dist, vec_r2, legend=false, markerstrokewidth=0.001, markeralpha=0.4)

    using Polynomials
    f = Polynomials.fit(Float64.(vec_dist), Float64.(vec_r2), 2)
    Plots.plot!(f, 0, maximum(vec_dist), label="Fit")

end




### GWAS AND GP
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

function GWALPHA(syncx::String, py_phenotype::String, init::Int64, term::Int64, maf::Float64, penalty::Bool=true, out::String="")::String
    # syncx = "/home/jeffersonfparil/Documents/poolgen/test/test_3.syncx"
    # py_phenotype = "/home/jeffersonfparil/Documents/poolgen/test/test_3-pheno-any-filename.py"
    # file = open(syncx, "r")
    # seekend(file)
    # threads = 4
    # vec_positions = Int.(round.(collect(0:(position(file)/threads):position(file))))
    # close(file)
    # init = vec_positions[2] # init = 0
    # term = vec_positions[3] # file = open(syncx, "r"); seekend(file); term = position(file);  close(file)
    # maf = 0.001
    # penalty = true
    # out = ""
    ### Output syncx
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-GWAlpha.tsv")
    end
    file_out = open(out, "a")
    ### Extract phenotype information
    phenotype_name, sigma, min, max, perc, q, bins = LOAD_PY(py_phenotype)
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
    # vec_alpha = []
    # vec_chr = []
    # vec_pos = []
    # vec_allele = []
    # vec_1 = []
    # vec_pA = []
    vec_alleles = ["A", "T", "C", "G", "INS", "DEL", "N"]
    i = init
    while (i < term) & (!eof(file))
        locus = PARSE([PARSE(SyncxLine(i, readline(file)))])
        i = position(file)
        freqs = (locus.cou ./ sum(locus.cou, dims=1))'
        freqs[isnan.(freqs)] .= 0.0
	    allele_freqs = sum(freqs .* bins, dims=1)
        if (sum(locus.cou) != 0.0)
            if (minimum(allele_freqs[allele_freqs .!= 0.0]) >= maf) & (maximum(allele_freqs) < (1.0 - maf)) #locus filtering by mean MAF
                for allele in 1:7
                    # allele = 1
                    if (allele_freqs[allele] > 0.0) & (maximum(freqs[:,allele]) < 0.999999)  #allele filtering remove alleles with no counts and that the number of pools with allele frequency close to one should not occur even once!
                        freqA = freqs[:, allele]
                        pA = sum(freqA .* bins)
                        pB = 1 - pA
                        binA = (freqA .* bins) ./ pA
                        binB = ( (1 .- freqA) .* bins ) ./ (1-pA)
                        percA = cumsum(binA)
                        percB = cumsum(binB)
                        ### optimize (minimize) -log-likelihood of these major allele frequencies modelled as a beta distribution
                        # using Nelder-Mead optimization or Box minimisation (try-catch if one or the other fails with preference to Nelder-Mead)
                        lower_limits = [1e-20, 1e-20, 1e-20, 1e-20]
                        upper_limits = [1.0, 1.0, 1.0, 1.0]
                        initial_values = [0.1, 0.1, 0.1, 0.1]
                        b = try
                            Optim.optimize(beta->NLL_BETA(beta, percA, percB), initial_values, NelderMead())
                        catch
                            try
                                Optim.optimize(beta->NLL_BETA(beta, percA, percB), lower_limits, upper_limits, initial_values)
                            catch ### lower limits of 1e-20 to 1e-6 causes beta dist parameter values to shrink to zero somehow - so we're setting lower limits to 1e-5 instead
                                lower_limits = [1e-5, 1e-5, 1e-5, 1e-5]
                                Optim.optimize(beta->NLL_BETA(beta, percA, percB), lower_limits, upper_limits, initial_values)
                            end
                        end
                        muA = min + ((max-min)*b.minimizer[1]/(b.minimizer[1]+b.minimizer[2]))
                        muB = min + ((max-min)*b.minimizer[3]/(b.minimizer[3]+b.minimizer[4]))
                        ### compute alpha
                        if penalty
                            w_penalty = 2*sqrt(pA*pB)
                        else
                            w_penalty = 1.0
                        end
                        a = w_penalty*(muA - muB) / sigma
                        line = join([locus.chr[1], locus.pos[1], vec_alleles[allele], 1, pA, a], "\t")
                        # append!(vec_alpha, a)
                        # append!(vec_chr, locus.chr)
                        # append!(vec_pos, locus.pos)
                        # append!(vec_allele, allele)
                        # append!(vec_1, 1)
                        # append!(vec_pA, pA)
                    else
                        ## for explicit zero effects for null (low to none) frequency alleles
                        line = join([locus.chr[1], locus.pos[1], vec_alleles[allele], 0, allele_freqs[allele], 0.0], "\t")
                        # append!(vec_alpha, 0.0)
                        # append!(vec_chr, locus.chr)
                        # append!(vec_pos, locus.pos)
                        # append!(vec_allele, allele)
                        # append!(vec_1, 0)
                        # append!(vec_pA, allele_freqs[allele])
                    end
                    write(file_out, string(line, "\n"))
                end
            end
        end
    end
    close(file)
    close(file_out)
    return(out)
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
        pval = repeat([1.0], inner=length(vec_b))
        lod = pval .- 1.0
    else
        pval = [Distributions.cdf(D, x) <= Distributions.ccdf(D, x) ? 2*Distributions.cdf(D, x) : 2*Distributions.ccdf(D, x) for x in vec_b]
        lod = -log.(10, pval)
    end
    return(lod)
end

function ESTIMATE_LOD(vec_b::Vector{Float64}, D)::Vector{Float64}
    # D = [Distributions.Bernoulli, Distributions.Beta, Distributions.Binomial, Distributions.Categorical,
    #      Distributions.DiscreteUniform, Distributions.Exponential, Distributions.Normal, Distributions.Gamma,
    #      Distributions.Geometric, Distributions.Laplace, Distributions.Pareto, Distributions.Poisson,
    #      Distributions.InverseGaussian, Distributions.Uniform][7]
    distribution = Distributions.fit_mle(D, vec_b)
    pval = [Distributions.cdf(distribution, x) <= Distributions.ccdf(distribution, x) ? 2*Distributions.cdf(distribution, x) : 2*Distributions.ccdf(distribution, x) for x in vec_b]
    lod = -log.(10, pval)
    return(lod)
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

function GWAS_PLOT(tsv_gwalpha::String, estimate_empirical_lod::Bool)::Plots.Plot{Plots.GRBackend}
    # tsv_gwalpha = "/home/jeffersonfparil/Documents/poolgen/test/test_3-GWAlpha-maf_0.001.tsv"
    # estimate_empirical_lod = true
    vec_chr, vec_pos, vec_allele, vec_freq, vec_alpha = LOAD_OUT(tsv_gwalpha)
    if estimate_empirical_lod
        vec_lod = ESTIMATE_LOD(vec_alpha)
    end
    p1 = GWAS_PLOT_MANHATTAN(vec_chr, vec_pos, vec_lod)
    p2 = GWAS_PLOT_QQ(vec_lod)
    p3 = Plots.plot(p1, p2, layout=(2,1))
    return(p3)
end

function GWAS_SIMPLREG(syncx::String, init::Int64, term::Int64, maf::Float64, phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_col::Int=1, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], out::String="")::String
    # syncx = "/home/jeffersonfparil/Documents/poolgen/test/test_3.syncx"
    # file = open(syncx, "r")
    # seekend(file)
    # threads = 4
    # vec_positions = Int.(round.(collect(0:(position(file)/threads):position(file))))
    # close(file)
    # init = vec_positions[1] # init = 0
    # term = vec_positions[end] # file = open(syncx, "r"); seekend(file); term = position(file);  close(file)
    # maf = 0.001
    # phenotype = "/home/jeffersonfparil/Documents/poolgen/test/test_3-pheno-any-filename.csv"
    # delimiter = ","
    # header = true
    # id_col = 3 ### pool IDs
    # phenotype_col = 9 ### survival rate
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # out = ""
    ### Output syncx
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-GWASIMPLEREG.tsv")
    end
    file_out = open(out, "a")
    ### Load phenotype data and standard normalise so we don't need to fit an intercept for simplicity
    P = LOAD(phenotype, delimiter, header, id_col, [phenotype_col], missing_strings)
    y = P.phe[:,1]
    y = (y .- mean(y)) ./ std(y)
    n = length(y)
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

    vec_chr = []
    vec_pos = []
    vec_allele = []
    vec_b = []
    vec_σb = []
    vec_t = []
    vec_pval = []
    vec_alleles = ["A", "T", "C", "G", "INS", "DEL", "N"]
    i = init
    while (i < term) & (!eof(file))
        locus = PARSE([PARSE(SyncxLine(i, readline(file)))])
        i = position(file)
        X = locus.cou ./ sum(locus.cou, dims=1)
        X = X'
        if  (minimum(X[X .!= 0.0]) .>= maf) & (maximum(X[X .!= 0.0]) .<= (1 - maf))
            idx = collect(1:7)[(var(X, dims=1) .> 0.0)[1, :]]
            freqs = mean(X[:, idx], dims=1)[1,:]
            idx = idx[freqs .!= maximum(freqs)]
            X = X[:, idx]
            n, k = size(X)
            a = vec_alleles[idx]
            # Remove completey correlated alleles
            if k > 1
                C = LinearAlgebra.triu(cor(X), 1)
                idx = .!(maximum(C, dims=1) .≈ 1.0)[1, :]
                X = X[:, idx]
                a = a[idx]
            end
            n, k = size(X)
            for i in 1:k
                x = X[:,i]
                V =try
                    inv(x' * x)
                catch
                    pinv(x' * x)
                end
                b = V * X' * y
                ε = y - (X * b)
                Vε = (ε' * ε) / (n-k)
                if ((Vε < 0.0) | (Vε == Inf)) == false
                    Vb = Vε * V
                    Vb = Vb .- minimum(Vb)
                    σb = []
                    for j in 1:k
                        # j = 1
                        append!(σb, sqrt(Vb[j,j]))
                    end
                    # ### Correcting for the likelihood orrank correlation with small number of pools
                    # b = b .* (1 - (((n^2 - 1) / (2*n^2))^n))
                    t = b ./ σb
                    pval = Distributions.ccdf(Distributions.Chisq(k), t.^2)
                    append!(vec_chr, repeat(locus.chr, inner=length(a)))
                    append!(vec_pos, repeat(locus.pos, inner=length(a)))
                    append!(vec_allele, a)
                    append!(vec_b, b)
                    append!(vec_σb, σb)
                    append!(vec_t, t)
                    append!(vec_pval, pval)
                end
                # using HypothesisTests
                # for i in 1:k
                #     # i = 1
                #     # append!(vec_b, cor(y, X[:,i]))
                #     c = HypothesisTests.CorrelationTest(y, X[:,i])
                #     append!(vec_b, c.r)
                #     # append!(vec_b, (c.r)*(std(X[:,1])/std(y)))
                #     append!(vec_t, c.t)
                #     append!(vec_pval, HypothesisTests.pvalue(c))
                # end
            end
        end
    end
    
    Plots.histogram(vec_t)
    vec_TLOD = ESTIMATE_LOD(abs.(vec_t))
    GWAS_PLOT_MANHATTAN(String.(vec_chr), Int64.(vec_pos), vec_TLOD)
    GWAS_PLOT_QQ(vec_TLOD)

    vec_NLOD = ESTIMATE_LOD(Float64.(vec_b), Distributions.Normal)
    GWAS_PLOT_MANHATTAN(String.(vec_chr), Int64.(vec_pos), vec_NLOD)
    GWAS_PLOT_QQ(vec_NLOD)

    vec_ELOD = ESTIMATE_LOD((vec_b), Distributions.Exponential)
    GWAS_PLOT_MANHATTAN(String.(vec_chr), Int64.(vec_pos), vec_ELOD)
    GWAS_PLOT_QQ(vec_ELOD)

    vec_lod = ESTIMATE_LOD(Float64.(vec_b))
    # vec_lod = ESTIMATE_LOD(abs.(vec_b))
    # vec_lod = -log10.(Float64.(vec_pval))
    vec_lod[isnan.(vec_lod)] .= 0.0
    vec_lod[vec_lod .== Inf] .= maximum(vec_lod[vec_lod .!= Inf])
    vec_pval = 10 .^(-vec_lod)
    GWAS_PLOT_MANHATTAN(String.(vec_chr), Int64.(vec_pos), vec_lod)
    GWAS_PLOT_QQ(vec_lod)

    return(out)
end

function GWAS_MULTIREG(y::Vector{Float64}, X::Matrix{Float64})::Vector{Float64}
    X = hcat(ones(size(X, 1)), Float64.(X))
    # b = X \ Y; ### Slow
    b = X' * inv(X * X') * y; ### Wayyyy faster! ~12 times faster!
    return(b[2:end])
end

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
