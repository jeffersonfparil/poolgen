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
using StatsBase
using Dates
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

function SAVE(phenotype::Phenotype, filename::String, delimiter::String=",", header::Vector{String}=[""], missing_string::String="NA")
    # phenotype = Phenotype(["1", "2", "3", "4", "5"], ["y"], rand(5, 1))
    # filename = "test.csv"
    # delimiter = ","
    # header = ["id", "y"]
    # missing_string = "NA"
    file = open(filename, "a")
    if header != [""]
        write(file, string(join(header, ",")), "\n")
    end
    n, m = size(phenotype.phe)
    for i in 1:n
        phe = phenotype.phe[i, :]
        phe[ismissing.(phe)] .= missing_string
        write(file, string(join(vcat(phenotype.iid[i], phe), ",")), "\n")
    end
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

function LOAD_OUT(tsv_gwalpha::String, include_all_sites::Bool)::Tuple{Vector{String}, Vector{Int64}, Vector{String}, Vector{Float64}, Vector{Float64}}
    # tsv_gwalpha = "/home/jeffersonfparil/Documents/poolgen/test/test_3-GWAlpha-maf_0.001.tsv"
    # include_all_sites = false
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

function LOAD_OUT(tsv::String)::Tuple{Vector{String}, Vector{Int64}, Vector{String}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    # tsv = "/home/jeffersonfparil/Documents/poolgen/test/test_5-output.tsv"
    # @time vec_chr, vec_pos, vec_allele, vec_freq, vec_beta, vec_Vbeta, vec_pval = LOAD_OUT(tsv)
    file = open(tsv, "r")
    vec_chr = []
    vec_pos = []
    vec_allele = []
    vec_freq = []
    vec_beta = []
    vec_pval = []
    while !eof(file)
        line = split(readline(file), "\t")
        push!(vec_chr, line[1])
        push!(vec_pos, parse(Int, line[2]))
        push!(vec_allele, line[3])
        push!(vec_freq, parse(Float64, line[4]))
        push!(vec_beta, parse(Float64, line[5]))
        try
            push!(vec_pval, parse(Float64, line[6]))
        catch
            push!(vec_pval, 1.00)
        end
    end
    close(file)
    return(vec_chr, vec_pos, vec_allele, vec_freq, vec_beta, vec_pval)
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

function FILTER(χ::Window, maf::Float64=0.001, δ::Float64=1e-10, remove_insertions::Bool=false, remove_correlated_alleles::Bool=false, θ::Float64=0.99, centre::Bool=true)::Tuple{Matrix{Float64}, Vector{String}, Vector{Int64}, Vector{String}}
    # χ = LOAD("../test/test_3.syncx", false)
    # maf = 0.001
    # δ = 1e-10
    # remove_insertions = true
    # remove_correlated_alleles = true
    # θ = 0.99
    # centre = true
    X = Float64.(χ.cou')
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
    _X = reshape(sum(X, dims=1), 7, Int(size(X,2)/7))
    n, m = size(_X)
    vec_idx = []
    for j in 1:m
        x = _X[:, j]
        idx = x .!= 0
        idx_allele = collect(1:7)[idx]
        x = x[idx_allele]
        idx_allele = (j-1)*7 .+ idx_allele
        if sum(x .== minimum(x)) > 1
            saved_min_allele = idx_allele[x .== minimum(x)][2:end]
        else
            saved_min_allele = []
        end
        append!(vec_idx, append!(idx_allele, saved_min_allele))
    end
    X = X[:, vec_idx]
    ### Filter by minimum allele frequency
    idx = (minimum(X, dims=1)[1,:] .>= maf) .& (maximum(X, dims=1)[1,:] .<= 1-maf)
    X = X[:, idx]
    vec_idx = vec_idx[idx] ### update index
    ### Remove non-polymorphic loci
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
    vec_ale = repeat(["A", "T", "C", "G", "INS", "DEL", "N"], outer=m)[vec_idx]
    return(X, vec_chr, vec_pos, vec_ale)
end

### SIMULATION (using Int64 for the genotype encoding so we can include multi-allelic loci in the future)
function BUILD_FOUNDING_HETEROZYGOUS_GENOMES(n::Int64, m::Int64, l::Int64, k::Int64, ϵ::Int64=Int(1e+15), a::Int64=2, vec_chr_lengths::Vector{Int64}=([0]), vec_chr_names::Vector{String}=[""])::Tuple{Vector{String}, Vector{Int64}, Vector{Int64}, Array{Int64, 3}}
    # n = 2                 ### number of founders
    # m = 10_000            ### number of loci
    # l = 2_300_000         ### total genome size
    # k = 7                 ### number of chromosomes
    # ϵ = Int(1e+15)        ### an arbitrarily large Int64 number to indicate no LD, i.e. the distance between the termini of 2 adjacent chromosomes is infinitely large LD-wise but we don;t want to use Inf as it is not Int64
    # a = 2                 ### number of alleles per locus (biallelic or a=2 by default)
    # vec_chr_lengths = [0] ### optional vector of chromosome lengths
    # vec_chr_names = [""]  ### optional vector of chromosome names
    ### Define lengths of each chromosome
    if vec_chr_lengths == [0]
        vec_chr_lengths = repeat([Int(floor(l / k))], k)
        if sum(vec_chr_lengths)  < l
            vec_chr_lengths[end] = vec_chr_lengths[end] + (l - sum(vec_chr_lengths))
        end
    end
    ### Initiate the chromosome names output vector
    vec_chr = []
    ### Generate founders with two sets of chromosomes (randomly dsitributed 1's and 0's)
    B = Distributions.Binomial(a-1, 0.5)
    F = zeros(Int64, 2, n, m)
    F[1, :, :] = rand(B, 1, n, m)
    F[2, :, :] = abs.((a-1) .- F[1, :, :])
    ### Reset m if m > l
    m = minimum([m, l])
    ### Sample SNP (or loci) locations from a uniform distrbution
    vec_pos = sort(sample(collect(1:l), m, replace=false))
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
        if vec_chr_names == [""]
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

function SIMULATE_GENOMES(G::Array{Int64, 3}, vec_dist::Vector{Int64}, dist_noLD::Int64=10_000, o::Int64=10_000, t::Int64=1)::Array{Int64, 3}
    # n = 5                   ### number of founders
    # m = 10_000              ### number of loci
    # l = 20_000       ### total genome length
    # k = 1                   ### number of chromosomes
    # ϵ = Int(1e+15)          ### some arbitrarily large number to signify the distance at which LD is nil
    # a = 2                   ### number of alleles per locus
    # vec_chr_lengths = []    ### chromosome lengths
    # vec_chr_names = []      ### chromosome names 
    # vec_chr, vec_pos, vec_dist, G = BUILD_FOUNDING_HETEROZYGOUS_GENOMES(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names)
    # dist_noLD = 10_000     ### distance at which LD is nil (related to ϵ)
    # o = 1_000               ### total number of simulated individuals
    # t = 10
    ### Simulate t genetations of random mating starting from the founders with a constant population size of o per generation
    _, n, m = size(G)
    P = zeros(Int64, 2, o, m)
    @showprogress string("Simulating ", t, " generations:") for _ in 1:t
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

function LD(P::Array{Int64, 3}, vec_chr::Vector{String}, vec_pos::Vector{Int64}, chr::String="", window_size::Int64=1_000_000, n_pairs::Int64=10_000)::Tuple{Vector{Float64}, Vector{Int64}}
    # n = 10                  ### number of biallelic founders
    # m = 10_000              ### number of loci
    # l = 135_000_000         ### total genome length
    # k = 5                   ### number of chromosomes
    # ϵ = Int(1e+15)          ### some arbitrarily large number to signify the distance at which LD is nil
    # a = 4                   ### number of alleles per locus
    # vec_chr_length = []     ### chromosome lengths
    # vec_chr_names = []      ### chromosome names
    # dist_noLD = 250_000     ### distance at which LD is nil (related to ϵ)
    # o = 1_000               ### total number of simulated individuals
    # t = 10                  ### number of o-sized random mating generations to simulate
    # vec_chr, vec_pos, vec_dist, F = BUILD_FOUNDING_HETEROZYGOUS_GENOMES(n, m, l, k, ϵ, a, vec_chr_length, vec_chr_names) 
    # P = SIMULATE_GENOMES(F, vec_dist, dist_noLD, o, t)
    # chr = ""
    # window_size = 2*dist_noLD
    # n_pairs = 2_000
    ### Calculate LD on the first chromosome in a 2*dist_noLD window
    _, n, m = size(P)
    if chr == ""
        vec_idx = collect(1:m   )[vec_chr .== vec_chr[1]]
    else
        vec_idx = collect(1:m   )[vec_chr .== chr]
    end
    vec_idx = vec_idx[vec_pos[vec_idx] .<= window_size]
    min_idx = minimum(vec_idx) ### assumes correctly that the positions are sorted
    max_idx = maximum(vec_idx) ### assumes correctly that the positions are sorted
    mat_pairs = zeros(2, 1)
    while size(mat_pairs, 2) < n_pairs
        mat_pairs = Int.(ceil.((rand(2, 2*n_pairs) .* (max_idx - min_idx)) .+ min_idx))
        idx = abs.(mat_pairs[1,:] .- mat_pairs[2,:]) .> 0.00
        mat_pairs = mat_pairs[:, idx]
        mat_pairs = mat_pairs[:, 1:n_pairs]
    end
    vec_r2 = []
    vec_dist = []
    for i in 1:n_pairs
        # i = 19
        idx = mat_pairs[:, i]
        P_2_adjacent_loci = P[:, :, idx]
        locus_1 = P_2_adjacent_loci[:, :, 1]
        alleles_1 = sort(unique(locus_1))
        P_alleles_1 = [sum(a .== locus_1) for a in alleles_1] ./ prod(size(locus_1))

        locus_2 = P_2_adjacent_loci[:, :, 2]
        alleles_2 = sort(unique(locus_2))
        P_alleles_2 = [sum(a .== locus_2) for a in alleles_2] ./ prod(size(locus_2))
        
        vec_2_loci_genotypes = [join(P_2_adjacent_loci[1, i, :]) for i in 1:n]
        append!(vec_2_loci_genotypes, [join(P_2_adjacent_loci[2, i, :]) for i in 1:n])
        genotypes = []
        for a1 = alleles_1
            for a2 in alleles_2
                push!(genotypes, string(a1, a2))
            end
        end
        P_2_loci_genotypes = [sum(g .== vec_2_loci_genotypes) for g in genotypes] ./ length(vec_2_loci_genotypes)
        vec_D = []
        for g in genotypes
            # g = genotypes[4]
            _alleles_ = parse.(Int64, String.(split(g, "")))
            
            P_AxB = P_alleles_1[alleles_1 .== _alleles_[1]][end] * P_alleles_2[alleles_2 .== _alleles_[2]][end]
            P_AB  = P_2_loci_genotypes[genotypes .== g][end]
            append!(vec_D, abs(P_AB - P_AxB))
        end
        append!(vec_r2, sqrt(prod(vec_D)) / (prod(P_alleles_1)*prod(P_alleles_2)))
        append!(vec_dist, abs(diff(vec_pos[idx])[end]))
    end
    # Plots.scatter(vec_dist, vec_r2, legend=false, markerstrokewidth=0.001, markeralpha=0.4)
    # using Polynomials
    # f = Polynomials.fit(Float64.(vec_dist), Float64.(vec_r2), 10)
    # Plots.plot!(f, 0, maximum(vec_dist), label="Fit")
    return(vec_r2, vec_dist)
end

function SIMULATE(n::Int64, m::Int64, l::Int64, k::Int64, ϵ::Int64=Int(1e+15), a::Int64=2, vec_chr_lengths::Vector{Int64}=Int64.([0]), vec_chr_names::Vector{String}=[""], dist_noLD::Int64=10_000, o::Int64=1_000, t::Int64=10, nQTL::Int64=10, heritability::Float64=0.5, LD_chr::String="", LD_n_pairs::Int64=10_000, plot_LD::Bool=true)::Tuple{Vector{String}, Vector{Int64}, Matrix{Int64}, Vector{Float64}, Vector{Float64}}
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
    ### Instatiante founder genome/s
    vec_chr, vec_pos, vec_dist, G = BUILD_FOUNDING_HETEROZYGOUS_GENOMES(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names)
    ### Simulate random mating with constatnt population sizes for t genrations
    P = SIMULATE_GENOMES(G, vec_dist, dist_noLD, o, t) ### A 3-dimensional array of Int64 where the largest number, e.g.x, means that we have a macimum of x+1 alleles per locus
    _, n, m  = size(P)
    ### Assess LD
    if plot_LD
        LD_window_size = 2*dist_noLD
        vec_r2, vec_dist = LD(P, vec_chr, vec_pos, LD_chr, LD_window_size, LD_n_pairs)
        p = Plots.scatter(vec_dist, vec_r2, legend=false, xlab="Distance (bp)", ylab="r²")
        time_id = replace(join(collect(string(time()))[1:12], ""), "."=> "")
        Plots.savefig(p, string("Simulated_genome_LD-", time_id, ".png"))
    end
    ### Define genotype counts (fixed alleles/loci are all kept)
    vec_allele_counts_minus_one = maximum(P, dims=[1, 2])[1, 1, :]
    vec_alleles = collect(0:maximum(vec_allele_counts_minus_one)) ### e.g. [0, 1] for biallelic and [0, 1, 2] for triallelic
    m_new = m + sum(vec_allele_counts_minus_one[vec_allele_counts_minus_one .> 1] .- 1)
    X = zeros(Int64, n, m_new)
    m_new_counter = 1
    for i in 1:m
        # i = 1
        ### We count starting from the second allele, i.e. at index 2 for brevity since we only need n_alleles - 1 degrees of freedom anyway
        alleles = collect(2: vec_allele_counts_minus_one[i] + 1)
        ### Skip loci fixed for the first allele, i.e. "0" since we initiate X with zeros
        if alleles == []
            m_new_counter += 1
        else
            ### Else, iterate across n_alleles - 1 and count
            for j in alleles
                X[:, m_new_counter] = sum(P[:, :, i] .== vec_alleles[j], dims=1)[1, :]
                m_new_counter += 1
            end
        end
    end
    ### Update loci coordinates
    vec_chr_updated = []
    vec_pos_updated = []
    for i in 1:m
        append!(vec_chr_updated, repeat([vec_chr[i]], maximum([1, vec_allele_counts_minus_one[i]])))
        append!(vec_pos_updated, repeat([vec_pos[i]], maximum([1, vec_allele_counts_minus_one[i]])))
    end
    ### Simulate QTL effects
    QTL_idx = sort(Int.(ceil.(rand(nQTL) .* m_new)))
    QTL_effects = rand(Distributions.Normal(5, 1), nQTL)
    b = zeros(Float64, m_new)
    b[QTL_idx] = QTL_effects
    ### Simulate phenotypes
    g = X * b
    V_additive = var(g)
    V_residual = (V_additive / heritability) - V_additive
    e = rand(Distributions.Normal(0, sqrt(V_residual)), n)
    y = g + e
    return(vec_chr_updated, vec_pos_updated, X, y, b)
end

function POOL(X::Matrix{Int64}, y::Vector{Float64}, npools::Int64=5)::Tuple{Matrix{Float64}, Vector{Float64}}
    # n=5; m=100_000; l=135_000_000; k=5; ϵ=Int(1e+15); a=2; vec_chr_lengths=[]; vec_chr_names=[]; dist_noLD=500_000; o=100; t=10; nQTL=5; heritability=0.9
    # @time vec_chr, vec_pos, X, y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability)
    # npools = 5
    ### Sort the individuals into equally-sized npools pools
    idx = sortperm(y)
    y = y[idx]
    X = X[idx, :]
    n, m = size(X)
    vec_idx_pools = repeat([Int(floor(n/npools))], npools)
    add_to_last_pool = n - sum(vec_idx_pools)
    vec_idx_pools[end] = vec_idx_pools[end] + add_to_last_pool
    vec_idx_pools = append!([0], cumsum(vec_idx_pools))
    ### Generate the matrix of genotype frequencies per pool across m loci and vector of mean phenotypes per pool
    G = zeros(Float64, npools, m)
    p = zeros(Float64, npools)
    for i in 1:npools
        init = vec_idx_pools[i] + 1
        term = vec_idx_pools[i+1]
        G[i, :] = mean(X[init:term, :], dims=1) ./ 2
        p[i] = mean(y[init:term])
    end
    return(G, p)
end

function EXPORT_SIMULATED_DATA(vec_chr::Vector{String}, vec_pos::Vector{Int64}, X::Matrix{Int64}, y::Vector{Float64}, out_geno::String="", out_pheno::String="")::Tuple{String, String, String, String}
    # n=5; m=100_000; l=135_000_000; k=5; ϵ=Int(1e+15); a=2; vec_chr_lengths=[]; vec_chr_names=[]; dist_noLD=500_000; o=100; t=10; nQTL=5
    # @time vec_chr, vec_pos, X, y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL)
    # out_geno = ""
    # out_pheno = ""
    #####################
    ###               ###
    ### INDI-SEQ DATA ###
    ###               ###
    #####################
    ### Output syncx and csv files
    if (out_geno=="") | (out_pheno=="")
        id = replace(join(collect(string(time()))[1:12], ""), "."=> "")
    end
    if out_geno==""
        out_map = string.("Simulated-", id, ".map")
        out_bim = string.("Simulated-", id, ".bim")
        out_ped = string.("Simulated-", id, ".ped")
    else
        out_map = out_geno[1]
        out_bim = out_geno[2]
        out_ped = out_geno[3]
    end
    if out_pheno==""
        out_fam = string("Simulated-", id, ".fam")
    else
        out_fam = out_pheno
    end
    n, m = size(X)
    ################################
    ###@@@ Map and Bim files @@@####
    ################################
    map = open(out_map, "a")
    bim = open(out_bim, "a")
    mat_alleles = string.(zeros(Int64, m , 2))
    alleles = ["A", "T", "C", "G"]
    chromosomes = unique(vec_chr)
    chromosomes_idx = collect(1:length(chromosomes))
    for i in 1:m
        # i = 1
        chromosome_code = chromosomes_idx[chromosomes .== vec_chr[i]][1]
        variant_id = string(vec_chr[i], "-", vec_pos[i])
        position_in_morgans = "0"
        base_pair_coordinate = vec_pos[i]
        allele_1, allele_2 = sample(alleles, 2, replace=false)
        mat_alleles[i, :] = [allele_1, allele_2]
        map_line = string(join((chromosome_code,
                                variant_id,
                                position_in_morgans,
                                base_pair_coordinate), "\t"), "\n")
        bim_line = string(join((chromosome_code,
                                variant_id,
                                position_in_morgans,
                                base_pair_coordinate,
                                allele_1,
                                allele_2), "\t"), "\n")
        write(map, map_line)
        write(bim, bim_line)
    end
    close(map)
    close(bim)
    ###############################
    ###@@@ Fam and Ped files @@@###
    ###############################
    fam = open(out_fam, "a")
    ped = open(out_ped, "a")
    @showprogress for i in 1:n
        # i = 1
        family_id = i
        within_family_id = i
        within_family_id_father = 0
        within_family_id_mother = 0
        sex_code = 0
        phenotype_value = y[i]
        fam_line = string(join((family_id,
                                within_family_id,
                                within_family_id_father,
                                within_family_id_mother,
                                sex_code,
                                phenotype_value), "\t"), "\n")
        biallelic_geno = string.(zeros(Int64, 2, m))
        idx = X[i, :] .== 0
        biallelic_geno[:, idx] .= vcat(reshape(mat_alleles[idx, 1], (1, sum(idx))), reshape(mat_alleles[idx, 1], (1, sum(idx)))) ### homozygous for allele 1
        idx = X[i, :] .== 1
        biallelic_geno[:, idx] .= vcat(reshape(mat_alleles[idx, 1], (1, sum(idx))), reshape(mat_alleles[idx, 2], (1, sum(idx)))) ### heterozygote (allele 1)
        idx = X[i, :] .== 2
        biallelic_geno[:, idx] .= vcat(reshape(mat_alleles[idx, 2], (1, sum(idx))), reshape(mat_alleles[idx, 2], (1, sum(idx)))) ### homozygous for allele 2
        biallelic_geno = reshape(biallelic_geno, 2*m, :)
        ped_line = string(join((family_id,
                                within_family_id,
                                within_family_id_father,
                                within_family_id_mother,
                                sex_code,
                                phenotype_value,
                                join(biallelic_geno, "\t")), "\t"), "\n")
        write(fam, fam_line)
        write(ped, ped_line)
    end
    close(fam)
    close(ped)
    # ### test in plink and gemma
    # wget 'https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip'
    # wget 'https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz'
    # unzip plink_linux_x86_64_20220402.zip; rm plink_linux_x86_64_20220402.zip
    # gunzip gemma-0.98.5-linux-static-AMD64.gz; mv gemma-0.98.5-linux-static-AMD64 gemma
    # chmod +x plink
    # chmod +x gemma
    # time ./plink --file "Simulated-16644341691" --make-bed --out "Simulated-16644341691"
    # time ./gemma -bfile "Simulated-16644341691" -lm 1 -o "Simulated-16644341691"
    return(out_map, out_bim, out_ped, out_fam)
end

function EXPORT_SIMULATED_DATA(vec_chr::Vector{String}, vec_pos::Vector{Int64}, G::Matrix{Float64}, p::Vector{Float64}, out_geno::String="", out_pheno::String="")::Tuple{String, String}
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
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, X, y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 5
    # @time X, y = POOL(X, y, npools)
    # out_geno = ""
    # out_pheno = ""
    #####################
    ###               ###
    ### POOL-SEQ DATA ###
    ###               ###
    #####################
    ### Output syncx and csv files
    if (out_geno=="") | (out_pheno=="")
        id = replace(join(collect(string(time()))[1:12], ""), "."=> "")
    end
    if out_geno==""
        out_geno = string.("Simulated-", id, ".syncx")
    end
    if out_pheno==""
        out_pheno = string("Simulated-", id, ".csv")
    end
    n, m = size(G)
    ### Write simulated Pool-seq data into a syncx file
    geno = open(out_geno, "a")
    vec_unique_loci = unique(string.(vec_chr, "-", vec_pos))
    n_unique_loci = length(vec_unique_loci)
    ###@@@ Iterate per locus
    for i in 1:n_unique_loci
        # i = 1
        chr = split(vec_unique_loci[i], "-")[1]
        pos = parse(Int64, split(vec_unique_loci[i], "-")[2])
        idx = (vec_chr .== chr) .& (vec_pos .== pos)
        g_alt = G[:, idx]
        g_ref = 1 .- g_alt
        g = hcat(g_ref, g_alt)
        _, a = size(g)
        vec_counts_per_pool = repeat([""], n)
        ###@@@ Iterate per pool per locus
        for j in 1:n
            # j = 1
            vec_counts = zeros(Int64, 7)
            ###@@@ Iterate per allele per pool per locus
            d = Int(ceil(Distributions.rand(Distributions.Exponential(100), 1)[1])) ### simulate the total number of depth
            for k in 1:a
                # k = 1
                vec_counts[k] = Int(ceil(g[j, k] * d))
            end
            vec_counts_per_pool[j] = join(vec_counts, ":")
        end
        line = join(append!([chr, pos], vec_counts_per_pool), "\t")
        write(geno, string(line, "\n"))
    end
    close(geno)
    ### Write simulated phenotype data
    pheno = open(out_pheno, "a")
    header = string("Name,Phenotype\n")
    write(pheno, header)
    for i in 1:n
        # i = 1
        line = string("Pool-", i, ",", p[i], "\n")
        write(pheno, line)
    end
    close(pheno)

    return(out_geno, out_pheno)
end

### MODELS
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

function GENERATE_COVARIATE(syncx::String, nloci::Int64=1_000, θ::Float64=0.95, df::Int64=1, covariate::String=["XTX", "COR"][2])::Matrix{Float64}
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
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 5
    # @time X, y = POOL(_X, _y, npools)
    # @time syncx, csv = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X, y)
    # nloci = 1_000       ### number of equally-consecutively-spaced loci to include
    # θ = 0.95 ### maximum correlation between alleles/locus (Pearson's product moment correlation)
    # df = 1
    ### Extract allele frequencies
    ### Include nloci where the default is zero corresponding to all loci
    file = open(syncx, "r")
    m = countlines(file) ### total number of loci
    close(file)
    if (nloci == 0) | (nloci >= m)
        χ = LOAD(syncx, false)
    else
        _step_ = Int(round(m / nloci))
        file = open(syncx, "r")
        loci = []
        idx_lines = collect(1:_step_:m); counter = 1
        for i in 1:m
            if (i == idx_lines[counter])
                push!(loci, PARSE(SyncxLine(1, readline(file))))
                counter = minimum([nloci, counter + 1])
            else
                continue
            end
        end
        close(file)
        χ = PARSE(convert(Vector{LocusAlleleCounts}, loci))
    end
    ### Filter
    maf = 0.0001
    δ = 1e-10
    remove_insertions = true
    remove_correlated_alleles = true
    θ = 0.99
    centre = true
    X, vec_chr, vec_pos, vec_ref = FILTER(χ, maf, δ, remove_insertions, remove_correlated_alleles, θ, centre)
    n, p = size(X)
    ### Calculate covariate
    if covariate == "XTX"
        ### (1) X'*X/n covariate
        C = (X * X') ./ p
    elseif covariate == "COR"
        ### (2) Pearson's product-moment correlation
        C = cor(X')
    elseif covariate == "DIST"
        C = zeros(Int, n, n)
        for i in 1:n
            for j in 1:n
                C[i,j] = abs(vec_pos[i] - vec_pos[j])
            end
        end
    end
    ### Use the PCs if the we're asking for columns (i.e. df) less than the number of columns in C
    if (df < n)
        C = MultivariateStats.projection(MultivariateStats.fit(PCA, C; maxoutdim=df))
    end
    return(C)
end

function INVERSE(A::Array{T})::Matrix{Float64} where T <: Number
    try
        inv(A)
    catch
        pinv(A)
    end
end

function OLS(X::Array{T}, y::Array{T}, FE_method::String)::Vector{Float64} where T <: Number
    ### OLS using various methods and only outputting the effect estimates
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
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, X, y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # X = hcat(ones(o), X)
    # FE_method = ["CANONICAL", "MULTIALGORITHMIC", "N<<P"][3]
    # @time β̂ = OLS(X, y, FE_method)
    # p1 = Plots.scatter(b, title="True", legend=false);; p2 = Plots.scatter(abs.(β̂), title="Estimated", legend=false);; Plots.plot(p1, p2, layout=(2,1))
    if isa(X, Vector)
        X = reshape(X, (length(X), 1))
    end
    if FE_method == "CANONICAL"
        β̂ = INVERSE(X' * X) * (X' * y)
    elseif FE_method == "MULTIALGORITHMIC"
        β̂ = X \ y
    elseif FE_method == "N<<P"
        β̂ = X' * INVERSE(X * X') * y
    else
        println("Sorry. Invalid OLS method.")
        return(1)
    end
    return(β̂)
end

function OLS(X::Array{T}, y::Array{T})::Tuple{Vector{Float64}, Matrix{Float64}, Float64} where T <: Number
    ### Canonical OLS outputting estimates of the effects and variances
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
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, X, y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # X = hcat(ones(o), X)
    # @time β̂, Vβ̂, σ2ϵ̂ = OLS(X, y)
    # p1 = Plots.scatter(b, title="True", legend=false);; p2 = Plots.scatter(abs.(β̂), title="Estimated", legend=false);; Plots.plot(p1, p2, layout=(2,1))
    n, p = size(X)
    V = INVERSE(X' * X)
    β̂ = V * (X' * y)
    ε̂ = y - (X * β̂)
    σ2ϵ̂ = (ε̂' * ε̂) / (n-p)
    Vβ̂ = σ2ϵ̂ * V
    return(β̂, Vβ̂, σ2ϵ̂)
end

function GLMNET(X::Array{T}, y::Array{T}, alpha::Float64=1.0)::Tuple{Float64, Vector{T}} where T <: Number
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
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 5
    # @time X, y = POOL(_X, _y, npools)
    # X = hcat(ones(npools), X)
    # alpha = 1.0
    # β̂ = GLMNET(X, y, alpha)
    # p1 = Plots.scatter(b, title="True", legend=false);; p2 = Plots.scatter(abs.(β̂), title="Estimated", legend=false);; Plots.plot(p1, p2, layout=(2,1))
    ### Elastic net regularisation, where ridge if alpha=0.0, and LASSO if alpha=1.0
    n, _ = size(X)
    ### Remove inercept because we want the intercept to be unpenalised and to do that we let GLMNet.glmnetcv to automatically fit the intercept
    if X[:,1] == repeat([1], n)
        X = X[:,2:end]
    end
    ### GLMNet fit where if there isn't enough variation to fit the predictors which inhibits cross-validation to estimate mean loss, we resort to simply identifying the lambda which maximised the R2, i.e. dev_ratio
    β̂ = try
            cv = GLMNet.glmnetcv(X, y, intercept=true, alpha=alpha, tol=1e-7)
            idx = argmin(cv.meanloss)
            append!([cv.path.a0[idx]], cv.path.betas[:, idx])
        catch
            path = GLMNet.glmnet(X, y, intercept=true, alpha=alpha, tol=1e-7)
            idx = argmax(path.dev_ratio)
            append!([path.a0[idx]], path.betas[:, idx])
        end
    return(β̂[1], β̂[2:end])
end

function MM(X::Array{T}, y::Array{T}, Z::Array{T}, D::Array{T}, R::Array{T}, FE_method::String)::Tuple{Vector{Float64}, Vector{Float64}} where T <: Number
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
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 5
    # @time X, y = POOL(_X, _y, npools)
    # @time syncx, csv = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X, y)
    # @time K = GENERATE_COVARIATE(syncx, 1_000, 0.95, npools)
    # X = hcat(ones(npools), X)
    # ### GBLUP: y = Xβ + g + ϵ,
    # ###     where g = Zμ = μ,
    # ###         where Z = I(nxn), and (μ==g)~MVN(0, D),
    # ###             where D = σ2u * K
    # ###                 where K ≈ (X'X)/n
    # ###     and ϵ~MVN(0, R)
    # ###         where R = σ2e * I
    # ### SNP effects are fixed and we're controlling for random genotype effects.
    # Z = diagm(repeat([1.0], npools))    
    # σ2u = 2.0                           ### random effects variance
	# σ2e = 1.0                           ### error variance
    # D = σ2u * K                         ### random effects variance-covariance matrix
    # R = diagm(repeat([σ2e], npools))    ### homoscedastic error variance-covariance matrix
    # FE_method = "N<<P"
    # @time β̂, μ̂ = MM(X, y, Z, D, R, FE_method)
    # p1 = Plots.scatter(b, title="True", legend=false);; p2 = Plots.scatter(abs.(β̂), title="Estimated", legend=false);; Plots.plot(p1, p2, layout=(2,1))
    # ### RR-BLUP: y = Xβ + Zμ + ϵ,
    # ###     where Z are the SNPs,
    # ###     and μ~MVN(0, D),
    # ###         where D = σ2u * I
    # ###     and ϵ~MVN(0, R)
    # ###         where R = σ2e * I
    # ### Kinship effects are fixed and we're interested in random SNP effects with spherical variance.
    # Z = X[:, 2:end]                     ### SNPs with random effects
    # X = hcat(X[:,1], K)                 ### Intercept and kinship with fixed effects
    # σ2u = 2.0                           ### random effects variance
	# σ2e = 1.0                           ### error variance
    # D = diagm(repeat([σ2u], m))         ### homoscedastic random SNP effects variance-covariance matrix
    # R = diagm(repeat([σ2e], npools))    ### homoscedastic error variance-covariance matrix
    # FE_method = "N<<P"
    # @time β̂, μ̂ = MM(X, y, Z, D, R, FE_method)
    # p2 = Plots.scatter(b, title="True", legend=false);; p2 = Plots.scatter(abs.(μ̂), title="Estimated", legend=false);; Plots.plot(p1, p2, layout=(2,1))
    ### Variance-covariance matrix of y (combined inverse of R and D in Henderson's mixed model equations)
    V = (Z * D * Z') + R
	VI = INVERSE(V)
	### Mixed model equations
	###@@@ Fixed effects:
    if FE_method == "CANONICAL"
        ### More resource-intensive than X<<P
        β̂ = INVERSE(X' * VI * X) * (X' * VI * y)
    elseif FE_method == "N<<P"
        β̂ = (X' * VI) * INVERSE(X * X' * VI) * y
    else
        println("Sorry. Invalid OLS method.")
        return(1)
    end
	###@@@ Random effects:
    μ̂ = (D * Z' * VI) * (y - (X*β̂))
    return(β̂, μ̂)
end

function MM(X::Array{T}, y::Array{T}, Z::Array{T}, D::Array{T}, R::Array{T})::Tuple{Vector{Float64}, Vector{Float64}, Array{T}} where T <: Number
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
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 5
    # @time X, y = POOL(_X, _y, npools)
    # @time syncx, csv = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X, y)
    # @time K = GENERATE_COVARIATE(syncx, 1_000, 0.95, npools)
    # X = hcat(ones(npools), X[:, Int(round(rand()*size(X,2)))])
    # ### GBLUP: y = Xβ + g + ϵ,
    # ###     where g = Zμ = μ,
    # ###         where Z = I(nxn), and (μ==g)~MVN(0, D),
    # ###             where D = σ2u * K
    # ###                 where K ≈ (X'X)/n
    # ###     and ϵ~MVN(0, R)
    # ###         where R = σ2e * I
    # ### SNP effects are fixed and we're controlling for random genotype effects.
    # Z = diagm(repeat([1.0], npools))    
    # σ2u = 2.0                           ### random effects variance
	# σ2e = 1.0                           ### error variance
    # D = σ2u * K                         ### random effects variance-covariance matrix
    # R = diagm(repeat([σ2e], npools))    ### homoscedastic error variance-covariance matrix
    # @time β̂, μ̂, Σ̂ = MM(X, y, Z, D, R)
    # p1 = Plots.scatter(b, title="True", legend=false);; p2 = Plots.scatter(abs.(β̂), title="Estimated", legend=false);; Plots.plot(p1, p2, layout=(2,1))
    # ### RR-BLUP: y = Xβ + Zμ + ϵ,
    # ###     where Z are the SNPs,
    # ###     and μ~MVN(0, D),
    # ###         where D = σ2u * I
    # ###     and ϵ~MVN(0, R)
    # ###         where R = σ2e * I
    # ### Kinship effects are fixed and we're interested in random SNP effects with spherical variance.
    # Z = X[:, 2:end]                     ### SNPs with random effects
    # X = hcat(X[:,1], K)                 ### Intercept and kinship with fixed effects
    # σ2u = 2.0                           ### random effects variance
	# σ2e = 1.0                           ### error variance
    # D = diagm(repeat([σ2u], m))         ### homoscedastic random SNP effects variance-covariance matrix
    # R = diagm(repeat([σ2e], npools))    ### homoscedastic error variance-covariance matrix
    # _method_ = "N<<P"
    # @time β̂, μ̂, Σ̂ = MM(X, y, Z, D, R)
    # p2 = Plots.scatter(b, title="True", legend=false);; p2 = Plots.scatter(abs.(μ̂), title="Estimated", legend=false);; Plots.plot(p1, p2, layout=(2,1))
    ### Linear mixed model fitting which outputs the fixed effects variances for iterative regression
    ### Variance-covariance matrix of y (combined inverse of R and D in Henderson's mixed model equations)
    V = (Z * D * Z') + R
	VI = INVERSE(V)
    Σ̂ = INVERSE(X' * VI * X) ### Using canonical equations to estimate fixed effects so we estimate the fixed effect variances
	### Mixed model equations
    β̂ = Σ̂ * (X' * VI * y)
    μ̂ = (D * Z' * VI) * (y - (X*β̂))
    return(β̂, μ̂, Σ̂)
end

function NLL_MM(θ::Vector{T}, X::Array{T}, y::Array{T}, Z::Array{T}, K::Array{T}=zeros(2,2), FE_method::String=["CANONICAL", "N<<P"][2], method::String=["ML", "REML"][1])::Float64 where T <: Number
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
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 5
    # @time X, y = POOL(_X, _y, npools)
    # @time syncx, csv = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X, y)
    # K = zeros(2,2)
    # X = hcat(ones(npools), X)
    # ### GBLUP: y = Xβ + g + ϵ,
    # ###     where g = Zμ = μ,
    # ###         where Z = I(nxn), and (μ==g)~MVN(0, D),
    # ###             where D = σ2u * K
    # ###                 where K ≈ (X'X)/n
    # ###     and ϵ~MVN(0, R)
    # ###         where R = σ2e * I
    # ### SNP effects are fixed and we're controlling for random genotype effects.
    # Z = diagm(repeat([1.0], npools))
    # K = GENERATE_COVARIATE(syncx, 1_000, 0.95, npools)
    # FE_method = "N<<P"
    # method = "ML"
    # θ = [2.0, 1.0]
    # @time ll = NLL_MM(θ, X, y, Z, K, FE_method, method)
    # ### RR-BLUP: y = Xβ + Zμ + ϵ, (NOTE: CANNOT BE USED WITH REML SINCE Z IS NOT SQUARE!)
    # ###     where Z are the SNPs,
    # ###     and μ~MVN(0, D),
    # ###         where D = σ2u * I
    # ###     and ϵ~MVN(0, R)
    # ###         where R = σ2e * I
    # ### Kinship effects are fixed and we're interested in random SNP effects with spherical variance.
    # Z = X[:, 2:end]                     ### SNPs with random effects
    # X = hcat(X[:,1], K)                 ### Intercept and kinship with fixed effects
    # K = zeros(2, 2)
    # FE_method = "N<<P"
    # method = "ML"
    # θ = [2.0, 1.0]
    # @time ll = NLL_MM(θ, X, y, Z, K, FE_method, method)
    ### Homoscedastic error and random effect variances
	σ2u = θ[1]          # variance of the other random effects (assuming homoscedasticity)
    σ2e = θ[2]          # variance of the error effects (assuming homoscedasticity)
	n = length(y)		# number of individual samples
	l = size(X, 2)		# number of fixed effects
	nz, pz = size(Z)
    ### Random effects variance-covariance matrix (homoscedastic)
    if K == zeros(2, 2)
        D = diagm(repeat([σ2u], pz))
    else
        D = σ2u .* K
    end
    ### Error variance-covariance matrix (homoscedastic)
    R = diagm(repeat([σ2e], n))
    ### Variance-covariance matrix of y (combined inverse of R and D in Henderson's mixed model equations)
    V = (Z * D * Z') + R
    β̂, μ̂ = MM(X, y, Z, D, R, FE_method)
    ### Calculation negative log-likelihoods of variance θ
    if method == "ML"
        ### The negative log-likelihood function y given σ2e and σ2u
        μy = (X * β̂)
        neg_log_lik = 0.5 * ( log(abs(det(V))) + ((y - μy)' * INVERSE(V) * (y - μy)) + (n*log(2*pi)) )
    elseif method == "REML"
        if nz == pz
            ### NOTE: Z MUST BE SQUARE!
            M = try
                    INVERSE(LinearAlgebra.cholesky(σ2u*Z + R)) ### Cholesky decomposition
                catch
                    INVERSE(LinearAlgebra.lu(σ2u*Z + R).L) ### LU decomposition
                end
            y_new = M' * y
            intercept_new = sum(M', dims=2)
            V_new = M' * V * M
            n = length(y_new)
            ### Negative log-likelihood of σ2e and σ2u given X, y, and Z
            neg_log_lik = 0.5 * ( log(abs(det(V_new))) .+ ((y_new - intercept_new)' * INVERSE(V_new) * (y_new - intercept_new)) .+ (n*log(2*pi)) )[1,1]
        else
            if K == zeros(2, 2)
                println("Sorry, REML cannot be used with random SNP effects (e.g. RR-BLUP). Try GBLUP where the SNP effects are fixed.")
            else
                println("Sorry, REML is not a valid for a non-square Z matrix.")
            end
            return(2)
        end
    else
        println(string("Sorry. ", method, " is not a valid method of estimating the variances of the random effects effects. Please pick ML or REML."))
        return(1)
    end
    return(neg_log_lik)
end

function OPTIM_MM(X::Array{T}, y::Array{T}, Z::Array{T}, K::Array{T}=zeros(2,2), FE_method::String=["CANONICAL", "N<<P"][2], method::String=["ML", "REML"][1], inner_optimizer=[LBFGS(), BFGS(), SimulatedAnnealing(), GradientDescent(), NelderMead()][1], optim_trace::Bool=false)::Tuple{Float64, Float64} where T <: Number
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
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 5
    # @time X, y = POOL(_X, _y, npools)
    # @time syncx, csv = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X, y)
    # G = LOAD(syncx, false)
    # X, idx = FILTER(G)
    # vec_chr = repeat(G.chr, inner=7)[idx]
    # vec_pos = repeat(G.pos, inner=7)[idx]
    # ϕ = LOAD(csv, ",", true, 1, [2])
    # y = Float64.(ϕ.phe)[:,1]
    # @time C = GENERATE_COVARIATE(syncx, 1_000, 0.95, npools)
    # X = hcat(ones(npools), X)
    # ### GBLUP: y = Xβ + g + ϵ,
    # ###     where g = Zμ = μ,
    # ###         where Z = I(nxn), and (μ==g)~MVN(0, D),
    # ###             where D = σ2u * K
    # ###                 where K ≈ (X'X)/n
    # ###     and ϵ~MVN(0, R)
    # ###         where R = σ2e * I
    # ### SNP effects are fixed and we're controlling for random genotype effects.
    # Z = diagm(repeat([1.0], npools))    
    # K = C   
    # FE_method = "N<<P"
    # method = "ML"
    # optim_trace = true
    # inner_optimizer = LBFGS()
    # @time β̂, μ̂ = OPTIM_MM(X, y, Z, K, FE_method, method, inner_optimizer, optim_trace)
    # p1 = Plots.scatter(b, title="True", legend=false);; p2 = Plots.scatter(abs.(β̂[2:end]), title="Estimated", legend=false, markerstrokewidth=0.001, markeralpha=0.4);; Plots.plot(p1, p2, layout=(2,1))
    # ### RR-BLUP: y = Xβ + Zμ + ϵ,
    # ###     where Z are the SNPs,
    # ###     and μ~MVN(0, D),
    # ###         where D = σ2u * I
    # ###     and ϵ~MVN(0, R)
    # ###         where R = σ2e * I
    # ### Kinship effects are fixed and we're interested in random SNP effects with spherical variance.
    # Z = X[:, 2:end]                     ### SNPs with random effects
    # X = hcat(X[:,1], K)                 ### Intercept and kinship with fixed effects
    # K = zeros(2,2)
    # FE_method = "N<<P"
    # method = "ML"
    # inner_optimizer = LBFGS()
    # optim_trace = true
    # @time σ2u, σ2e = OPTIM_MM(X, y, Z, K, FE_method, method, inner_optimizer, optim_trace)
    # p2 = Plots.scatter(b, title="True", legend=false);; p2 = Plots.scatter(abs.(μ̂), title="Estimated", legend=false, markerstrokewidth=0.001, markeralpha=0.4);; Plots.plot(p1, p2, layout=(2,1))
    ### Find the optimum β̂ and μ̂
    lower_limits = [1e-5, 1e-5]
    upper_limits = [1e+5, 1e+5]
    initial_values = [1.0, 1.0]
    θ = Optim.optimize(parameters->NLL_MM(parameters, X, y, Z, K, FE_method, method),
                       lower_limits,
                       upper_limits,
                       initial_values,
                       Fminbox(inner_optimizer),
                       Optim.Options(f_tol = 1e-20,
                                     g_tol = 1e-10,
                                     iterations = 1_000,
                                     store_trace = false,
                                     show_trace = optim_trace,
                                     show_every=1, 
                                     time_limit=NaN))
    σ2u = θ.minimizer[1] # variance of the other random effects (assuming homoscedasticity)
    σ2e = θ.minimizer[2] # variance of the error effects (assuming homoscedasticity)
    # ### Output messages
    @show θ
    if (θ.f_converged) | (θ.g_converged)
        println("CONVERGED! 😄")
    else
        println("DID NOT CONVERGE! 😭")
    end
    return(σ2u, σ2e)
end

function BOOTSTRAP(ρ, x::Array{T}, y::Array{T}, F::Function, F_params="")::Tuple{Int64, Int64} where T <: Number
    # n = 5
    # x = rand(n)
    # # x = rand(n, 1)
    # y = rand(n)
    # F = cor
    # F_params = ""
    # # F = OLS
    # # F_params = "CANONICAL"
    # ρ = F(x, y, F_params)
    n = length(x)
    niter = minimum([100, Int64(n*(n-1)/2)])
    vec_ρ = []
    for i in 1:niter
        x_rand = sample(x, n; replace=false)
        y_rand = sample(y, n; replace=false)
        if F_params != ""
            if isa(F_params, Vector)
                append!(vec_ρ, F(x_rand, y_rand, F_params...))
            else
                append!(vec_ρ, F(x_rand, y_rand, F_params))
            end
        else
            append!(vec_ρ, F(x_rand, y_rand))
        end
    end
    positive_cases = sum(abs.(vec_ρ) .>= abs(ρ))
    return(positive_cases, niter)
end

function BOOTSTRAP_PVAL(ρ, x::Array{T}, y::Array{T}, F::Function, F_params="", nburnin::Int64=10_000, δ::Float64=1e-10, maxiter::Int64=1_000)::Vector{Float64} where T <: Number
    # n = 5
    # x = rand(n)
    # y = x * 10
    # F = cor
    # ρ = F(x, y)
    # nburnin = 10_000
    # δ = 1e-10
    # maxiter = 1_000
    pval = [0.50]
    positive_cases = 0
    total_cases = 0
    ### Burn-in steps
    for i in 1:nburnin
        _p, _t = BOOTSTRAP(ρ, x, y, F, F_params)
        positive_cases += _p
        total_cases += _t
        append!(pval, positive_cases / total_cases)
    end
    ### Bootstrap until convergence or maximum number of iterations is reached
    iter = 0
    pval_prev, pval_curr = pval[(end-1):end]
    while (mean(diff(pval[(end-10):end])) >= δ) | (iter >= maxiter)
        iter += 1
        _p, _t = BOOTSTRAP(ρ, x, y, F, F_params)
        positive_cases += _p
        total_cases += _t
        append!(pval, positive_cases / total_cases)
        pval_prev, pval_curr = pval[(end-1):end]
    end
    ### Output
    return(pval)
end

function OLS_ITERATIVE(syncx::String, init::Int64, term::Int64, maf::Float64, phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_col::Int=1, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], covariate::String=["", "XTX", "COR"][2], covariate_df::Int64=1, out::String="")::String
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
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 5
    # @time X, y = POOL(_X, _y, npools)
    # @time syncx, phenotype = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X, y)
    # file = open(syncx, "r")
    # seekend(file)
    # threads = 1
    # vec_positions = Int.(round.(collect(0:(position(file)/threads):position(file))))
    # close(file)
    # init = vec_positions[1] # init = 0
    # term = vec_positions[end] # file = open(syncx, "r"); seekend(file); term = position(file);  close(file)
    # maf = 0.001
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # covariate = "COR"
    # covariate_df = 1
    # out = ""
    # tsv = OLS_ITERATIVE(syncx, init, term, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, covariate, covariate_df, out)
    # vec_chr, vec_pos, vec_allele, vec_freq, vec_beta, vec_pval = LOAD_OUT(tsv)
    # p1 = Plots.scatter(b, title="True", legend=false);; p2 = Plots.scatter(-log10.(vec_pval), title="Estimated", legend=false, markerstrokewidth=0.001, markeralpha=0.4);; Plots.plot(p1, p2, layout=(2,1))
    # # syncx = "/home/jeffersonfparil/Documents/poolgen/test/test_3.syncx"
    # # file = open(syncx, "r")
    # # seekend(file)
    # # threads = 4
    # # vec_positions = Int.(round.(collect(0:(position(file)/threads):position(file))))
    # # close(file)
    # # init = vec_positions[1] # init = 0
    # # term = vec_positions[end] # file = open(syncx, "r"); seekend(file); term = position(file);  close(file)
    # # maf = 0.001
    # # phenotype = "/home/jeffersonfparil/Documents/poolgen/test/test_3-pheno-any-filename.csv"
    # # delimiter = ","
    # # header = true
    # # id_col = 3 ### pool IDs
    # # phenotype_col = 9 ### survival rate
    # # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # # covariate = "COR"
    # # covariate_df = 1
    # # out = ""
    ### Output tab-delimeted file including the (1) chromosome name, (2) position, (3) allele, (4) allele frequency, (5) allele effect, an (6) p-value
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-OLS_ITERATIVE.tsv")
    end
    file_out = open(out, "a")
    ### Load phenotype data and standard normalise so we don't need to fit an intercept for simplicity
    ϕ = LOAD(phenotype, delimiter, header, id_col, [phenotype_col], missing_strings)
    y = ϕ.phe[:,1]
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
    ### Prepare covariate
    if covariate != ""
        ρ = GENERATE_COVARIATE(syncx, 1_000, 0.95, covariate_df, covariate)
    end
    ### Regress
    # vec_chr = []
    # vec_pos = []
    # vec_allele = []
    # vec_β̂ = []
    # vec_σβ̂ = []
    # vec_t = []
    # vec_pval = []
    vec_alleles = ["A", "T", "C", "G", "INS", "DEL", "N"]
    pos = init
    while (pos < term) & (!eof(file))
        # @show position(file)
        locus = PARSE([PARSE(SyncxLine(pos, readline(file)))])
        pos = position(file)
        X = locus.cou ./ sum(locus.cou, dims=1)
        X = X'
        ### Continue if allele frequences are larger that the minimum allele frequence (and vice-versa for 1-maf); and if the locus is polymorphic
        if (minimum(X[X .!= 0.0]) .>= maf) & (maximum(X[X .!= 0.0]) .<= (1 - maf)) & (sum(var(X, dims=1)) > 0.0)
            ### Keep p-1 alleles where p is the number of polymorphic alleles
            idx = collect(1:7)[(var(X, dims=1) .> 0.0)[1, :]]
            freqs = mean(X[:, idx], dims=1)[1,:]
            idx = idx[freqs .!= maximum(freqs)]
            X = X[:, idx]
            n, k = size(X)
            a = vec_alleles[idx]
            # Remove completey correlated alleles
            if k > 1
                A = LinearAlgebra.triu(cor(X), 1)
                idx = .!(maximum(A, dims=1) .≈ 1.0)[1, :]
                X = X[:, idx]
                a = a[idx]
            end
            n, k = size(X)
            for i in 1:k
                x = reshape(X[:,i], (n, 1))
                p = 1 ### the number of parameters to estimate which is equal to 1 since we're not including an intercept as we centered y
                ### Prepend the covariate
                if covariate != ""
                    x = hcat(ρ, x)
                    _, p = size(x)
                end
                ### OLS using the canonical method and outputting the estimates of the effects and variances
                β̂, Vβ̂, σ2ϵ̂ = OLS(x, y)
                ### Test if we have reasonable residual variance estimate then proceed to output
                if ((σϵ̂ < 0.0) | (σϵ̂ == Inf)) == false
                    σβ̂ = []
                    for j in 1:p
                        # j = 1
                        append!(σβ̂, sqrt(Vβ̂[j,j]))
                    end
                    t = β̂ ./ σβ̂
                    pval = Distributions.ccdf(Distributions.Chisq(k), t.^2)
                    line = join([locus.chr[1],
                                 locus.pos[1],
                                 a[i],
                                 mean(x[:,end]),
                                 β̂[end],
                                 pval[end]], "\t")
                    write(file_out, string(line, "\n"))
                    # append!(vec_chr, locus.chr[1])
                    # append!(vec_pos, locus.pos[1])
                    # append!(vec_allele, a)
                    # append!(vec_β̂, β̂[end])
                    # append!(vec_σβ̂, σβ̂[end])
                    # append!(vec_t, t[end])
                    # append!(vec_pval, pval[end])
                end
            end
        end
    end
    close(file)
    close(file_out)
    return(out)
end

function BOO_ITERATIVE(syncx::String, init::Int64, term::Int64, maf::Float64, phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_col::Int=1, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], F::Function=OLS, F_params="", nburnin::Int64=1_000, δ::Float64=1e-10, maxiter::Int64=1_000, out::String="")::String
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
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time _vec_chr, _vec_pos, _X, _y, _β = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 5
    # @time X, y = POOL(_X, _y, npools)
    # @time syncx, phenotype = EXPORT_SIMULATED_DATA(_vec_chr, _vec_pos, X, y)
    # file = open(syncx, "r")
    # seekend(file)
    # threads = 4
    # vec_positions = Int.(round.(collect(0:(position(file)/threads):position(file))))
    # close(file)
    # init = vec_positions[1] # init = 0
    # term = vec_positions[end] # file = open(syncx, "r"); seekend(file); term = position(file);  close(file)
    # maf = 0.001
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # # F = cor
    # F = OLS
    # # F_params=""
    # F_params="CANONICAL"
    # out = ""
    # nburnin = 10_000
    # δ = 1e-10
    # maxiter = 1_000
    # @time tsv = BOO_ITERATIVE(syncx, init, term, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, F, F_params, nburnin, δ, maxiter, out)
    # vec_chr, vec_pos, vec_allele, vec_freq, vec_beta, vec_pval = LOAD_OUT(tsv)
    # p1 = Plots.scatter(_β, title="True", legend=false);; p2 = Plots.scatter(-log10.(vec_pval), title="Estimated", legend=false, markerstrokewidth=0.001, markeralpha=0.4);; Plots.plot(p1, p2, layout=(2,1))
    ### Output tab-delimeted file including the (1) chromosome name, (2) position, (3) allele, (4) allele frequency, (5) allele effect, an (6) p-value
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-BOO_ITERATIVE.tsv")
    end
    file_out = open(out, "a")
    ### Load phenotype data and standard normalise so we don't need to fit an intercept for simplicity
    ϕ = LOAD(phenotype, delimiter, header, id_col, [phenotype_col], missing_strings)
    y = ϕ.phe[:,1]
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
    ### Correlate
    # vec_chr = []
    # vec_pos = []
    # vec_allele = []
    # vec_ρ = []
    # vec_pval = []
    vec_alleles = ["A", "T", "C", "G", "INS", "DEL", "N"]
    pos = init
    while (pos < term) & (!eof(file))
        # @show position(file)
        locus = PARSE([PARSE(SyncxLine(pos, readline(file)))])
        pos = position(file)
        X = locus.cou ./ sum(locus.cou, dims=1)
        X = X'
        ### Continue if allele frequences are larger that the minimum allele frequence (and vice-versa for 1-maf); and if the locus is polymorphic
        if (minimum(X[X .!= 0.0]) .>= maf) & (maximum(X[X .!= 0.0]) .<= (1 - maf)) & (sum(var(X, dims=1)) > 0.0)
            ### Keep p-1 alleles where p is the number of polymorphic alleles
            idx = collect(1:7)[(var(X, dims=1) .> 0.0)[1, :]]
            freqs = mean(X[:, idx], dims=1)[1,:]
            idx = idx[freqs .!= maximum(freqs)]
            X = X[:, idx]
            n, k = size(X)
            a = vec_alleles[idx]
            # Remove completey correlated alleles
            if k > 1
                A = LinearAlgebra.triu(cor(X), 1)
                idx = .!(maximum(A, dims=1) .≈ 1.0)[1, :]
                X = X[:, idx]
                a = a[idx]
            end
            n, k = size(X)
            for i in 1:k
                # i = 1
                x = X[:,i:i]
                if F_params != ""
                    if isa(F_params, Vector)
                        ρ = F(x, y, F_params...)[1]
                    else
                        ρ = F(x, y, F_params)[1]
                    end
                else
                    ρ = F(x, y)[1]
                end
                pval = BOOTSTRAP_PVAL(ρ, x, y, F, F_params, nburnin, δ, maxiter)[end]
                line = join([locus.chr[1],
                                 locus.pos[1],
                                 a[i],
                                 mean(x[:,end]),
                                 ρ,
                                 pval], "\t")
                write(file_out, string(line, "\n"))
                # append!(vec_chr, locus.chr[1])
                # append!(vec_pos, locus.pos[1])
                # append!(vec_allele, a)
                # append!(vec_ρ, ρ)
                # append!(vec_pval, pval)
            end
        end
    end
    close(file)
    close(file_out)
    return(out)
end

function LMM_ITERATIVE(syncx::String, init::Int64, term::Int64, maf::Float64, phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_col::Int=1, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], covariate::String=["", "XTX", "COR"][2], model::String=["GBLUP", "RRBLUP"][1], method::String=["ML", "REML"][1], FE_method::String=["CANONICAL", "N<<P"][2], inner_optimizer=[LBFGS(), BFGS(), SimulatedAnnealing(), GradientDescent(), NelderMead()][1], optim_trace::Bool=false, out::String="")::String
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
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 5
    # @time X, y = POOL(_X, _y, npools)
    # @time syncx, phenotype = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X, y)
    # file = open(syncx, "r")
    # seekend(file)
    # threads = 1
    # vec_positions = Int.(round.(collect(0:(position(file)/threads):position(file))))
    # close(file)
    # init = vec_positions[1] # init = 0
    # term = vec_positions[end] # file = open(syncx, "r"); seekend(file); term = position(file);  close(file)
    # maf = 0.001
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # covariate = "COR"
    # model = ["GBLUP", "RRBLUP"][1]
    # method = ["ML", "REML"][1]
    # FE_method = ["CANONICAL", "N<<P"][2]
    # inner_optimizer=[LBFGS(), BFGS(), SimulatedAnnealing(), GradientDescent(), NelderMead()][1]
    # optim_trace = true
    # out = ""
    # @time tsv = LMM_ITERATIVE(syncx, init, term, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, covariate, model, method, FE_method, inner_optimizer, optim_trace, out)
    # vec_chr, vec_pos, vec_allele, vec_freq, vec_beta, vec_pval = LOAD_OUT(tsv)
    # p1 = Plots.scatter(b, title="True", legend=false);; p2 = Plots.scatter(-log10.(vec_pval), title="Estimated", legend=false, markerstrokewidth=0.001, markeralpha=0.4);; Plots.plot(p1, p2, layout=(2,1))
    ### Output tab-delimeted file including the (1) chromosome name, (2) position, (3) allele, (4) allele frequency, (5) allele effect, an (6) p-value
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-LMM_ITERATIVE.tsv")
    end
    file_out = open(out, "a")
    ### Load phenotype data and standard normalise so we don't need to fit an intercept for simplicity
    ϕ = LOAD(phenotype, delimiter, header, id_col, [phenotype_col], missing_strings)
    y = ϕ.phe[:,1]
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
    ### Prepare covariate and variance-covariance matrix
    C = GENERATE_COVARIATE(syncx, 1_000, 0.95, n, "COR")
    if model == "GBLUP"
        Z = diagm(repeat([1.0], n))    ### Genotypes have random effects...
        K = C                               ### ... and are distributed normally μ=0, and Σ=σ2g*K
    elseif model == "RRBLUP"
        K = zeros(2,2)                      ### Kinships have random effects and are spherically distributed proportional to σ2g
    else
        println(string("Sorry ", model, " is not implemented."))
        println("Please choose from 'GBLUP' and 'RRBLUP'.")
    end
    ### Regress
    # vec_chr = []
    # vec_pos = []
    # vec_allele = []
    # vec_β̂ = []
    # vec_σβ̂ = []
    # vec_t = []
    # vec_pval = []
    vec_alleles = ["A", "T", "C", "G", "INS", "DEL", "N"]
    pos = init
    while (pos < term) & (!eof(file))
        # @show position(file)
        locus = PARSE([PARSE(SyncxLine(pos, readline(file)))])
        pos = position(file)
        X = locus.cou ./ sum(locus.cou, dims=1)
        X = X'
        ### Continue if allele frequences are larger that the minimum allele frequence (and vice-versa for 1-maf); and if the locus is polymorphic
        if (minimum(X[X .!= 0.0]) .>= maf) & (maximum(X[X .!= 0.0]) .<= (1 - maf)) & (sum(var(X, dims=1)) > 0.0)
            ### Keep p-1 alleles where p is the number of polymorphic alleles
            idx = collect(1:7)[(var(X, dims=1) .> 0.0)[1, :]]
            freqs = mean(X[:, idx], dims=1)[1,:]
            idx = idx[freqs .!= maximum(freqs)]
            X = X[:, idx]
            n, k = size(X)
            a = vec_alleles[idx]
            # Remove completey correlated alleles
            if k > 1
                A = LinearAlgebra.triu(cor(X), 1)
                idx = .!(maximum(A, dims=1) .≈ 1.0)[1, :]
                X = X[:, idx]
                a = a[idx]
            end
            n, k = size(X)
            for i in 1:k
                # i = 1
                x = reshape(X[:,i], (n, 1))
                p = 1 ### the number of parameters to estimate which is equal to 1 since we're not including an intercept as we centered y
                ### Add the intercept and covariate
                if model == "GBLUP"
                    _X_ = hcat(ones(n), x)  ### SNPs have fixed effects
                elseif model == "RRBLUP"
                    Z = x                   ### SNPs have random effects
                    _X_ = hcat(ones(n), C)  ### Kinships have effects
                end
                _, p = size(_X_)
                ### Linear mixed model fitting using the canonical method and outputting the estimates of the effects and variances
                σ2u, σ2e = OPTIM_MM(_X_, y, Z, K, FE_method, method, inner_optimizer, optim_trace)
                ### Random effects variance-covariance matrix (homoscedastic)
                if model == "GBLUP"
                    D = σ2u .* K
                elseif model == "RRBLUP"
                    D = diagm(repeat([σ2u], size(Z,2)))
                end
                ### Error variance-covariance matrix (homoscedastic)
                R = diagm(repeat([σ2e], n))
                ### Solve the mixed model equations
                β̂, μ̂, Σ̂ = MM(_X_, y, Z, D, R)
                if model == "GBLUP"
                    b = β̂[end]
                    W = b^2 / Σ̂[end, end] ### Wald's test statistic
                elseif model == "RRBLUP"
                    b = μ̂[end]
                    W = b^2 / D[end, end] ### Wald's test statistic
                end
                pval = Distributions.ccdf(Distributions.Chisq(1), W)
                ### Output
                line = join([locus.chr[1],
                             locus.pos[1],
                             a[i],
                             mean(x[:,end]),
                             b,
                             pval], "\t")
                write(file_out, string(line, "\n"))
                # append!(vec_chr, locus.chr[1])
                # append!(vec_pos, locus.pos[1])
                # append!(vec_allele, a)
                # append!(vec_β̂, β̂[end])
                # append!(vec_σβ̂, σβ̂[end])
                # append!(vec_t, t[end])
                # append!(vec_pval, pval[end])
            end
        end
    end
    close(file)
    close(file_out)
    return(out)
end

function OLS_MULTIVAR(syncx::String, maf::Float64, phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_col::Int=1, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], FE_method::String=["CANONICAL", "N<<P"][2], out::String="")::String
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
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 5
    # @time X, y = POOL(_X, _y, npools)
    # @time syncx, phenotype = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X, y)
    # maf = 0.001
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # FE_method = ["CANONICAL", "N<<P"][2]
    # out = ""
    # @time tsv = OLS_MULTIVAR(syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, FE_method, out)
    ### Output tab-delimeted file including the (1) chromosome name, (2) position, (3) allele, (4) allele frequency, (5) allele effect, an (6) p-value
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-OLS_MUTIVAR.tsv")
    end
    ### Load genotype data
    χ = LOAD(syncx, false)
    ### Filter
    δ = 1e-10
    remove_insertions = true
    remove_correlated_alleles = true
    θ = 0.99
    centre = false
    X, vec_chr, vec_pos, vec_ale = FILTER(χ, maf, δ, remove_insertions, remove_correlated_alleles, θ, centre)
    vec_frq = mean(X, dims=1)
    n, p = size(X)
    ### Load phenotype data
    ϕ = LOAD(phenotype, delimiter, header, id_col, [phenotype_col], missing_strings)
    y = ϕ.phe[:,1]
    y = (y .- mean(y)) ./ std(y)
    ### Check is we have the same number of individuals in the genotype and phenotype data
    @assert n == length(y) "Genotype and phenotype data mismatch!"
    ### Fit (Note: not intercept as we normalised the y)
    β̂ = OLS(X, y, FE_method)
    ### Output
    file_out = open(out, "a")
    for i in 1:p
        line = join([vec_chr[i],
                     vec_pos[i],
                     vec_ale[i],
                     vec_frq[i],
                     β̂[i],
                     "NA"], "\t")
        write(file_out, string(line, "\n"))
    end
    close(file_out)
    return(out)
end

function ELA_MULTIVAR(syncx::String, maf::Float64, phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_col::Int=1, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], alpha::Float64=1.0, out::String="")::String
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
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 5
    # @time X, y = POOL(_X, _y, npools)
    # @time syncx, phenotype = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X, y)
    # maf = 0.001
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # alpha = 1.0
    # out = ""
    # @time tsv = ELA_MULTIVAR(syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, alpha, out)
    ### Output tab-delimeted file including the (1) chromosome name, (2) position, (3) allele, (4) allele frequency, (5) allele effect, an (6) p-value
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-ELA_MUTIVAR_alpha_", round(alpha, digits=2),".tsv")
    end
    ### Load genotype data
    χ = LOAD(syncx, false)
    ### Filter
    δ = 1e-10
    remove_insertions = true
    remove_correlated_alleles = true
    θ = 0.99
    centre = false
    X, vec_chr, vec_pos, vec_ale = FILTER(χ, maf, δ, remove_insertions, remove_correlated_alleles, θ, centre)
    vec_frq = mean(X, dims=1)
    n, p = size(X)
    ### Load phenotype data
    ϕ = LOAD(phenotype, delimiter, header, id_col, [phenotype_col], missing_strings)
    y = ϕ.phe[:,1]
    y = (y .- mean(y)) ./ std(y)
    ### Check is we have the same number of individuals in the genotype and phenotype data
    @assert n == length(y) "Genotype and phenotype data mismatch!"
    ### Fit (Note: not intercept as we normalised the y)
    a0, β̂ = GLMNET(X, y, alpha)
    ### Output
    file_out = open(out, "a")
    for i in 1:p
        line = join([vec_chr[i],
                     vec_pos[i],
                     vec_ale[i],
                     vec_frq[i],
                     β̂[i],
                     "NA"], "\t")
        write(file_out, string(line, "\n"))
    end
    close(file_out)
    return(out)
end

function LMM_MULTIVAR(syncx::String, maf::Float64, phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_col::Int=1, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], covariate::String=["", "XTX", "COR"][2], model::String=["GBLUP", "RRBLUP"][1], method::String=["ML", "REML"][1], FE_method::String=["CANONICAL", "N<<P"][2], inner_optimizer=[LBFGS(), BFGS(), SimulatedAnnealing(), GradientDescent(), NelderMead()][1], optim_trace::Bool=false, out::String="")::String
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
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 5
    # @time X, y = POOL(_X, _y, npools)
    # @time syncx, phenotype = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X, y)
    # maf = 0.001
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # covariate = "COR"
    # model = "GBLUP"
    # model = "RRBLUP"
    # method = "ML"
    # FE_method = "N<<P"
    # inner_optimizer = LBFGS()
    # optim_trace = true
    # out = ""
    # @time tsv = LMM_MULTIVAR(syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, covariate, model, method, FE_method, inner_optimizer, optim_trace, out)
    # vec_chr, vec_pos, vec_allele, vec_freq, vec_beta, vec_pval = LOAD_OUT(tsv)
    ### Output tab-delimeted file including the (1) chromosome name, (2) position, (3) allele, (4) allele frequency, (5) allele effect, an (6) p-value
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-LMM_MUTIVAR_", method,".tsv")
    end
    ### Load genotype data
    χ = LOAD(syncx, false)
    ### Filter
    δ = 1e-10
    remove_insertions = true
    remove_correlated_alleles = true
    θ = 0.99
    centre = false
    X, vec_chr, vec_pos, vec_ale = FILTER(χ, maf, δ, remove_insertions, remove_correlated_alleles, θ, centre)
    vec_frq = mean(X, dims=1)
    n, p = size(X)
    ### Load phenotype data
    ϕ = LOAD(phenotype, delimiter, header, id_col, [phenotype_col], missing_strings)
    y = ϕ.phe[:,1]
    y = (y .- mean(y)) ./ std(y)
    ### Check is we have the same number of individuals in the genotype and phenotype data
    @assert n == length(y) "Genotype and phenotype data mismatch!"
    ### Prepare the design matrices for the fixed and random effects; and variance-covariance matrix of the random effects
    if model == "GBLUP"
        Z = diagm(repeat([1.0], n))    ### Genotypes have random effects...
        K = GENERATE_COVARIATE(syncx, 1_000, 0.95, n, "COR") ### ... and are distributed normally μ=0, and Σ=σ2g*K
    elseif model == "RRBLUP"
        Z = X
        if covariate == ""
            X = ones(n, 1)
        else
            X = GENERATE_COVARIATE(syncx, 1_000, 0.95, n, "COR") ### ... and are distributed normally μ=0, and Σ=σ2g*K
        end
        K = zeros(2,2)                 ### Kinships have random effects and are spherically distributed proportional to σ2g
    else
        println(string("Sorry ", model, " is not implemented."))
        println("Please choose from 'GBLUP' and 'RRBLUP'.")
    end
    ### Linear mixed model fitting using the canonical method and outputting the estimates of the effects and variances
    σ2u, σ2e = OPTIM_MM(X, y, Z, K, FE_method, method, inner_optimizer, optim_trace)
    ### Random effects variance-covariance matrix (homoscedastic)
    if model == "GBLUP"
        D = σ2u .* K
    elseif model == "RRBLUP"
        D = diagm(repeat([σ2u], size(Z,2)))
    end
    ### Error variance-covariance matrix (homoscedastic)
    R = diagm(repeat([σ2e], n))
    ### Solve the mixed model equations
    β̂, μ̂ = MM(X, y, Z, D, R, FE_method)
    if model == "GBLUP"
        b = β̂
    elseif model == "RRBLUP"
        b = μ̂
    end
    ### Output
    file_out = open(out, "a")
    for i in 1:p
        line = join([vec_chr[i],
                     vec_pos[i],
                     vec_ale[i],
                     vec_frq[i],
                     b[i],
                     "NA"], "\t")
        write(file_out, string(line, "\n"))
    end
    close(file_out)
    return(out)
end

### CROSS-VALIDATION
function PREDICT(tsv::String, syncx_validate::String)::Vector{Float64}
    # n = 5                 ### number of founders
    # m = 10_000            ### number of loci
    # l = 135_000_000       ### total genome length
    # k = 5                 ### number of chromosomes
    # ϵ = Int(1e+15)        ### some arbitrarily large number to signify the distance at which LD is nil
    # a = 2                 ### number of alleles per locus
    # vec_chr_lengths = [0] ### chromosome lengths
    # vec_chr_names = [""]  ### chromosome names 
    # dist_noLD = 500_000   ### distance at which LD is nil (related to ϵ)
    # o = 1_000             ### total number of simulated individuals
    # t = 10                ### number of random mating constant population size generation to simulate
    # nQTL = 10             ### number of QTL to simulate
    # heritability = 0.5    ### narrow(broad)-sense heritability as only additive effects are simulated
    # LD_chr = ""           ### chromosome to calculate LD decay from
    # LD_n_pairs = 10_000   ### number of randomly sampled pairs of loci to calculate LD
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 50
    # @time X, y = POOL(_X, _y, npools)
    # idx = collect(1:npools)
    # idx_training = sort(sample(idx, Int(round(npools*0.75)), replace=false))
    # idx_validate = idx[[sum(idx_training .== x)==0 for x in idx]]
    # @time syncx_training, phenotype_training = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X[idx_training, :], y[idx_training])
    # @time syncx_validate, phenotype_validate = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X[idx_validate, :], y[idx_validate])
    # maf = 0.001
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # FE_method = ["CANONICAL", "N<<P"][2]
    # alpha = 1.0
    # covariate = "COR"; model = "GBLUP"
    # # covariate = "";    model = "RRBLUP"
    # method = "ML"
    # FE_method = "N<<P"
    # inner_optimizer = LBFGS()
    # optim_trace = true
    # out = ""
    # # @time tsv = OLS_MULTIVAR(syncx_training, maf, phenotype_training, delimiter, header, id_col, phenotype_col, missing_strings, FE_method, out)
    # # @time tsv = ELA_MULTIVAR(syncx_training, maf, phenotype_training, delimiter, header, id_col, phenotype_col, missing_strings, alpha, out)
    # @time tsv = LMM_MULTIVAR(syncx_training, maf, phenotype_training, delimiter, header, id_col, phenotype_col, missing_strings, covariate, model, method, FE_method, inner_optimizer, optim_trace, out)
    # @time ŷ = PREDICT(tsv, syncx_validate)
    # cor(ŷ, y[idx_validate])
    # Plots.scatter(ŷ, y[idx_validate])
    ### Load allele effects
    vec_chr, vec_pos, vec_ale, vec_freq, vec_beta, _ = LOAD_OUT(tsv)
    ### Genotype
    χ = LOAD(syncx_validate, false)
    p = length(χ.chr)
    ### Consolidate alleles and loci across trained effects and validation genotype data
    idx = []
    vec_alleles = ["A", "T", "C", "G", "INS", "DEL", "N"]
    for i in 1:p
        # i = 1
        idx_chr = χ.chr[i] .== vec_chr
        idx_pos = χ.pos[i] .== vec_pos
        if (sum(idx_chr)>0) & (sum(idx_pos)>0)
            idx_ale = vec_ale[idx_chr .& idx_pos] .== vec_alleles
            append!(idx, collect(((7*(i-1))+1):(7*i))[idx_ale])
        end            
    end
    ### It is required that we capture all the alleles and loci in the estimated effects file (tsv) are present in the validation genotype data (syncx_validate)
    if length(idx) < length(vec_beta)
        println("ERROR: Missing loci in prediction dataset.")
        println("Please, consider havng the same loci covered in both training and prediction sets, i.e.")
        println(string("in '", tsv, "' and '", syncx_validate, "' files."))
        return(1)
    end
    X = χ.cou'[:, idx]
    ŷ = X * vec_beta
    return(ŷ)
end

function CV_METRICS(y::Vector{T}, ŷ::Vector{T}, y_training::Vector{T})::Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Plots.Plot{Plots.GRBackend}} where T <: Number
    # n = 5                 ### number of founders
    # m = 10_000            ### number of loci
    # l = 135_000_000       ### total genome length
    # k = 5                 ### number of chromosomes
    # ϵ = Int(1e+15)        ### some arbitrarily large number to signify the distance at which LD is nil
    # a = 2                 ### number of alleles per locus
    # vec_chr_lengths = [0] ### chromosome lengths
    # vec_chr_names = [""]  ### chromosome names 
    # dist_noLD = 500_000   ### distance at which LD is nil (related to ϵ)
    # o = 1_000             ### total number of simulated individuals
    # t = 10                ### number of random mating constant population size generation to simulate
    # nQTL = 10             ### number of QTL to simulate
    # heritability = 0.5    ### narrow(broad)-sense heritability as only additive effects are simulated
    # LD_chr = ""           ### chromosome to calculate LD decay from
    # LD_n_pairs = 10_000   ### number of randomly sampled pairs of loci to calculate LD
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 50
    # @time X, y = POOL(_X, _y, npools)
    # idx = collect(1:npools)
    # idx_training = sort(sample(idx, Int(round(npools*0.75)), replace=false))
    # idx_validate = idx[[sum(idx_training .== x)==0 for x in idx]]
    # @time syncx_training, phenotype_training = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X[idx_training, :], y[idx_training])
    # @time syncx_validate, phenotype_validate = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X[idx_validate, :], y[idx_validate])
    # maf = 0.001
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # FE_method = ["CANONICAL", "N<<P"][2]
    # alpha = 1.0
    # covariate = "COR"; model = "GBLUP"
    # # covariate = "";    model = "RRBLUP"
    # method = "ML"
    # FE_method = "N<<P"
    # inner_optimizer = LBFGS()
    # optim_trace = true
    # out = ""
    # @time tsv = OLS_MULTIVAR(syncx_training, maf, phenotype_training, delimiter, header, id_col, phenotype_col, missing_strings, FE_method, out)
    # # @time tsv = ELA_MULTIVAR(syncx_training, maf, phenotype_training, delimiter, header, id_col, phenotype_col, missing_strings, alpha, out)
    # # @time tsv = LMM_MULTIVAR(syncx_training, maf, phenotype_training, delimiter, header, id_col, phenotype_col, missing_strings, covariate, model, method, FE_method, inner_optimizer, optim_trace, out)
    # @time ŷ = PREDICT(tsv, syncx_validate)
    # ϕ_validate = LOAD(phenotype_validate, delimiter, header, id_col, [phenotype_col], missing_strings)
    # y = Float64.(ϕ_validate.phe[:,1])
    # ϕ_training = LOAD(phenotype_training, delimiter, header, id_col, [phenotype_col], missing_strings)
    # y_training = Float64.(ϕ_training.phe[:,1])
    # correlation_pearson, correlation_spearman, correlation_kendall, R2, R2_adj, MAE, MBE, RAE, MSE, RMSE, RRMSE, RMSLE, p = CV_METRICS(y, ŷ, y_training)
    ### Transform the predictions into the range of the training data
    ŷ = (ŷ .* std(y_training)) .+ mean(y_training)
    ### Correlation metrics
    correlation_pearson = cor(y, ŷ)
    correlation_spearman = corspearman(y, ŷ)
    correlation_kendall = corkendall(y, ŷ)
    ### Coefficients of determination
    n = length(y)
    p = 2 ### intercept and y as predictors
    _β̂, _Vβ̂, σ2ϵ̂ = OLS(hcat(ones(n), y), ŷ)
    R2 = 1 - (σ2ϵ̂ / var(ŷ; corrected=false))
    R2_adj = 1 - (1 - R2) * ((n-1) / (n-p))
    ### Mean absolute error
    MAE = mean(abs.(y .- ŷ))
    ### Mean bias error
    MBE = mean(y .- ŷ)
    ### Relative absolute error
    RAE = sum(abs.(y .- ŷ)) / sum(abs.(y .- mean(y)))
    # RAE = (MAE * n) / sum(abs.(y .- mean(y)))
    ### Mean square error
    MSE = mean((y .- ŷ).^2)
    ### Root mean square error
    RMSE = sqrt(MSE)
    ### Relative root mean square error
    RRMSE = sqrt(MSE / sum(ŷ.^2))
    ### Root mean squared logarithmic error
    min_y = [minimum(y)<0 ? abs(minimum(y)) : 0][1] ### to get rid of negative values for 
    min_ŷ = [minimum(ŷ)<0 ? abs(minimum(ŷ)) : 0][1] ### to get rid of negative values for 
    RMSLE = sqrt(mean((log.(y .+ min_y .+ 1) .- log10.(ŷ .+ min_ŷ .+ 1)).^2))
    ### Plot
    x = vcat(y, ŷ)
    δ = std(x) / sqrt(length(x))
    limits = [minimum(x)-δ, maximum(x)+δ]
    p = Plots.scatter(y, ŷ, xlims=limits, ylims=limits, legend=false, markerstrokewidth=0.001, markeralpha=0.4, title="")
    Plots.plot!(p, [0,1], [0 ,1], seriestype=:straightline, linecolor=:gray, legend=false)
    ### Output metrics and plot
    return(correlation_pearson, correlation_spearman, correlation_kendall,
           R2, R2_adj,
           MAE, MBE, RAE,
           MSE, RMSE, RRMSE, RMSLE,
           p)
end

function CV_OLS_MULTIVAR(nfold::Int64, nrep::Int64, syncx::String, maf::Float64, phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_col::Int=1, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], FE_method::String=["CANONICAL", "N<<P"][2], out::String="")::String
    # n = 5                 ### number of founders
    # m = 10_000            ### number of loci
    # l = 135_000_000       ### total genome length
    # k = 5                 ### number of chromosomes
    # ϵ = Int(1e+15)        ### some arbitrarily large number to signify the distance at which LD is nil
    # a = 2                 ### number of alleles per locus
    # vec_chr_lengths = [0] ### chromosome lengths
    # vec_chr_names = [""]  ### chromosome names 
    # dist_noLD = 500_000   ### distance at which LD is nil (related to ϵ)
    # o = 1_000               ### total number of simulated individuals
    # t = 10                ### number of random mating constant population size generation to simulate
    # nQTL = 10             ### number of QTL to simulate
    # heritability = 0.5    ### narrow(broad)-sense heritability as only additive effects are simulated
    # LD_chr = ""           ### chromosome to calculate LD decay from
    # LD_n_pairs = 10_000   ### number of randomly sampled pairs of loci to calculate LD
    # plot_LD = false       ### simulate# plot simulated LD decay
    # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 50
    # @time X, y = POOL(_X, _y, npools)
    # nfold = 10
    # nrep = 3
    # @time syncx, phenotype = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X, y)
    # maf = 0.001
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # FE_method = ["CANONICAL", "N<<P"][2]
    # out = ""
    # cv_tsv = CV_OLS_MULTIVAR(nfold, nrep, syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, FE_method, out)
    # ###########################
    # nfold = 10
    # nrep = 1
    # syncx = "Simulated-16663168544.syncx"
    # maf = 0.0001
    # phenotype = "Simulated-16663168544.csv"
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # FE_method = "N<<P"
    # out = ""
    # poolgen.user_functions.functions.CV_OLS_MULTIVAR(nfold, nrep, syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, FE_method, out)
    # ###########################
    ### Output tab-delimeted file including the (1) chromosome name, (2) position, (3) allele, (4) allele frequency, (5) allele effect, an (6) p-value
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-OLS_MULTIVAR_CV.tsv")
    end
    ### Load genotype data
    χ = LOAD(syncx, true)
    p, n = size(χ.cou)
    ### Load phenotype data
    ϕ = LOAD(phenotype, delimiter, header, id_col, [phenotype_col], missing_strings)
    y = Float64.(ϕ.phe[:,1])
    ### Check is we have the same number of individuals in the genotype and phenotype data
    @assert n == length(y) "Genotype and phenotype data mismatch!"
    ### Cross-validate for nrep randomisations and k-fold cross-validation
    vec_rep = []
    vec_fold = []
    vec_correlation_pearson = []
    vec_correlation_spearman = []
    vec_correlation_kendall = []
    vec_R2 = []
    vec_R2_adj = []
    vec_MAE = []
    vec_MBE = []
    vec_RAE = []
    vec_MSE = []
    vec_RMSE = []
    vec_RRMSE = []
    vec_RMSLE = []
    vec_p = []
    vec_idx = collect(1:n)
    for i in 1:nrep
        # i = 1
        println("############################################")
        println(string("Time: ", Dates.format(now(), "Y-U-d"), " | UTC: ", Dates.format(now(), "HH:MM:SS")))
        println(string("Replication ", i, " of ", nrep))
        ### Randomise observations
        vec_idx = sample(vec_idx, n, replace=false)
        ### Divide the observations into nfold partitions
        vec_fld = repeat(collect(1:nfold), inner=Int(floor(n/nfold)))
        if length(vec_fld) < n
            append!(vec_fld, repeat([nfold], n-length(vec_fld)))
        end
        # @showprogress string(nfold, "-fold cross-validation: ") for j in 1:nfold
        #     # j= 1
        #     idx_training = vec_fld .!= j
        #     idx_validate = vec_fld .== j
        #     syncx_training, phenotype_training = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X[idx_training, :], y[idx_training], string(syncx, "-CV-rep_", i, "-fold_", j, "-TRAINING.syncx"), string(syncx, "-CV-rep_", i, "-fold_", j, "-TRAINING.csv"))
        #     syncx_validate, phenotype_validate = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X[idx_validate, :], y[idx_validate], string(syncx, "-CV-rep_", i, "-fold_", j, "-VALIDATE.syncx"), string(syncx, "-CV-rep_", i, "-fold_", j, "-VALIDATE.csv"))
        #     tsv = OLS_MULTIVAR(syncx_training, maf, phenotype_training, delimiter, header, id_col, phenotype_col, missing_strings, FE_method, out)
        #     ŷ = PREDICT(tsv, syncx_validate)
        #     correlation_pearson, correlation_spearman, correlation_kendall, R2, R2_adj, MAE, MBE, RAE, MSE, RMSE, RRMSE, RMSLE, p = CV_METRICS(y[idx_validate], ŷ, y[idx_training])
        #     append!(vec_rep, rep)
        #     append!(vec_fold, fold)
        #     append!(vec_correlation_pearson, correlation_pearson)
        #     append!(vec_correlation_spearman, correlation_spearman)
        #     append!(vec_correlation_kendall, correlation_kendall)
        #     append!(vec_R2, R2)
        #     append!(vec_R2_adj, R2_adj)
        #     append!(vec_MAE, MAE)
        #     append!(vec_MBE, MBE)
        #     append!(vec_RAE, RAE)
        #     append!(vec_MSE, MSE)
        #     append!(vec_RMSE, RMSE)
        #     append!(vec_RRMSE, RRMSE)
        #     append!(vec_RMSLE, RMSLE)
        #     append!(vec_p, p)
        # end
        ### Parallel execution
        @time vec_vec_metrics = @sync @showprogress @distributed (push!) for j in 1:nfold
            # j= 1
            idx_training = vec_fld .!= j
            idx_validate = vec_fld .== j

            syncx_training = string(syncx, "-CV-rep_", i, "-fold_", j, "-TRAINING.syncx")
            pheno_training = string(syncx, "-CV-rep_", i, "-fold_", j, "-TRAINING.csv")
            syncx_validate = string(syncx, "-CV-rep_", i, "-fold_", j, "-VALIDATE.syncx")
            pheno_validate = string(syncx, "-CV-rep_", i, "-fold_", j, "-VALIDATE.csv")

            SAVE(Window(χ.chr, χ.pos, χ.ref, χ.cou[:, idx_training], zeros(1,1)),  syncx_training)
            SAVE(Phenotype(ϕ.iid[idx_training], [ϕ.tid[1]], ϕ.phe[idx_training, 1:1]), pheno_training, delimiter, ["id", ϕ.tid[1]])
            SAVE(Window(χ.chr, χ.pos, χ.ref, χ.cou[:, idx_validate], zeros(1,1)),  syncx_validate)
            SAVE(Phenotype(ϕ.iid[idx_validate], [ϕ.tid[1]], ϕ.phe[idx_validate, 1:1]), pheno_validate, delimiter, ["id", ϕ.tid[1]])

            tsv = OLS_MULTIVAR(syncx_training, maf, pheno_training, delimiter, header, id_col, phenotype_col, missing_strings, FE_method, out)
            ŷ = PREDICT(tsv, syncx_validate)
            correlation_pearson, correlation_spearman, correlation_kendall, R2, R2_adj, MAE, MBE, RAE, MSE, RMSE, RRMSE, RMSLE, p = CV_METRICS(y[idx_validate], ŷ, y[idx_training])
            [correlation_pearson, correlation_spearman, correlation_kendall, R2, R2_adj, MAE, MBE, RAE, MSE, RMSE, RRMSE, RMSLE, p]
        end
        ### Consolidate metrics
        for x in eachindex(vec_vec_metrics)
            correlation_pearson, correlation_spearman, correlation_kendall, R2, R2_adj, MAE, MBE, RAE, MSE, RMSE, RRMSE, RMSLE, p = vec_vec_metrics[x]
            append!(vec_rep, rep)
            append!(vec_fold, fold)
            append!(vec_correlation_pearson, correlation_pearson)
            append!(vec_correlation_spearman, correlation_spearman)
            append!(vec_correlation_kendall, correlation_kendall)
            append!(vec_R2, R2)
            append!(vec_R2_adj, R2_adj)
            append!(vec_MAE, MAE)
            append!(vec_MBE, MBE)
            append!(vec_RAE, RAE)
            append!(vec_MSE, MSE)
            append!(vec_RMSE, RMSE)
            append!(vec_RRMSE, RRMSE)
            append!(vec_RMSLE, RMSLE)
            append!(vec_p, p)
        end
    end
    ### Output
    file_out = open(out, "a")
    line = string(join(["rep", "fold", "correlation_pearson", "correlation_spearman", "correlation_kendall", "R2", "R2_adj", "MAE", "MBE", "RAE", "MSE", "RMSE", "RRMSE", "RMSLE"], "\t"), "\n")
    write(file_out, line)
    for i in eachindex(vec_rep)
        line = string(join([vec_rep[i], vec_fold[i], vec_correlation_pearson[i], vec_correlation_spearman[i], vec_correlation_kendall[i], vec_R2[i], vec_R2_adj[i], vec_MAE[i], vec_MBE[i], vec_RAE[i], vec_MSE[i], vec_RMSE[i], vec_RRMSE[i], vec_RMSLE[i]], "\t"), "\n")
        write(file_out, line)
    end
    close(file_out)
    return(out)
end






### HYPOTHESIS TESTING
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

function GWAS_PLOT_MANHATTAN(vec_chr::Vector{String}, vec_pos::Vector{Int64}, vec_lod::Vector{Float64}, title::String="")::Tuple{Plots.Plot{Plots.GRBackend}, Vector{String}, Vector{Int64}, Vector{Int64}}
    # n=5; m=100_000; l=135_000_000; k=5; ϵ=Int(1e+15); a=2; vec_chr_lengths=[]; vec_chr_names=[]; dist_noLD=500_000; o=100; t=10; nQTL=5; heritability=0.9
    # @time vec_chr, vec_pos, X, y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability)
    # npools = 5
    # G, p = POOL(X, y, npools)
    # b_hat = G \ p
    # vec_lod = ESTIMATE_LOD(b_hat)
    # title = ""
    ### Find the cummulative locus coordinates across chromosomes for the whole genome and the total number of loci
    chr_names = sort(unique(vec_chr))
    pos_init = [1]
    pos_term = []
    for chr in chr_names
        # chr = sort(unique(vec_chr))[1]
        idx = vec_chr .== chr
        pos = vec_pos[idx]
        if length(pos_term) > 0
            append!(pos_init, pos_term[end] + 1)
        end
        append!(pos_term, (pos_init[end]-1) + maximum(pos))
    end
    l = pos_term[end]
    ### Plot
    p1 = Plots.scatter([0], [0], xlims=[0, l], ylims=[minimum(vec_lod), maximum(vec_lod)], legend=false, markersize=0, markerstrokewidth=0, title=title)
    for chr in chr_names
        # chr = chr_names[1]
        idx = vec_chr .== chr
        pos = vec_pos[idx]
        x = (pos_init[chr_names .== chr][end] - 1) .+ pos
        y = vec_lod[idx]
        Plots.scatter!(x, y, legend=false,
                    markerstrokewidth=0.001,
                    markeralpha=0.4)
    end
    m = length(vec_lod)
    LOD_threshold = -log10(0.05/m)
    Plots.plot!([0,1], [LOD_threshold,LOD_threshold], seriestype=:straightline, legend=false)
    return(p1, chr_names, pos_init, pos_term)
end

function GWAS_PLOT_SIMULATED_QTL!(p::Plots.Plot{Plots.GRBackend}, vec_chr_QTL::Vector{String}, vec_pos_QTL::Vector{Int64}, b::Vector{Float64}, chr_names::Vector{String}, pos_init::Vector{Int64}, pos_term::Vector{Int64})::Plots.Plot{Plots.GRBackend}
    # n=5; m=100_000; l=135_000_000; k=5; ϵ=Int(1e+15); a=2; vec_chr_lengths=[]; vec_chr_names=[]; dist_noLD=500_000; o=100; t=10; nQTL=5; heritability=0.9
    # @time vec_chr, vec_pos, X, y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability)
    # npools = 5
    # G, p = POOL(X, y, npools)
    # b_hat = G \ p
    # vec_lod = ESTIMATE_LOD(b_hat)
    # p, chr_names, pos_init, pos_term = GWAS_PLOT_MANHATTAN(vec_chr, vec_pos, vec_lod)
    # idx = b .!= 0.0
    # vec_chr_QTL = vec_chr[idx]
    # vec_pos_QTL = vec_pos[idx]
    # b = b[idx]
    ### Append positions of the QTL into the manhattan plot
    q = length(vec_chr_QTL)
    for i in 1:q
        # i = 1
        chr = vec_chr_QTL[i]
        pos = vec_pos_QTL[i]
        idx = chr .== chr_names
        x = minimum([(pos_init[idx][end] - 1) .+ pos, pos_term[idx][end]])
        Plots.plot!(p, [x, x], [0, 1], seriestype=:straightline, legend=false)
        Plots.annotate!(p, x, 0, Plots.text(string("β", i, "\n", round(b[i])), 0, 7, :bottom))
    end
    return(p)
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
    p1, chr_names, pos_init, pos_term = GWAS_PLOT_MANHATTAN(vec_chr, vec_pos, vec_lod)
    p2 = GWAS_PLOT_QQ(vec_lod)
    p3 = Plots.plot(p1, p2, layout=(2,1))
    return(p3)
end

### IMPUTE
function IMPUTE!(window::Window, model::String=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance::Bool=true)::Window
    # syncx = "/home/jeffersonfparil/Documents/poolgen/test/test_4.syncx"
    # file = open(syncx, "r")
    # model = "OLS"
    # distance = true
    # window = []
    # for i in 1:100
    #     locus = PARSE(SyncxLine(i, readline(file)))
    #     i = position(file)
    #     push!(window, locus)
    # end
    # window = PARSE(Array{LocusAlleleCounts}(window))
    # model = ["Mean", "OLS", "RR", "LASSO", "GLMNET"][2]
    # distance = true
    # close(file)
    _, p = size(window.cou)
    ### Find the indices of pools with missing data.
    ### These will be used independently and iteratively as our response variables
    idx_pools = sum(ismissing.(window.cou), dims=1)[1,:] .> 0
    ### If we have at least one pool with no missing data, then we proceed with imputation
    if (sum(.!idx_pools) >= 1) & (sum(idx_pools) >= 1)
        ### Explanatory variables
        X = Int.(window.cou[:, .!idx_pools])
         ### Distance covariate (only add if the window is within a single chromosome)
         if (distance) & (length(unique(window.chr))==1) & (model != "Mean")
            m = length(window.pos)
            D = zeros(Int, m, m)
            for i in 1:m
                for j in 1:m
                    D[i,j] = abs(window.pos[i] - window.pos[j])
                end
            end
            Z = MultivariateStats.projection(MultivariateStats.fit(PCA, repeat(D, inner=(7,1)); maxoutdim=1))
            X = hcat(X, (X .!= 0) .* Z) ### multiply Z by (X .!= 0) so we get rid of covariate effects when X_ij is zero
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
                    OLS(hcat(ones(nf), X_train), y_train, "MULTIALGORITHMIC")
                catch
                    try
                        OLS(hcat(ones(nf), X_train), y_train, "CANONICAL")
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
                y_imputed = ceil.(hcat(ones(sum(idx_loci)), X_valid) * β)
                ni = length(y_imputed)
                for i in 1:ni
                    ### For when the predicted counts are too high for Int64 to handle
                    y_imputed[i] =  try
                                        Int64(y_imputed[i])
                                    catch
                                        Int64(maximum(y_train))
                                    end
                end
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
    # syncx = "/home/jeffersonfparil/Documents/poolgen/test/test_3.syncx"
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

end
