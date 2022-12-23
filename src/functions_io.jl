### INPUT AND OUTPUT FUNCTIONS

####### TEST ########
# include("structs.jl")
# using .structs: PileupLine, SyncxLine, LocusAlleleCounts, Window, PhenotypeLine, Phenotype
#####################

### INPUT FUNCTIONS
function PARSE(line::PileupLine, minimum_quality=20)::LocusAlleleCounts
    ####### TEST ########
    # fname = "../test/test.pileup"
    # file = open(fname, "r")
    # line = PileupLine(1, readline(file))
    # close(file)
    # minimum_quality=20
    #####################
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
            if (str_state == "+") || (str_state == "-")
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
            elseif ((str_state==".") || 
                    (str_state==",") || 
                    (str_state ∈ ["A", "T", "C", "G"])
                   ) && (q>=minimum_quality)
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
    ####### TEST ########
    # fname = "../test/test.syncx"
    # file = open(fname, "r")
    # line = SyncxLine(1, readline(file))
    # close(file)
    #####################
    lin = map(x -> x=="0:0:0:0:0:0:0" ? "missing:missing:missing:missing:missing:missing:missing" : string(x), split(line.lin, "\t"))
    chr = lin[1]
    pos = parse(Int, lin[2])
    vco = map(x -> x=="missing" ? missing : parse(Int, x), vcat(split.(lin[3:end], ":")...))
    cou = (Base.convert(Matrix{Any}, reshape(vco, 7, Int(length(vco)/7))))
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
    ####### TEST ########
    # fname = "../test/test.pileup"
    # file = open(fname, "r")
    # window = []
    # while !eof(file)
    # # for i in 1:10
    #     push!(window, PARSE(PileupLine(1, readline(file))))
    # end
    # close(file)
    # window = Base.convert(Vector{LocusAlleleCounts}, window)
    #####################
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
    ####### TEST ########
    # fname = "../test/test.csv"
    # file = open(fname, "r")
    # header = PhenotypeLine(1, readline(file), ",", 1, [2])
    # line = PhenotypeLine(1, readline(file), ",", 1, [2])
    # close(file)
    # trait_names = ["Phenotype"]
    # missing_strings = ["NA", ""]
    #####################
    lin = split(line.lin, line.dlm)
    y = map(x -> sum(x .== missing_strings)>0 ? missing : parse(Float64, x), lin[line.trc])
    ids = [string(lin[line.idc])]
    if trait_names==[""]
        trait_names = repeat([""], length(y))
    end
    Y = Base.convert(Array{Any}, reshape(y, 1, length(y)))
    return(Phenotype(ids, trait_names, Y))
end

function PARSE(lines::Vector{Phenotype}, rename_traits::Vector{String}=[""])::Phenotype
    ####### TEST ########
    # fname = "../test/test.csv"
    # file = open(fname, "r")
    # header = PhenotypeLine(1, readline(file), ",", 1, [2])
    # lines = []
    # while !eof(file)
    #     push!(lines, PARSE(PhenotypeLine(1, readline(file), ",", 1, [2])))
    # end
    # close(file)
    # rename_traits = [""]
    #####################
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
    ####### TEST ########
    # fname = "../test/test.pileup"
    # file = open(fname, "r")
    # window = []
    # for i in 1:10
    #     push!(window, PARSE(PileupLine(1, readline(file))))
    # end
    # close(file)
    # window = PARSE(Base.convert(Vector{LocusAlleleCounts}, window))
    # locus = 1
    #####################
    Window([window.chr[locus]],
           [window.pos[locus]],
           [window.ref[locus]],
           window.cou[(7*(locus-1))+1:(locus*7), :],
           window.imp[(7*(locus-1))+1:(locus*7), :])
end

function EXTRACT(window::Window, loci::UnitRange{Int})::Window
    ####### TEST ########
    # fname = "../test/test.pileup"
    # file = open(fname, "r")
    # window = []
    # for i in 1:10
    #     push!(window, PARSE(PileupLine(1, readline(file))))
    # end
    # close(file)
    # window = PARSE(Base.convert(Vector{LocusAlleleCounts}, window))
    # loci = 1:5
    #####################
    Window(window.chr[loci],
           window.pos[loci],
           window.ref[loci],
           window.cou[(7*(loci.start-1))+1:(loci.stop*7), :],
           window.imp[(7*(loci.start-1))+1:(loci.stop*7), :])
end

function SLIDE!(window::Window, locus::LocusAlleleCounts)::Window
    ####### TEST ########
    # fname = "../test/test.pileup"
    # file = open(fname, "r")
    # window = []
    # for i in 1:10
    #     push!(window, PARSE(PileupLine(1, readline(file))))
    # end
    # window = PARSE(Base.convert(Vector{LocusAlleleCounts}, window))
    # locus = PARSE(PileupLine(1, readline(file)))
    # close(file)
    #####################
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

function LOAD(syncx::String, count::Bool)::Window
    ####### TEST ########
    # syncx = "../test/test.syncx"
    # count = false
    #####################
    file = open(syncx, "r")
    loci = []
    while !eof(file)
        push!(loci, PARSE(SyncxLine(1, readline(file))))
    end
    close(file)
    GENOTYPE = PARSE(Base.convert(Vector{LocusAlleleCounts}, loci))
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

function LOAD(phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_cols::Vector{Int}=[2], missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""])::Phenotype
    ####### TEST ########
    # phenotype = "../test/test.csv"
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_cols = [2]
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    #####################
    file = open(phenotype, "r")
    if header
        header_line = Base.convert(Vector{String}, split(readline(file), delimiter))
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
    return(PARSE(Base.convert(Vector{Phenotype}, lines), header_line[phenotype_cols]))
end

function LOAD_OUT(tsv_gwalpha::String, include_all_sites::Bool)::Tuple{Vector{String}, Vector{Int64}, Vector{String}, Vector{Float64}, Vector{Float64}}
    ####### TEST ########
    # tsv_gwalpha = "../test/test-GWAlpha.tsv"
    # include_all_sites = false
    #####################
    file = open(tsv_gwalpha, "r")
    vec_chr = []
    vec_pos = []
    vec_allele = []
    vec_freq = []
    vec_alpha = []
    while !eof(file)
        line = split(readline(file), "\t")
        if include_all_sites || ((include_all_sites == false) && (line[4] != "0"))
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
    ####### TEST ########
    # tsv = "../test/test.tsv"
    #####################
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
    ####### TEST ########
    # py_phenotype = "../test/test.py"
    #####################
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

function LOAD_GFF(gff::String)::Tuple{Vector{String}, Vector{Int64}, Vector{Int64}, Vector{String}}
    ####### TEST ########
    # gff = "../test/test_Lr.gff"
    #####################
    ### Extract annotations
    file = open(gff, "r")
    vec_gff_chr = []
    vec_gff_ini = []
    vec_gff_fin = []
    vec_gff_ann = []
    while !eof(file)
        line = readline(file)
        if line[1] == '#'
            continue
        end
        line = split(line, "\t")
        chr = line[1]
        pos_ini = parse.(Int64, line[4])
        pos_fin = parse.(Int64, line[5])
        ann = line[end]
        if (match(Regex(";genome=chromosome;"), ann) == nothing) && ( match(Regex(";model_evidence="), ann) != nothing)
            push!(vec_gff_chr, chr)
            push!(vec_gff_ini, pos_ini)
            push!(vec_gff_fin, pos_fin)
            push!(vec_gff_ann, ann)
        end
    end
    close(file)
    # hcat(vec_gff_chr, vec_gff_ini, vec_gff_fin, vec_gff_ann)
    return(vec_gff_chr, vec_gff_ini, vec_gff_fin, vec_gff_ann)
end  

function CONVERT(syncx_or_sync::String, init::Int, term::Int, out::String="")::String
    ####### TEST ########
    # syncx_or_sync = "../test/test.syncx"
    # # syncx_or_sync = "../test/test.sync"
    # file = open(syncx_or_sync, "r")
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
        extension = [split(syncx_or_sync, ".")[end][end] == 'x' ? "-converted.sync" : "-converted.syncx"][1]
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
        idx = [1, 2, 3, 4, 7, 6] ### Remove column 5, i.e. INSERTION COUNT COLUMN: from A:T:C:G:I:DEL:N into A:T:C:G:N:DEL
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
    elseif p == 6
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
    else
        println("Unrecognised file format. Please use sync or syncx.")
    end
    close(file)
    close(file_out)
    return(out)
end

### SPLIT AND MERGE FILES
function SPLIT(threads::Int, pileup::String)::Tuple{Vector{Int64}, Vector{Int64}}
    ###### TEST ########
    # threads = 4
    # pileup = "../test/test.pileup"
    #####################
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
    ###### TEST ########
    # threads = 4
    # pileup = "../test/test.pileup"
    # window_size = 20
    #####################
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
end

function MERGE(filenames_out::Vector{String}, out::String)::String
    ###### TEST ########
    # threads = 4
    # syncx = "../test/test.syncx"
    # positions_init, positions_term = SPLIT(threads, syncx)
    # filenames_out = []
    # for i in 1:length(positions_init)
    #     # i = 1
    #     init = positions_init[i]
    #     term = positions_term[i]
    #     tmp = string(syncx, "-CONVERT2SYNCXOR2SYNC-", i, ".tmp")
    #     push!(filenames_out, CONVERT(syncx, init, term, tmp))
    # end
    # out = "../test/test-MERGED.sync"
    #####################
    ### Sort the output files from parallel processing
    sort!(filenames_out)
    ### Merge
    file_out = open(out, "a")
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
    ###### TEST ########
    # threads = 4
    # window_size = 10
    # syncx = "../test/test.syncx"
    # positions_init, positions_term = SPLIT(threads, syncx, window_size)
    # filenames_out = []
    # for i in 1:length(positions_init)
    #     # i = 1
    #     init = positions_init[i]
    #     term = positions_term[i]
    #     tmp = string(syncx, "-CONVERT2SYNCXOR2SYNC-", i, ".tmp")
    #     push!(filenames_out, CONVERT(syncx, init, term, tmp))
    # end
    # out = "../test/test-MERGED-FOR_OVERLAPPING_WINDOW.sync"
    #####################
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
            while (!eof(file_in)) && (j < max_line)
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

### SAVING FUNCTIONS
function SAVE(line::PileupLine, filename::String)
    ####### TEST ########
    # fname = "../test/test.pileup"
    # file = open(fname, "r")
    # line = PileupLine(1, readline(file))
    # filename = "../test/test-SAVE.pileup"
    #####################
    file = open(filename, "a")
    write(file, string(line.lin, '\n'))
    close(file)
end

function SAVE(window::Window, filename::String)
    ####### TEST ########
    # fname = "../test/test.pileup"
    # file = open(fname, "r")
    # window = []
    # for i in 1:10
    #     push!(window, PARSE(PileupLine(1, readline(file))))
    # end
    # close(file)
    # window = PARSE(Base.convert(Vector{LocusAlleleCounts}, window))
    # filename = "../test/test-SAVE_WINDOW.syncx"
    #####################
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
    ####### TEST ########
    # phenotype = Phenotype(["individual_1", "individual_2", "individual_3", "individual_4", "individual_5"],
    #                       ["Random_phenotype"],
    #                       rand(5, 1))
    # filename = "../test/test-SAVE.csv"
    # delimiter = ","
    # header = ["id", "Random_phenotype"]
    # missing_string = "NA"
    #####################
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

function AVOID_FILE_OVERWRITE(out::String)::String
    ####### TEST ########
    # out = "../test/test.tsv"
    #####################
    if isfile(out)
        out_basename = join(split(out, ".")[1:(end-1)], ".")
        out_extension = split(out, ".")[end]
        out = string(out_basename, "-", Dates.now(Dates.UTC), ".", out_extension)
    end
    return(out)
end