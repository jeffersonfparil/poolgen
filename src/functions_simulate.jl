### SIMULATION FUNCTIONS
### Note: using Int64 for the genotype encoding so we can include multi-allelic loci in the future

####### TEST ########
# using Distributions, ProgressMeter, Plots
#####################

function BUILD_FOUNDING_HETEROZYGOUS_GENOMES(n::Int64, m::Int64, l::Int64, k::Int64, ϵ::Int64=Int(1e+15), a::Int64=2, vec_chr_lengths::Vector{Int64}=([0]), vec_chr_names::Vector{String}=[""])::Tuple{Vector{String}, Vector{Int64}, Vector{Int64}, Array{Int64, 3}}
    ####### TEST ########
    # n = 2                 ### number of founders
    # m = 10_000            ### number of loci
    # l = 2_300_000         ### total genome size
    # k = 7                 ### number of chromosomes
    # ϵ = Int(1e+15)        ### an arbitrarily large Int64 number to indicate no LD, i.e. the distance between the termini of 2 adjacent chromosomes is infinitely large LD-wise but we don;t want to use Inf as it is not Int64
    # a = 2                 ### number of alleles per locus (biallelic or a=2 by default)
    # vec_chr_lengths = [0] ### optional vector of chromosome lengths
    # vec_chr_names = [""]  ### optional vector of chromosome names
    #####################
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
    ####### TEST ########
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
    #####################
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
    # Plots.scatter(vec_idx) ### should show blocks of 0's and 1's
    ### Define the resulting allele
    vec_gamete = zeros(Int64, m)
    vec_gamete[  vec_idx] = X[1,   vec_idx] ### homologous chromosome 1
    vec_gamete[.!vec_idx] = X[2, .!vec_idx] ### homologous chromosome 2
    return(vec_gamete)
end

function SIMULATE_GENOMES(G::Array{Int64, 3}, vec_dist::Vector{Int64}, dist_noLD::Int64=10_000, o::Int64=10_000, t::Int64=1)::Array{Int64, 3}
    ####### TEST ########
    # n = 5                   ### number of founders
    # m = 10_000              ### number of loci
    # l = 20_000              ### total genome length
    # k = 1                   ### number of chromosomes
    # ϵ = Int(1e+15)          ### some arbitrarily large number to signify the distance at which LD is nil
    # a = 2                   ### number of alleles per locus
    # vec_chr_lengths = [0]   ### chromosome lengths
    # vec_chr_names = [""]    ### chromosome names 
    # vec_chr, vec_pos, vec_dist, G = BUILD_FOUNDING_HETEROZYGOUS_GENOMES(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names)
    # dist_noLD = 10_000     ### distance at which LD is nil (related to ϵ)
    # o = 1_000               ### total number of simulated individuals
    # t = 10
    #####################
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
    ####### TEST ########
    # n = 10                  ### number of biallelic founders
    # m = 10_000              ### number of loci
    # l = 135_000_000         ### total genome length
    # k = 5                   ### number of chromosomes
    # ϵ = Int(1e+15)          ### some arbitrarily large number to signify the distance at which LD is nil
    # a = 4                   ### number of alleles per locus
    # vec_chr_length = [0]    ### chromosome lengths
    # vec_chr_names = [""]    ### chromosome names
    # dist_noLD = 250_000     ### distance at which LD is nil (related to ϵ)
    # o = 1_000               ### total number of simulated individuals
    # t = 10                  ### number of o-sized random mating generations to simulate
    # vec_chr, vec_pos, vec_dist, F = BUILD_FOUNDING_HETEROZYGOUS_GENOMES(n, m, l, k, ϵ, a, vec_chr_length, vec_chr_names)
    # P = SIMULATE_GENOMES(F, vec_dist, dist_noLD, o, t)
    # chr = ""
    # window_size = 2*dist_noLD
    # n_pairs = 2_000
    # using UnicodePlots
    #####################
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
    # UnicodePlots.scatterplot(Float64.(vec_dist), Float64.(vec_r2))
    return(vec_r2, vec_dist)
end

function SIMULATE(n::Int64, m::Int64, l::Int64, k::Int64, ϵ::Int64=Int(1e+15), a::Int64=2, vec_chr_lengths::Vector{Int64}=Int64.([0]), vec_chr_names::Vector{String}=[""], dist_noLD::Int64=10_000, o::Int64=1_000, t::Int64=10, nQTL::Int64=10, heritability::Float64=0.5, LD_chr::String="", LD_n_pairs::Int64=10_000, plot_LD::Bool=true)::Tuple{Vector{String}, Vector{Int64}, Matrix{Int64}, Vector{Float64}, Vector{Float64}}
    ####### TEST ########
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
    #####################
    ### Instatiante founder genome/s
    vec_chr, vec_pos, vec_dist, G = BUILD_FOUNDING_HETEROZYGOUS_GENOMES(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names)
    ### Simulate random mating with constatnt population sizes for t genrations
    P = SIMULATE_GENOMES(G, vec_dist, dist_noLD, o, t) ### A 3-dimensional array of Int64 where the largest number, e.g.x, means that we have a macimum of x+1 alleles per locus
    _, n, m  = size(P)
    ### Assess LD
    if plot_LD
        time_id = replace(join(collect(string(time()))[1:12], ""), "."=> "")
        LD_window_size = 2*dist_noLD
        vec_r2, vec_dist = LD(P, vec_chr, vec_pos, LD_chr, LD_window_size, LD_n_pairs)
        p = Plots.scatter(vec_dist, vec_r2, legend=false, xlab="Distance (bp)", ylab="r²")
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
    ####### TEST ########
    # n=5; m=100_000; l=135_000_000; k=5; ϵ=Int(1e+15); a=2; vec_chr_lengths=[0]; vec_chr_names=[""]; dist_noLD=500_000; o=100; t=10; nQTL=5; heritability=0.9
    # @time vec_chr, vec_pos, X, y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability)
    # npools = 5
    ####################
    ### Sort the individuals into equally-sized npools pools
    idx = sortperm(y)
    y = y[idx]
    X = X[idx, :]
    n, m = size(X)
    vec_idx_pools = repeat([Int(floor(n/npools))], npools)
    add_to_last_pool = n - sum(vec_idx_pools)
    vec_idx_pools[end] = vec_idx_pools[end] + add_to_last_pool
    vec_idx_pools = append!([0], cumsum(vec_idx_pools))
    ### Setup the exponentially distributed allele counting errors during pooling
    λ = log(n/npools)
    D = Distributions.Exponential(λ)
    ### Generate the matrix of genotype frequencies per pool across m loci and vector of mean phenotypes per pool
    G = zeros(Float64, npools, m)
    p = zeros(Float64, npools)
    for i in 1:npools
        init = vec_idx_pools[i] + 1
        term = vec_idx_pools[i+1]
        ϵ = rand(D)
        G[i, :] = mean(X[init:term, :], dims=1) ./ 2
        p[i] = mean(y[init:term])
    end
    return(G, p)
end

function EXPORT_SIMULATED_DATA(vec_chr::Vector{String}, vec_pos::Vector{Int64}, X::Matrix{Int64}, y::Vector{Float64}, out_geno::String="", out_pheno::String="")::Tuple{String, String, String, String}
    ####### TEST ########
    # n=5; m=100_000; l=135_000_000; k=5; ϵ=Int(1e+15); a=2; vec_chr_lengths=[0]; vec_chr_names=[""]; dist_noLD=500_000; o=100; t=10; nQTL=5
    # @time vec_chr, vec_pos, X, y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL)
    # out_geno = ""
    # out_pheno = ""
    #####################
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
    ####### TEST ########
    # n=5; m=100_000; l=135_000_000; k=5; ϵ=Int(1e+15); a=2; vec_chr_lengths=[0]; vec_chr_names=[""]; dist_noLD=500_000; o=100; t=10; nQTL=5
    # @time vec_chr, vec_pos, X, y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL)
    # @time X, y = POOL(X, y, npools)
    # out_geno = ""
    # out_pheno = ""
    #####################
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
