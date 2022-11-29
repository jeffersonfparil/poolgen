### GENOME-WIDE ASSOCIATION

####### TEST ########
# include("structs.jl")
# using .structs: PileupLine, SyncxLine, LocusAlleleCounts, Window, PhenotypeLine, Phenotype, MinimisationError
# include("functions_io.jl")
# using ProgressMeter, Distributions
# include("functions_filterTransform.jl")
# using LinearAlgebra, MultivariateStats, GLMNet, Optim
# include("functions_linearModel.jl")
# using Plots
# include("functions_simulate.jl")
#####################

### PARALLELISABLE WRAPPER FUNCTIONS
function OLS_ITERATIVE(syncx::String, init::Int64, term::Int64, maf::Float64, phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_col::Int=2, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], covariate::String=["", "XTX", "COR"][2], covariate_df::Int64=1, out::String="")::String
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
    # maf = 0.01
    # phenotype = "../test/test.csv"
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # covariate::String=["", "XTX", "COR"][2]
    # covariate_df::Int64=1
    # out = ""
    #####################
    ### Output tab-delimeted file including the (1) chromosome name, (2) position, (3) allele, (4) allele frequency, (5) allele effect, an (6) p-value
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-OLS_ITERATIVE.tsv")
    end
    file_out = open(out, "a")
    ### Load phenotype data and standard normalise so we don't need to fit an intercept for simplicity
    ϕ = LOAD(phenotype, delimiter, header, id_col, [phenotype_col], missing_strings)
    y = Float64.(ϕ.phe[:,1])
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
    # init = 114646 ### Passes MAF filter
    seek(file, init)
    ### Prepare covariate
    if covariate != ""
        nloci = 1_000
        maf = 0.001
        δ = 1e-10
        remove_insertions = true
        remove_minor_alleles = false
        remove_correlated_alleles = false
        θ = 0.95
        centre = true
        ρ = try
                GENERATE_COVARIATE(syncx, nloci, covariate_df, covariate, maf, δ, remove_insertions, remove_minor_alleles, remove_correlated_alleles, θ, centre)
            catch
                GENERATE_COVARIATE(syncx, 100, covariate_df, covariate, maf, δ, remove_insertions, remove_minor_alleles, remove_correlated_alleles, θ, centre)
            end
    end
    ### Regress
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
                ### Prepend the covariate
                if covariate != ""
                    x = hcat(ρ, x) ### Keep the allele in the last column
                end
                ### OLS using the canonical method and outputting the estimates of the effects and variances
                x = hcat(ones(n), x) ### Keep the allele in the last column
                n, p = size(x)
                β̂, Vβ̂, Vϵ̂ = OLS(x, y) ### In β̂ the last element refers to the estimated allele effect
                ### Test if we have reasonable residual variance estimate then proceed to output
                if ((Vϵ̂ < 0.0) | (Vϵ̂ == Inf)) == false
                    σβ̂ = []
                    for j in 1:p
                        # j = 1
                        append!(σβ̂, sqrt(Vβ̂[j,j]))
                    end
                    t = β̂ ./ σβ̂ ### Allele t-value is the last element
                    pval = Distributions.ccdf(Distributions.Chisq(p-1), t.^2)
                    line = join([locus.chr[1],
                                 locus.pos[1],
                                 a[i],
                                 mean(x[:,end]),
                                 β̂[end],
                                 pval[end]], "\t")
                    write(file_out, string(line, "\n"))
                end
            end
        end
    end
    close(file)
    close(file_out)
    return(out)
end

function LMM_ITERATIVE(syncx::String, init::Int64, term::Int64, maf::Float64, phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_col::Int=2, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], covariate::String=["", "XTX", "COR"][2], model::String=["GBLUP", "RRBLUP"][1], method::String=["ML", "REML"][1], FE_method::String=["CANONICAL", "N<<P"][2], inner_optimizer=[GradientDescent(), LBFGS(), BFGS(), SimulatedAnnealing(), NelderMead()][1], optim_trace::Bool=false, out::String="")::String
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
    # maf = 0.01
    # phenotype = "../test/test.csv"
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # covariate::String=["", "XTX", "COR"][2]
    # model = ["GBLUP", "RRBLUP"][1]
    # method = ["ML", "REML"][1]
    # FE_method = ["CANONICAL", "N<<P"][2]
    # inner_optimizer = [GradientDescent(), LBFGS(), BFGS(), SimulatedAnnealing(), NelderMead()][1]
    # optim_trace = true
    # out = ""
    #####################
    ### Output tab-delimeted file including the (1) chromosome name, (2) position, (3) allele, (4) allele frequency, (5) allele effect, an (6) p-value
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-LMM_ITERATIVE.tsv")
    end
    file_out = open(out, "a")
    ### Load phenotype data and standard normalise so we don't need to fit an intercept for simplicity
    ϕ = LOAD(phenotype, delimiter, header, id_col, [phenotype_col], missing_strings)
    y = Float64.(ϕ.phe[:,1])
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
    ### Prepare covariate or variance-covariance matrix
    nloci = 1_000
    maf = 0.001
    δ = 1e-10
    remove_insertions = true
    remove_minor_alleles = false
    remove_correlated_alleles = false
    θ = 0.95
    centre = true
    C = try
            GENERATE_COVARIATE(syncx, nloci, n, covariate, maf, δ, remove_insertions, remove_minor_alleles, remove_correlated_alleles, θ, centre)
        catch
            GENERATE_COVARIATE(syncx, 100, n, covariate, maf, δ, remove_insertions, remove_minor_alleles, remove_correlated_alleles, θ, centre)
        end
    ### Assign fixed or random or variance-covariance matrices
    if model == "GBLUP"
        # Assign within the loop: _X_ = hcat(ones(n), x)  ### SNPs have fixed effects
        Z = 1.0*I ### Genotypes have random effects...
        K = C ### ... and are distributed normally μ=0, and Σ=σ2g*K
    elseif model == "RRBLUP"
        _X_ = hcat(ones(n), C)  ### Kinships have effects
        # Assign within the loop: Z = x                   ### SNPs have random effects
        K = 1.0*I
    end
    ### Regress
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
                # p = 1 ### the number of parameters to estimate which is equal to 1 since we're not including an intercept as we centered y
                ### Add the intercept and covariate
                if model == "GBLUP"
                    _X_ = hcat(ones(n), x)  ### SNPs have fixed effects
                    # Assigned outside the loop: Z = 1.0*I ### Genotypes have random effects...
                    # Assigned outside the loop: K = C ### ... and are distributed normally μ=0, and Σ=σ2g*K
                elseif model == "RRBLUP"
                    # Assigned outside the loop: _X_ = hcat(ones(n), C)  ### Kinships have effects
                    Z = x                   ### SNPs have random effects
                    # Assigned outside the loop: K = 1.0*I
                end
                n, p = size(_X_)
                ### Linear mixed model fitting using the canonical method and outputting the estimates of the effects and variances
                σ2u, σ2e = OPTIM_MM(_X_, y, Z, K, FE_method, method, inner_optimizer, optim_trace)
                ### Random effects variance-covariance matrix
                D = σ2u .* K
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
                pval = Distributions.ccdf(Distributions.Chisq(p-1), W)
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

function GWALPHA(syncx::String, py_phenotype::String, init::Int64, term::Int64, maf::Float64, penalty::Bool=true, out::String="")::String
    ####### TEST ########
    # syncx = "../test/test.syncx"
    # py_phenotype = "../test/test.py"
    # file = open(syncx, "r")
    # seekend(file)
    # final_position = position(file)
    # close(file)
    # threads = 4
    # vec_positions = Int.(round.(collect(0:(final_position/threads):final_position)))
    # init = vec_positions[2]
    # term = vec_positions[3]
    # maf = 0.01
    # penalty = true
    # out = ""
    #####################
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
                        ### Optimize (minimize) -log-likelihood of these major allele frequencies modelled as a beta distribution
                        ### using Nelder-Mead optimization or Box minimisation (try-catch if one or the other fails with preference to Nelder-Mead)
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

### HYPOTHESIS TESTING
function BEST_FITTING_DISTRIBUTION(vec_b::Vector{Float64})
    ####### TEST ########
    # vec_b = rand(Distributions.Beta(1.0, 1.5), 100)
    # vec_b = rand(Distributions.Exponential(3.3), 100)
    # vec_b = rand(Distributions.Gamma(1.0, 0.5), 100)
    # vec_b = rand(Distributions.Poisson(1.0), 100)
    # vec_b = rand(Distributions.Normal(0.0, 4.2), 100)
    # vec_b = rand(Distributions.Laplace(0.0, 1.5), 100)
    # vec_b = rand(Distributions.InverseGaussian(0.5, 1.5), 100)
    # vec_b = rand(Distributions.Pareto(4.2, 6.9), 100)
    # vec_b = rand(Distributions.Uniform(-1.5, +2.5), 100)
    # using UnicodePlots
    # UnicodePlots.histogram(vec_b)
    #####################
    DIST_NAMES =   [Distributions.Beta, Distributions.Exponential, Distributions.Gamma, Distributions.Poisson,
                    Distributions.Normal, Distributions.Laplace, Distributions.InverseGaussian,
                    Distributions.Pareto, Distributions.Uniform]
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
    ####### TEST ########
    # syncx = "../test/test.syncx"
    # py_phenotype = "../test/test.py"
    # file = open(syncx, "r")
    # seekend(file)
    # final_position = position(file)
    # close(file)
    # tsv = GWALPHA(syncx, py_phenotype, 0, final_position, 0.001)
    # vec_chr, vec_pos, vec_allele, vec_freq, vec_alpha = LOAD_OUT(tsv, false)
    # vec_b = vec_alpha
    # using UnicodePlots
    # UnicodePlots.histogram(vec_b)
    # UnicodePlots.histogram(pval)
    # UnicodePlots.scatterplot(lod)
    #####################
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
    ####### TEST ########
    # syncx = "../test/test.syncx"
    # py_phenotype = "../test/test.py"
    # file = open(syncx, "r")
    # seekend(file)
    # final_position = position(file)
    # close(file)
    # tsv = GWALPHA(syncx, py_phenotype, 0, final_position, 0.001)
    # vec_chr, vec_pos, vec_allele, vec_freq, vec_alpha = LOAD_OUT(tsv, false)
    # vec_b = vec_alpha
    # D = [Distributions.Beta, Distributions.Exponential, Distributions.Gamma, Distributions.Poisson,
    #      Distributions.Normal, Distributions.Laplace, Distributions.InverseGaussian,
    #      Distributions.Pareto, Distributions.Uniform][5]
    # using UnicodePlots
    # UnicodePlots.histogram(vec_b)
    # UnicodePlots.histogram(pval)
    # UnicodePlots.scatterplot(lod)
    #####################
    distribution = Distributions.fit_mle(D, vec_b)
    pval = [Distributions.cdf(distribution, x) <= Distributions.ccdf(distribution, x) ? 2*Distributions.cdf(distribution, x) : 2*Distributions.ccdf(distribution, x) for x in vec_b]
    lod = -log.(10, pval)
    return(lod)
end

### PLOTTING FUNCTIONS
function GWAS_PLOT_MANHATTAN(vec_chr::Vector{String}, vec_pos::Vector{Int64}, vec_lod::Vector{Float64}, title::String="")::Tuple{Plots.Plot{Plots.GRBackend}, Vector{String}, Vector{Int64}, Vector{Int64}}
    ####### TEST ########
    # syncx = "../test/test.syncx"
    # py_phenotype = "../test/test.py"
    # file = open(syncx, "r")
    # seekend(file)
    # final_position = position(file)
    # close(file)
    # tsv = GWALPHA(syncx, py_phenotype, 0, final_position, 0.001)
    # vec_chr, vec_pos, vec_allele, vec_freq, vec_alpha = LOAD_OUT(tsv, false)
    # vec_b = vec_alpha
    # vec_lod = ESTIMATE_LOD(vec_b)
    # title = "test"
    #####################
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
    # plot_LD = false       ### simulated LD decay
    # vec_chr, vec_pos, X, y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 5
    # G, p = POOL(X, y, npools)
    # θ̂ = GLMNET(G, p, 1.00)
    # vec_b = θ̂[2:end]
    # vec_lod = ESTIMATE_LOD(vec_b)
    # p, chr_names, pos_init, pos_term = GWAS_PLOT_MANHATTAN(vec_chr, vec_pos, vec_lod)
    # idx = b .!= 0.0
    # vec_chr_QTL = vec_chr[idx]
    # vec_pos_QTL = vec_pos[idx]
    # b = b[idx]
    #####################
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
    # plot_LD = false       ### simulated LD decay
    # vec_chr, vec_pos, X, y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # npools = 5
    # G, p = POOL(X, y, npools)
    # θ̂ = GLMNET(G, p, 1.00)
    # vec_b = θ̂[2:end]
    # vec_lod = ESTIMATE_LOD(vec_b)
    # title = "test"
    #####################
    m = length(vec_lod)
    LOD_expected  = -log10.(cdf(Distributions.Uniform(0, 1), (collect(1:m) .- 0.5) ./ m))
    p2 = Plots.scatter(sort(LOD_expected), sort(vec_lod), markerstrokewidth=0.001, markeralpha=0.4, legend=false)
    Plots.plot!([0,1], [0,1], seriestype=:straightline, legend=false, title=title)
    return(p2)
end

function GWAS_PLOT(tsv_gwalpha::String, estimate_empirical_lod::Bool)::Plots.Plot{Plots.GRBackend}
    ####### TEST ########
    # syncx = "../test/test.syncx"
    # py_phenotype = "../test/test.py"
    # file = open(syncx, "r")
    # seekend(file)
    # final_position = position(file)
    # close(file)
    # tsv_gwalpha = GWALPHA(syncx, py_phenotype, 0, final_position, 0.001)
    # estimate_empirical_lod = true
    #####################
    vec_chr, vec_pos, vec_allele, vec_freq, vec_alpha = LOAD_OUT(tsv_gwalpha, false)
    if estimate_empirical_lod
        vec_lod = ESTIMATE_LOD(vec_alpha)
    end
    p1, chr_names, pos_init, pos_term = GWAS_PLOT_MANHATTAN(vec_chr, vec_pos, vec_lod)
    p2 = GWAS_PLOT_QQ(vec_lod)
    p3 = Plots.plot(p1, p2, layout=(2,1))
    return(p3)
end

### BOOTSTRAPPING AND BAYESIAN FUNCTIONS (NOTE: Experimental)
function BOO_ITERATIVE(syncx::String, init::Int64, term::Int64, maf::Float64, phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_col::Int=2, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], F::Function=OLS, F_params="", nburnin::Int64=1_000, δ::Float64=1e-10, maxiter::Int64=1_000, out::String="")::String
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
    # maf = 0.01
    # phenotype = "../test/test.csv"
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # F = [OLS, MM, BOOTSTRAP][1]
    # F_params="N<<P"
    # nburnin = 1_000
    # δ = 1e-10
    # maxiter = 1_000
    # out = ""
    #####################
    ### Output tab-delimeted file including the (1) chromosome name, (2) position, (3) allele, (4) allele frequency, (5) allele effect, an (6) p-value
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-BOO_ITERATIVE.tsv")
    end
    file_out = open(out, "a")
    ### Load phenotype data and standard normalise so we don't need to fit an intercept for simplicity
    ϕ = LOAD(phenotype, delimiter, header, id_col, [phenotype_col], missing_strings)
    y = Float64.(ϕ.phe[:,1])
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
                        ρ = F(x, y, F_params...)[end]
                    else
                        ρ = F(x, y, F_params)[end]
                    end
                else
                    ρ = F(x, y)[end]
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
