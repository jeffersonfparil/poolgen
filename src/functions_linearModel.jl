### LINEAR MODELS

####### TEST ########
# include("structs.jl")
# using .structs: PileupLine, SyncxLine, LocusAlleleCounts, Window, PhenotypeLine, Phenotype, MinimisationError
# include("functions_io.jl")
# using ProgressMeter, Distributions
# include("functions_filterTransform.jl")
# using LinearAlgebra, MultivariateStats, GLMNet, Optim
#####################

### HELPER FUNCTIONS
function GENERATE_COVARIATE(X::Array{T}, df::Int64=1, covariate::String=["XTX", "COR"][2])::Array{T} where T <: Number
    ####### TEST ########
    # syncx = "../test/test.syncx"
    # χ = LOAD(syncx, false)
    # X = Base.convert(Array{Float64}, χ.cou')
    # df = 1
    # covariate = "XTX"
    #####################
    n, p = size(X)
    ### Calculate covariate
    if covariate == "XTX"
        ### (1) X'*X/n covariate
        C = (X * X') ./ p
    elseif covariate == "COR"
        ### (2) Pearson's product-moment correlation
        C = cor(X')
    end
    ### Use the PCs if the we're asking for columns (i.e. df) less than the number of columns in C
    if (df < n)
        C = MultivariateStats.projection(MultivariateStats.fit(PCA, C; maxoutdim=df))
    end
    return(C)
end

function GENERATE_COVARIATE(syncx::String, nloci::Int64=1_000, df::Int64=1, covariate::String=["XTX", "COR"][2], maf::Float64=0.001, δ::Float64=1e-10, remove_insertions::Bool=true, remove_minor_alleles::Bool=false, remove_correlated_alleles::Bool=false, θ::Float64=0.95, centre::Bool=true)::Array{Float64}
    ####### TEST ########
    # syncx = "../test/test.syncx"
    # nloci = 100
    # df = 1
    # covariate = "XTX"
    # maf = 0.0001
    # δ = 1e-10
    # remove_insertions = true
    # remove_minor_alleles = false
    # remove_correlated_alleles = true
    # θ = 0.99
    # centre = true
    #####################
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
        χ = PARSE(Base.convert(Vector{LocusAlleleCounts}, loci))
    end
    ### Filter
    X, vec_chr, vec_pos, vec_ref = FILTER(χ, maf, δ, remove_insertions, remove_minor_alleles, remove_correlated_alleles, θ, centre)
    ### Generate covariate
    C = GENERATE_COVARIATE(X, df, covariate)
    return(C)
end

function INVERSE(A::Union{Array{T}, UniformScaling{T}, Cholesky{T, Matrix{T}}})::Union{Array{T}, UniformScaling{T}} where T <: Number
    ####### TEST ########
    # A = rand(10, 10)
    #####################
    try
        inv(A)
    catch
        pinv(A)
    end
end

### BASIC LINEAR MODEL FITTING FUNCTIONS
function OLS(X::Array{T}, y::Array{T}, FE_method::String)::Array{Float64} where T <: Number
    ####### TEST ########
    # syncx = "../test/test.syncx"
    # csv = "../test/test.csv"
    # χ = LOAD(syncx, false)
    # X, vec_chr, vec_pos, vec_ref = FILTER(χ)
    # ϕ = LOAD(csv, ",", true, 1, [2])
    # y = ϕ.phe[:,1]
    # FE_method = ["CANONICAL", "MULTIALGORITHMIC", "N<<P"][3]
    #####################
    if isa(X, Vector)
        X = reshape(X, length(X), 1)
    end
    n, p = size(X)
    if X[:,1] != ones(n)
        X = hcat(ones(n), X)
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

function OLS(X::Array{T}, y::Array{T})::Tuple{Array{Float64}, Array{Float64}, Float64} where T <: Number
    ####### TEST ########
    # syncx = "../test/test.syncx"
    # csv = "../test/test.csv"
    # χ = LOAD(syncx, false)
    # X, vec_chr, vec_pos, vec_ref = FILTER(χ)
    # X = X[:, 1:100]
    # ϕ = LOAD(csv, ",", true, 1, [2])
    # y = ϕ.phe[:,1]
    #####################
    if isa(X, Vector)
        X = reshape(X, length(X), 1)
    end
    n, p = size(X)
    if X[:,1] != ones(n)
        X = hcat(ones(n), X)
    end
    V = INVERSE(X' * X)
    β̂ = V * (X' * y)
    ε̂ = y - (X * β̂)
    Vϵ̂ = (ε̂' * ε̂) / (n-p)
    Vβ̂ = Vϵ̂ * V
    return(β̂, Vβ̂, Vϵ̂)
end

function GLMNET(X::Array{T}, y::Array{T}, alpha::Float64=1.0)::Vector{T} where T <: Number
    ####### TEST ########
    # syncx = "../test/test_Lr.syncx"
    # csv = "../test/test_Lr.csv"
    # χ = LOAD(syncx, false)
    # X, vec_chr, vec_pos, vec_ref = FILTER(χ, 0.01, 1e-3, true, false, true)
    # ϕ = LOAD(csv, ",", true, 1, [2])
    # y = Float64.(ϕ.phe[:,1])
    # alpha = 0.50
    #####################
    ### Elastic net regularisation, where ridge if alpha=0.0, and LASSO if alpha=1.0
    if isa(X, Vector)
        X = reshape(X, length(X), 1)
    end
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
    return(β̂)
end

function MM(X::Array{T}, y::Array{T}, Z::Union{Array{T}, UniformScaling{T}}, D::Union{Array{T}, UniformScaling{T}}, R::Union{Array{T}, UniformScaling{T}}, FE_method::String)::Tuple{Vector{Float64}, Vector{Float64}} where T <: Number
    ####### TEST ########
    # syncx = "../test/test_Lr.syncx"
    # csv = "../test/test_Lr.csv"
    # χ = LOAD(syncx, false)
    # G, vec_chr, vec_pos, vec_ref = FILTER(χ, 0.01, 1e-3, true, false, true)
    # ϕ = LOAD(csv, ",", true, 1, [2])
    # y = Float64.(ϕ.phe[:,1])
    # K = GENERATE_COVARIATE(G, size(G,1))
    # ### gBLUP
    # X = G
    # Z = 1.00*I
    # D = 0.10*K
    # R = 0.01*I
    # FE_method = ["CANONICAL", "N<<P"][2]
    # ### rrBLUP
    # X = hcat(ones(size(G,1)))
    # Z = G
    # D = 0.10*I
    # R = 0.01*I
    # FE_method = ["CANONICAL", "N<<P"][2]
    #####################
    if isa(X, Vector)
        X = reshape(X, length(X), 1)
    end
    n, p = size(X)
    if X[:,1] != ones(n)
        X = hcat(ones(n), X)
    end
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

function MM(X::Array{T}, y::Array{T}, Z::Union{Array{T}, UniformScaling{T}}, D::Union{Array{T}, UniformScaling{T}}, R::Union{Array{T}, UniformScaling{T}})::Tuple{Vector{Float64}, Vector{Float64}, Array{T}} where T <: Number
    ####### TEST ########
    # syncx = "../test/test_Lr.syncx"
    # csv = "../test/test_Lr.csv"
    # χ = LOAD(syncx, false)
    # G, vec_chr, vec_pos, vec_ref = FILTER(χ, 0.01, 1e-3, true, false, true)
    # ϕ = LOAD(csv, ",", true, 1, [2])
    # y = Float64.(ϕ.phe[:,1])
    # K = GENERATE_COVARIATE(G, size(G,1))
    # ### gBLUP
    # X = G
    # Z = 1.00*I
    # D = 0.10*K
    # R = 0.01*I
    # ### rrBLUP
    # X = hcat(ones(size(G,1)))
    # Z = G
    # D = 0.10*I
    # R = 0.01*I
    #####################
    ### Linear mixed model fitting which outputs the fixed effects variances for iterative regression
    if isa(X, Vector)
        X = reshape(X, length(X), 1)
    end
    n, p = size(X)
    if X[:,1] != ones(n)
        X = hcat(ones(n), X)
    end
    ### Variance-covariance matrix of y (combined inverse of R and D in Henderson's mixed model equations)
    V = (Z * D * Z') + R
	VI = INVERSE(V)
    Σ̂ = INVERSE(X' * VI * X) ### Using canonical equations to estimate fixed effects so we estimate the fixed effect variances
	### Mixed model equations
    β̂ = Σ̂ * (X' * VI * y)
    μ̂ = (D * Z' * VI) * (y - (X*β̂))
    return(β̂, μ̂, Σ̂)
end

function NLL_MM(θ::Vector{T}, X::Array{T}, y::Array{T}, Z::Union{Array{T}, UniformScaling{T}}, K::Union{Array{T}, UniformScaling{T}}, FE_method::String=["CANONICAL", "N<<P"][2], method::String=["ML", "REML"][1])::Float64 where T <: Number
    ####### TEST ########
    # syncx = "../test/test_Lr.syncx"
    # csv = "../test/test_Lr.csv"
    # θ = [0.1, 0.01]
    # χ = LOAD(syncx, false)
    # G, vec_chr, vec_pos, vec_ref = FILTER(χ, 0.01, 1e-3, true, false, true)
    # ϕ = LOAD(csv, ",", true, 1, [2])
    # y = Float64.(ϕ.phe[:,1])
    # K = GENERATE_COVARIATE(G, size(G,1))
    # FE_method = ["CANONICAL", "N<<P"][2]
    # method = ["ML", "REML"][1]
    # T = Float64
    # ### gBLUP
    # X = G
    # Z = 1.0*I
    # K = 1.0*I
    # ### rrBLUP
    # X = hcat(ones(size(G,1)))
    # Z = G
    # K = K
    #####################
    if isa(X, Vector)
        X = reshape(X, length(X), 1)
    end
    n, p = size(X)
    if X[:,1] != ones(n)
        X = hcat(ones(n), X)
    end
    ### Homoscedastic error and random effect variances
	σ2u = θ[1]          # variance of the other random effects (assuming homoscedasticity)
    σ2e = θ[2]          # variance of the error effects (assuming homoscedasticity)
	n = length(y)		# number of individual samples
	l = size(X, 2)		# number of fixed effects
	nz, pz = try
            size(Z)
        catch
            (n, n)
        end
    ### Random effects variance-covariance matrix
    D = σ2u .* K
    ### Error variance-covariance matrix (homoscedastic)
    # R = diagm(repeat([σ2e], n));
    R = σ2e*I;
    ### Variance-covariance matrix of y (combined inverse of R and D in Henderson's mixed model equations)
    V = (Z * D * Z') + R;
    if typeof(V) == UniformScaling{T}
        V = diagm(repeat([V.λ], n))
    end
    β̂, μ̂ = MM(X, y, Z, D, R, FE_method);
    ### Calculation negative log-likelihoods of variance θ
    if method == "ML"
        ### The negative log-likelihood function y given σ2e and σ2u
        μy = (X * β̂)
        neg_log_lik = 0.5 * ( log(abs(det(V))) + ((y - μy)' * INVERSE(V) * (y - μy)) + (n*log(2*pi)) )
    elseif method == "REML"
        if nz == pz
            ### NOTE: Z MUST BE SQUARE!
            A = σ2u*Z + R
            if typeof(A) == UniformScaling{T}
                A = diagm(repeat([A.λ], n))
            end
            M = try
                    INVERSE(LinearAlgebra.cholesky(A)) ### Cholesky decomposition
                catch
                    INVERSE(LinearAlgebra.lu(A).L) ### LU decomposition
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

function OPTIM_MM(X::Array{T}, y::Array{T}, Z::Union{Array{T}, UniformScaling{T}}, K, FE_method::String=["CANONICAL", "N<<P"][2], method::String=["ML", "REML"][1], inner_optimizer=[GradientDescent(), LBFGS(), BFGS(), SimulatedAnnealing(), NelderMead()][1], optim_trace::Bool=false)::Tuple{Float64, Float64} where T <: Number
    ####### TEST ########
    # syncx = "../test/test_Lr.syncx"
    # csv = "../test/test_Lr.csv"
    # θ = [0.1, 0.01]
    # χ = LOAD(syncx, false)
    # G, vec_chr, vec_pos, vec_ref = FILTER(χ, 0.01, 1e-3, true, false, true)
    # ϕ = LOAD(csv, ",", true, 1, [2])
    # y = Float64.(ϕ.phe[:,1])
    # K = GENERATE_COVARIATE(G, size(G,1))
    # FE_method = ["CANONICAL", "N<<P"][2]
    # method = ["ML", "REML"][1]
    # inner_optimizer = [LBFGS(), BFGS(), SimulatedAnnealing(), GradientDescent(), NelderMead()][1]
    # optim_trace = true
    # T = Float64
    # ### gBLUP
    # X = G
    # Z = 1.0*I
    # K = 1.0*I
    # ### rrBLUP
    # X = hcat(ones(size(G,1)))
    # Z = G
    # K = K
    #####################
    if isa(X, Vector)
        X = reshape(X, length(X), 1)
    end
    n, p = size(X)
    if X[:,1] != ones(n)
        X = hcat(ones(n), X)
    end
    ### Find the optimum β̂ and μ̂
    lower_limits = [1e-5, 1e-5]
    upper_limits = [1e+5, 1e+5]
    initial_values = [1.0, 1.0]
    θ = try
            Optim.optimize(parameters->NLL_MM(parameters, X, y, Z, K, FE_method, method),
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
                                        time_limit=60))
        catch
            MinimisationError(initial_values, false, false)
        end
    σ2u = θ.minimizer[1] # variance of the other random effects (assuming homoscedasticity)
    σ2e = θ.minimizer[2] # variance of the error effects (assuming homoscedasticity)
    # ### Output messages
    if optim_trace
        @show θ
        if (θ.f_converged) | (θ.g_converged)
            println("CONVERGED! 😄")
        else
            println("DID NOT CONVERGE! 😭")
        end
    end
    return(σ2u, σ2e)
end

function MM(X::Array{T}, y::Array{T}, model::String=["GBLUP", "RRBLUP"][1], method::String=["ML", "REML"][1], inner_optimizer=["GradientDescent", "LBFGS", "BFGS", "SimulatedAnnealing", "NelderMead"][1], optim_trace::Bool=false, FE_method::String=["CANONICAL", "N<<P"][2], GBLUP_K::String=["XTX", "COR"][2])::Array{T} where T <: Number
    ####### TEST ########
    # syncx = "../test/test_Lr.syncx"
    # csv = "../test/test_Lr.csv"
    # θ = [0.1, 0.01]
    # χ = LOAD(syncx, false)
    # X, vec_chr, vec_pos, vec_ref = FILTER(χ, 0.01, 1e-3, true, false, true)
    # ϕ = LOAD(csv, ",", true, 1, [2])
    # y = Float64.(ϕ.phe[:,1])
    # model = ["GBLUP", "RRBLUP"][2]
    # method = ["ML", "REML"][1]
    # inner_optimizer = ["GradientDescent", "LBFGS", "BFGS", "SimulatedAnnealing", "NelderMead"][1]
    # optim_trace = false
    # FE_method = ["CANONICAL", "N<<P"][2]
    # GBLUP_K = ["XTX", "COR"][2]
    # T = Float64
    #####################
    if isa(X, Vector)
        X = reshape(X, length(X), 1)
    end
    n, p = size(X)
    if X[:,1] != ones(n)
        X = hcat(ones(n), X)
    end
    if model == "GBLUP"
        Z = 1.0*I                             ### Genotypes have random effects...
        K = GENERATE_COVARIATE(X, n, GBLUP_K)    ### ... and are distributed normally μ=0, and Σ=σ2g*K
    elseif model == "RRBLUP"
        Z = X[:,2:end]                               ### SNPs have random effects...
        # K = diagm(repeat([1.0], size(Z,2))) ### ... and are spherically distributed proportional to σ2g
        K = 1.0*I
        X = X[:,1:1]                      ### Intercept is the only fixed effect
    else
        println(string("Sorry ", model, " is not implemented."))
        println("Please choose from 'GBLUP' and 'RRBLUP'.")
    end
    ### Define the inner optimiser
    if inner_optimizer == "LBFGS"
        inner_optimizer = Optim.LBFGS()
    elseif inner_optimizer == "BFGS"
        inner_optimizer = Optim.BFGS()
    elseif inner_optimizer == "SimulatedAnnealing"
        inner_optimizer = Optim.SimulatedAnnealing()
    elseif inner_optimizer == "GradientDescent"
        inner_optimizer = Optim.GradientDescent()
    elseif inner_optimizer == "NelderMead"
        inner_optimizer = Optim.NelderMead()
    end
    ### Linear mixed model fitting using the canonical method and outputting the estimates of the effects and variances
    σ2u, σ2e = OPTIM_MM(X, y, Z, K, FE_method, method, inner_optimizer, optim_trace)
    ### Random effects variance-covariance matrix
    D = σ2u*K
    ### Error variance-covariance matrix (homoscedastic)
    R = σ2e*I
    ### Solve the mixed model equations
    β̂, μ̂ = MM(X, y, Z, D, R, FE_method)
    if model == "GBLUP"
        θ̂ = β̂
    elseif model == "RRBLUP"
        θ̂ = vcat(β̂, μ̂) # include the fixed intercept effect
    end
    return(θ̂)
end

function NLL_BETA(beta::Array{T}, percA::Array{T}, percB::Array{T}) where T <: Number
    ####### TEST ########
    # beta = [1.0,
    #         1.0,
    #         1.0,
    #         0.5]
    # n = 5
    # freqA = sort(rand(Distributions.Beta(1.0, 0.5), n))
    # bins = repeat([0.2], n)
    # pA = sum(freqA .* bins)
    # pB = 1 - pA
    # binA = (freqA .* bins) ./ pA
    # binB = ( (1 .- freqA) .* bins ) ./ (1-pA)
    # percA = cumsum(binA)
    # percB = cumsum(binB)
    #####################
	-sum(
		 log.(10,
		 	 (Distributions.cdf.(Distributions.Beta(beta[1], beta[2]), percA)) .-
		 	 (Distributions.cdf.(Distributions.Beta(beta[1], beta[2]), append!([0.0], percA[1:(length(percA)-1)])))
		     )
		 )     -sum(
		 log.(10,
		 	 (Distributions.cdf.(Distributions.Beta(beta[3], beta[4]), percB)) .-
		 	 (Distributions.cdf.(Distributions.Beta(beta[3], beta[4]), append!([0.0], percB[1:(length(percB)-1)])))
		     )
		 )
end

### BOOTSTRAPPING AND BAYESIAN FUNCTIONS (NOTE: Experimental)
function BOOTSTRAP(ρ, x::Array{T}, y::Array{T}, F::Function, F_params="")::Tuple{Int64, Int64} where T <: Number
    ####### TEST ########
    # n = 5; p = 1
    # x = rand(n, p)
    # y = rand(n)
    # ### Pearson's product moment correlation
    # F = cor
    # F_params = ""
    # ρ = F(x, y)[end]
    # ### OLS
    # F = OLS
    # F_params = "CANONICAL"
    # ρ = F(x, y, F_params)[end]
    #####################
    n = length(x)
    niter = minimum([100, Int64(n*(n-1)/2)])
    vec_ρ = []
    for i in 1:niter
        x_rand = sample(x, n; replace=false)
        y_rand = sample(y, n; replace=false)
        if F_params != ""
            if isa(F_params, Vector)
                append!(vec_ρ, F(x_rand, y_rand, F_params...)[end])
            else
                append!(vec_ρ, F(x_rand, y_rand, F_params)[end])
            end
        else
            append!(vec_ρ, F(x_rand, y_rand)[end])
        end
    end
    positive_cases = sum(abs.(vec_ρ) .>= abs(ρ))
    return(positive_cases, niter)
end

function BOOTSTRAP_PVAL(ρ, x::Array{T}, y::Array{T}, F::Function, F_params="", nburnin::Int64=10_000, δ::Float64=1e-10, maxiter::Int64=1_000)::Vector{Float64} where T <: Number
    ####### TEST ########
    # n = 5
    # x = rand(n)
    # y = x * 10
    # F = cor
    # ρ = F(x, y)
    # nburnin = 10_000
    # δ = 1e-10
    # maxiter = 1_000
    #####################
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

function ACCEPT_OR_REJECT!(β̂::Array{T}, X::Array{T}, y::Array{T}, θ̂::Float64, j::Int64)::Array{T} where T <: Number
    ####### TEST ########
    # syncx = "../test/test.syncx"
    # csv = "../test/test.csv"
    # χ = LOAD(syncx, false)
    # ϕ = LOAD(csv, ",", true, 1, [2])
    # X, vec_chr, vec_pos, vec_ref = FILTER(χ)
    # y = ϕ.phe[:,1]
    # n, p = size(X)
    # θ̂ = 0.1
    # β̂ = rand(Distributions.Chisq(θ̂), p)
    # β̂[β̂ .< 1e-3] .= 0.00
    # j = 3
    #####################
    ϵ0 = mean((y - X * β̂).^2)
    b_original = β̂[j]
    b_proposed = rand(Distributions.Chisq(θ̂), 1)[1]
    β̂[j] = b_proposed
    ϵ1 = mean((y - (X * β̂)).^2)
    if ϵ0 < ϵ1
        β̂[j] = b_original
    else
        ϵ1 = ϵ0
    end
    return(β̂)
end

# function ABC(X::Array{T}, y::Array{T})::Array{T} where T <: Number
#     # n = 5                 ### number of founders
#     # m = 10_000            ### number of loci
#     # l = 135_000_000       ### total genome length
#     # k = 5                 ### number of chromosomes
#     # ϵ = Int(1e+15)        ### some arbitrarily large number to signify the distance at which LD is nil
#     # a = 2                 ### number of alleles per locus
#     # vec_chr_lengths = [0] ### chromosome lengths
#     # vec_chr_names = [""]  ### chromosome names 
#     # dist_noLD = 500_000   ### distance at which LD is nil (related to ϵ)
#     # o = 1_000               ### total number of simulated individuals
#     # t = 10                ### number of random mating constant population size generation to simulate
#     # nQTL = 10             ### number of QTL to simulate
#     # heritability = 0.5    ### narrow(broad)-sense heritability as only additive effects are simulated
#     # LD_chr = ""           ### chromosome to calculate LD decay from
#     # LD_n_pairs = 10_000   ### number of randomly sampled pairs of loci to calculate LD
#     # plot_LD = false       ### simulate# plot simulated LD decay
#     # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
#     # npools = 100
#     # @time X, y = POOL(_X, _y, npools)
#     # using UnicodePlots
#     # BLAS.set_num_threads(length(Sys.cpu_info())-1)
#     df = 1
#     nburnin = 100
#     niter = 1_000
#     n, p = size(X)

#     ### A single prior
#     θ̂ = 1.0
#     β̂ = rand(Distributions.Exponential(θ̂), p)
#     @showprogress for i in 1:(nburnin+niter)
#         @simd for j in 1:p
#             # j = 1
#             ACCEPT_OR_REJECT!(β̂, X, y, θ̂, j)
#         end
#         θ̂ = Distributions.fit_mle(Distributions.Exponential, β̂).θ
#     end

#     UnicodePlots.histogram(Float64.(β̂))

#     @time β̂_glmnet = GLMNET(X, y); UnicodePlots.scatterplot(Float64.(b), Float64.(β̂_glmnet[2:end])); cor(Float64.(b), Float64.(β̂_glmnet[2:end]))

#     UnicodePlots.scatterplot(Float64.(b))
#     UnicodePlots.scatterplot(Float64.(β̂))
#     UnicodePlots.scatterplot(Float64.(β̂_glmnet))
#     UnicodePlots.scatterplot(Float64.(b), Float64.(β̂))
#     UnicodePlots.scatterplot(Float64.(y), Float64.(X * β̂))
#     UnicodePlots.scatterplot(Float64.(b), Float64.(β̂_glmnet[2:end]))
#     UnicodePlots.scatterplot(Float64.(β̂), Float64.(β̂_glmnet[2:end]))

#     cor(b, β̂)
#     cor(b, β̂_glmnet[2:end])


#     ####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#     ####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#     ####@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
#     ### Multiple priors


#     function BEST_FITTING_DISTRIBUTION(vec_b::Vector{Float64})
#         DIST_NAMES =   [Distributions.Bernoulli, Distributions.Beta, Distributions.Binomial, Distributions.Categorical,
#                         Distributions.DiscreteUniform, Distributions.Exponential, Distributions.Normal, Distributions.Gamma,
#                         Distributions.Geometric, Distributions.Laplace, Distributions.Pareto, Distributions.Poisson,
#                         Distributions.InverseGaussian, Distributions.Uniform]
#         DIST_INSTANCES = [try Distributions.fit_mle(D, vec_b); catch nothing; end for D in DIST_NAMES]
#         NEG_LOGLIK = [try -sum(Distributions.logpdf.(D, vec_b)); catch nothing; end for D in DIST_INSTANCES]
#         DISTRIBUTIONS_DF = hcat((DIST_NAMES[NEG_LOGLIK .!= nothing],
#                                 DIST_INSTANCES[NEG_LOGLIK .!= nothing],
#                                 NEG_LOGLIK[NEG_LOGLIK .!= nothing])...)
#         D = try
#             (DISTRIBUTIONS_DF[argmin(DISTRIBUTIONS_DF[:,3]), 2], DISTRIBUTIONS_DF[argmin(DISTRIBUTIONS_DF[:,3]), 1])
#         catch
#             (nothing, "Failure to fit into any of the distributions tested.")
#         end
#         return(D)
#     end

    
#     θ̂ = 1.0
#     β̂ = rand(Distributions.Exponential(θ̂), p)
#     B = 0
#     @showprogress for i in 1:nburnin
#         β̂ = rand(Distributions.Chisq(df), p)
#         @simd for j in 1:p
#             # j = 1
#             ACCEPT_OR_REJECT!(β̂, X, y, θ̂, j)
#         end
#         θ̂ = Distributions.fit_mle(Distributions.Exponential, β̂).θ
#         if B == 0
#             B = β̂
#         else
#             B = hcat(B, β̂)
#         end
#     end

#     ### Setup prior distributions for each allele
#     vec_mus = std(B, dims=2)[:,1]
#     vec_std = std(B, dims=2)[:,1]
#     vec_dist = []
#     @showprogress for i in 1:p
#         # i = 1
#         if vec_std[i] == 0
#             d = missing
#         else
#             d = BEST_FITTING_DISTRIBUTION(B[i, :])[1]
#         end
#         push!(vec_dist, d)
#     end
#     UnicodePlots.scatterplot(Float64.(b))
#     UnicodePlots.scatterplot(Float64.(β̂))
#     UnicodePlots.histogram(Float64.(β̂))
    
#     UnicodePlots.scatterplot(Float64.(B[:,1]))
#     UnicodePlots.scatterplot(Float64.(B[:,2]))
#     UnicodePlots.scatterplot(Float64.(B[:,3]))
#     UnicodePlots.scatterplot(Float64.(B[:,4]))
#     UnicodePlots.scatterplot(Float64.(B[:,5]))
#     UnicodePlots.scatterplot(Float64.(B[:,6]))
#     UnicodePlots.scatterplot(Float64.(B[:,7]))


#     ### Prepare preliminary β̂
#     β̂ = []
#     @showprogress for j in 1:p
#         if ismissing(vec_dist[j])
#             append!(β̂, vec_mus[j])
#         else
#             append!(β̂, rand(vec_dist[j], 1)[1])
#         end
#     end

#     ### Gibbs sampling
#     vec_idx = collect(1:p)[.!ismissing.(vec_dist)]
#     @showprogress for i in 1:niter
#         # i = 1
#         for j in vec_idx
#             # j = 1
#             ϵ0 = mean((y - X * β̂).^2)
#             b_original = β̂[j]
#             b_proposed = rand(vec_dist[j], 1)[1]
#             β̂[j] = b_proposed
#             ϵ1 = mean((y - (X * β̂)).^2)
#             if ϵ0 < ϵ1
#                 β̂[j] = b_original
#             else
#                 append!(B[j, 2:end], β̂[j])
#                 vec_dist[j] = BEST_FITTING_DISTRIBUTION(B[i, :])[1] 
#             end
#         end
#     end
#     UnicodePlots.scatterplot(Float64.(y), Float64.(X * β̂))
#     UnicodePlots.histogram(Float64.(β̂))
#     UnicodePlots.histogram(Float64.(b))
#     UnicodePlots.scatterplot(Float64.(b), Float64.(β̂))
#     cor(Float64.(b), Float64.(β̂))
#     UnicodePlots.scatterplot(Float64.(b))
#     UnicodePlots.scatterplot(Float64.(β̂))
#     UnicodePlots.scatterplot(Float64.(β̂_glmnet))

    
#     @time β̂_glmnet = GLMNET(X, y); UnicodePlots.scatterplot(Float64.(b), Float64.(β̂_glmnet[2:end])); cor(Float64.(b), Float64.(β̂_glmnet[2:end]))



#     # using DecisionTree
#     # n, m = 10^3, 5
#     # features = randn(n, m)
#     # weights = rand(-2:2, m)
#     # labels = features * weights
#     # # train regression forest, using 2 random features, 10 trees,
#     # # averaging of 5 samples per leaf, and 0.7 portion of samples per tree
#     # n_random_features = 2
#     # n_trees = 10
#     # frac_samples_per_tree = 0.7
#     # n_samples_per_leaf = 5
#     # model = DecisionTree.build_forest(labels, features, n_random_features, n_trees, frac_samples_per_tree, n_samples_per_leaf)
#     # # # apply learned model
#     # # DecisionTree.apply_forest(model, [-0.9,3.0,5.1,1.9,0.0])
#     # # run 3-fold cross validation on regression forest, using 2 random features per split
#     # n_subfeatures=2; n_folds=3
#     # r2 = DecisionTree.nfoldCV_forest(labels, features, n_folds, n_subfeatures)

#     # # set of regression build_forest() parameters and respective default values
#     # # n_subfeatures: number of features to consider at random per split (default: -1, sqrt(# features))
#     # # n_trees: number of trees to train (default: 10)
#     # # partial_sampling: fraction of samples to train each tree on (default: 0.7)
#     # # max_depth: maximum depth of the decision trees (default: no maximum)
#     # # min_samples_leaf: the minimum number of samples each leaf needs to have (default: 5)
#     # # min_samples_split: the minimum number of samples in needed for a split (default: 2)
#     # # min_purity_increase: minimum purity needed for a split (default: 0.0)
#     # # keyword rng: the random number generator or seed to use (default Random.GLOBAL_RNG)
#     # #              multi-threaded forests must be seeded with an `Int`

#     # n = 5                 ### number of founders
#     # m = 10_000            ### number of loci
#     # l = 135_000_000       ### total genome length
#     # k = 5                 ### number of chromosomes
#     # ϵ = Int(1e+15)        ### some arbitrarily large number to signify the distance at which LD is nil
#     # a = 2                 ### number of alleles per locus
#     # vec_chr_lengths = [0] ### chromosome lengths
#     # vec_chr_names = [""]  ### chromosome names 
#     # dist_noLD = 500_000   ### distance at which LD is nil (related to ϵ)
#     # o = 1_000               ### total number of simulated individuals
#     # t = 10                ### number of random mating constant population size generation to simulate
#     # nQTL = 10             ### number of QTL to simulate
#     # heritability = 0.5    ### narrow(broad)-sense heritability as only additive effects are simulated
#     # LD_chr = ""           ### chromosome to calculate LD decay from
#     # LD_n_pairs = 10_000   ### number of randomly sampled pairs of loci to calculate LD
#     # plot_LD = false       ### simulate# plot simulated LD decay
#     # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
#     # npools = 100
#     # @time X, y = POOL(_X, _y, npools)
#     syncx = "../test/test_Lr.syncx"
#     phenotype = "../test/test_Lr.csv"
#     maf = 0.001
#     delimiter = ","
#     header = true
#     id_col = 1
#     phenotype_col = 2
#     missing_strings = ["NA", "NAN", "NaN", "missing", ""]
#     χ = LOAD(syncx, true)
#     ### Load phenotype data
#     ϕ = LOAD(phenotype, delimiter, header, id_col, [phenotype_col], missing_strings)
#     ### Remove missing observations
#     idx = .!ismissing.(ϕ.phe[:,1])
#     χ.cou = χ.cou[:, idx]
#     ϕ = Phenotype(ϕ.iid[idx], ϕ.tid, ϕ.phe[idx, 1:1])
#     X, vec_chr, vec_pos, vec_ale = FILTER(χ, maf)
#     y = ϕ.phe[:, 1]


#     using UnicodePlots
#     BLAS.set_num_threads(length(Sys.cpu_info())-1)

#     # features = X
#     # labels = y
#     # model = DecisionTree.build_forest(labels, features,
#     #                         n_subfeatures,
#     #                         n_trees)
#     # n_folds=10
#     # r2 =  DecisionTree.nfoldCV_forest(labels, features,
#     #                      n_folds, 
#     #                      n_subfeatures,
#     #                      n_trees)
#     #                     #  partial_sampling,
#     #                     #  max_depth,
#     #                     #  min_samples_leaf,
#     #                     #  min_samples_split,
#     #                     #  min_purity_increase;
#     #                     #  verbose = true,
#     #                     #  rng = seed)
#     # UnicodePlots.histogram(r2)
#     # mean(r2)


#     # n_subfeatures=-1 ### (default: -1 for sqrt(# features))
#     n_subfeatures=-1
#     n_trees=1_000
#     partial_sampling=0.7
#     max_depth=-1 ### (default: -1 for no maximum)
#     # min_samples_leaf=5
#     min_samples_leaf=10
#     # min_samples_split=2
#     min_samples_split=10
#     min_purity_increase=0.0
#     seed=3

#     n, p = size(X)
#     vec_idx = collect(1:n)
#     vec_y_true = []
#     vec_y_pred = []
#     @showprogress for i in 1:n
#         # i = 1
#         idx = vec_idx[vec_idx .!= i]
#         labels = y[idx]
#         features = X[idx, :]
#         model = DecisionTree.build_forest(labels, features,
#                             n_subfeatures,
#                             n_trees)
#                             #  partial_sampling,
#                             #  max_depth,
#                             #  min_samples_leaf,
#                             #  min_samples_split,
#                             #  min_purity_increase;
#                             #  rng = seed)
#         y_true = y[i]
#         y_pred = DecisionTree.apply_forest(model, X[i, :])
#         append!(vec_y_true, y_true)
#         append!(vec_y_pred, y_pred)
#     end
#     correlation_pearson, correlation_spearman, correlation_kendall, R2, R2_adj, MAE, MBE, RAE, MSE, RMSE, RRMSE, RMSLE, SMAE, SMBE, SRMSE, p = CV_METRICS(Float64.(vec_y_true), Float64.(vec_y_pred))



#     return(β̂)
# end
