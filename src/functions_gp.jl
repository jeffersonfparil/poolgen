### GENOMIC PREDICTION

####### TEST ########
# include("structs.jl")
# using .structs: PileupLine, SyncxLine, LocusAlleleCounts, Window, PhenotypeLine, Phenotype, MinimisationError
# include("functions_io.jl")
# using ProgressMeter, Distributions
# include("functions_filterTransform.jl")
# using LinearAlgebra, MultivariateStats, GLMNet, Optim
# include("functions_linearModel.jl")
# using StatsBase, Plots, Distributed, Dates
#####################

function FIT(syncx::String, maf::Float64, phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_col::Int=2, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], filter_genotype::Bool=true, transform_phenotype::Bool=true, standardise::Bool=false, model::Function=[OLS, GLMNET, MM][1], params=[["N<<P"], [0.5], ["RRBLUP", "ML", "GradientDescent", true, "N<<P", "XTX"]][1], out::String="")::String
    ####### TEST ########
    # syncx = "../test/test.syncx"
    # maf = 0.01
    # phenotype = "../test/test.csv"
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # filter_genotype = true
    # transform_phenotype = true
    # standardise = false
    # ### OLS
    # model = OLS
    # params = ["N<<P"]
    # out = ""
    # FIT(syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, filter_genotype, transform_phenotype, standardise, model, params, out)
    # ### GLMNET
    # model = GLMNET
    # params = [0.51]
    # out = ""
    # FIT(syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, filter_genotype, transform_phenotype, standardise, model, params, out)
    # ### MM
    # model = MM
    # MM_model = ["GBLUP", "ABLUP", "RRBLUP"][2]
    # MM_method = ["ML", "REML"][1]
    # MM_inner_optimizer = ["GradientDescent", "LBFGS", "BFGS", "SimulatedAnnealing", "NelderMead"][1]
    # MM_optim_trace = false
    # FE_method = ["CANONICAL", "N<<P"][2]
    # GBLUP_K = ["XTX", "COR"][2]
    # params = [MM_model, MM_method, MM_inner_optimizer, MM_optim_trace, FE_method, GBLUP_K]
    # out = ""
    # FIT(syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, filter_genotype, transform_phenotype, standardise, model, params, out)
    #####################
    if out==""
        if string(model) == "MM"
            out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-", string(model), "_", join(params[1:2], "_"), "_FIT.tsv")
        elseif string(model) == "GLMNET"
            out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-", string(model), "_", params[1], "_FIT.tsv")
        else
            out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-", string(model), "_FIT.tsv")
        end
    end
    if isfile(out)
        out_basename = join(split(out, ".")[1:(end-1)], ".")
        out_extension = split(out, ".")[end]
        out = string(out_basename, "-", Dates.now(Dates.UTC), ".", out_extension)
    end
    ### Load enotype data
    χ = LOAD(syncx, true)
    ### Load phenotype data
    ϕ = LOAD(phenotype, delimiter, header, id_col, [phenotype_col], missing_strings)
    ### Remove missing observations
    idx = .!ismissing.(ϕ.phe[:,1])
    χ.cou = χ.cou[:, idx]
    ϕ = Phenotype(ϕ.iid[idx], ϕ.tid, ϕ.phe[idx, 1:1])
    ### Filter genotype data
    if filter_genotype
        X, vec_chr, vec_pos, vec_ale = FILTER(χ, maf)
    else
        X = Float64.(χ.cou')
        vec_chr = repeat(χ.chr, inner=7)
        vec_pos = repeat(χ.pos, inner=7)
        vec_ale = repeat(["A", "T", "C", "G", "INS", "DEL", "N"], outer=length(χ.chr))
    end
    n, p = size(X)
    ### Transform phenotype data
    y = Float64.(ϕ.phe[:,1])
    if transform_phenotype
        y, vec_fun_skew, vec_max_skew = TRANSFORM(y)
    else
        vec_fun = []
        vec_max = []
    end
    ### Standardise X and y so we can get rid of the intercept, i.e. for one less parameter to estimate
    if standardise
        μX = mean(X)
        σX = std(X)
        μy = mean(y)
        σy = std(y)
        X = (X .- μX) ./ σX
        y = (y .- μy) ./ σy
    end
    ### Check is we have the same number of individuals in the genotype and phenotype data
    @assert n == length(y) "Genotype and phenotype data mismatch!"
    ### Fit model
    β̂ = model(X, y, params...)
    ### Save (Note: p-values are all set to 1.0)
    file_out = open(out, "a")
    line = join(["Intercept",
                    0,
                    "NA",
                    1.0,
                    β̂[1],
                    1.0], "\t")
    write(file_out, string(line, "\n"))
    for i in 1:(length(β̂)-1)
        line = join([vec_chr[i],
                     vec_pos[i],
                     vec_ale[i],
                     mean(X[:,i]),
                     β̂[i+1],
                     1.0], "\t")
        write(file_out, string(line, "\n"))
    end
    close(file_out)
    return(out)
end

function PREDICT(tsv::String, syncx_validate::String)::Vector{Float64}
    ####### TEST ########
    # syncx = "../test/test.syncx"
    # maf = 0.01
    # phenotype = "../test/test.csv"
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # ### OLS
    # model = OLS
    # params = ["N<<P"]
    # tsv = FIT(syncx, maf, phenotype, delimiter, header, id_col, phenotype_col)
    # syncx_validate = syncx
    # PREDICT(tsv, syncx)
    #####################
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
        # @show i
        idx_chr = χ.chr[i] .== vec_chr
        idx_pos = χ.pos[i] .== vec_pos
        if (sum(idx_chr)>0) & (sum(idx_pos)>0)
            idx_ale = [sum(vec_ale[idx_chr .& idx_pos] .== x)>0 for x in vec_alleles]
            append!(idx, collect(((7*(i-1))+1):(7*i))[idx_ale])
        end            
    end
    ### It is required that we capture all the alleles and loci in the estimated effects file (tsv) are present in the validation genotype data (syncx_validate)
    if length(idx) < length(vec_beta) - 1 ### Account for the intercept
        println("ERROR: Missing loci in prediction dataset.")
        println("Please, consider havng the same loci covered in both training and prediction sets, i.e.")
        println(string("in '", tsv, "' and '", syncx_validate, "' files."))
        return(1)
    end
    X = χ.cou'[:, idx]
    n, _ = size(X)
    ŷ = hcat(ones(n), X) * vec_beta
    return(ŷ)
end

function CV_METRICS(y::Vector{T}, ŷ::Vector{T})::Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Plots.Plot{Plots.GRBackend}} where T <: Number
    ###### TEST ########
    # syncx = "../test/test.syncx"
    # maf = 0.01
    # phenotype = "../test/test.csv"
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # filter_genotype = true
    # transform_phenotype = false
    # standardise = false
    # ### OLS
    # model = OLS
    # params = ["N<<P"]
    # tsv = FIT(syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, filter_genotype, transform_phenotype, standardise)
    # syncx_validate = syncx
    # ϕ = LOAD(phenotype, delimiter, header, id_col, [phenotype_col], missing_strings)
    # idx = .!ismissing.(ϕ.phe[:,1])
    # y = Float64.(ϕ.phe[idx,1])
    # ŷ = PREDICT(tsv, syncx)
    # CV_METRICS(y, ŷ)
    ####################
    min_y = [minimum(y)<0 ? abs(minimum(y)) : 0][1] ### to get rid of negative values for 
    min_ŷ = [minimum(ŷ)<0 ? abs(minimum(ŷ)) : 0][1] ### to get rid of negative values for 
    max_y = [maximum(y)<0 ? abs(maximum(y)) : 0][1] ### to get rid of negative values for 
    max_ŷ = [maximum(ŷ)<0 ? abs(maximum(ŷ)) : 0][1] ### to get rid of negative values for 
    ### Correlation metrics
    correlation_pearson = StatsBase.cor(y, ŷ)
    correlation_spearman = StatsBase.corspearman(y, ŷ)
    correlation_kendall = StatsBase.corkendall(y, ŷ)
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
    ### Scaled mean absolute error
    SMAE = mean(abs.(y .- ŷ)) / max_y
    ### Scaled mean bias error
    SMBE = mean(y .- ŷ) / max_y
    ### Scaled root mean square error
    SRMSE = sqrt(MSE) / max_y
    ### Plot
    x = vcat(y, ŷ)
    δ = std(x) / sqrt(length(x))
    limits = [minimum(x)-δ, maximum(x)+δ]
    p = Plots.scatter(y, ŷ, xlabel="True", ylab="Predicted", xlims=limits, ylims=limits,
                      markercolor=RGB(0.5,0.5,0.9), markerstrokewidth=0.001, markeralpha=0.4, title="", legend=false);
    Plots.plot!(p, [0,1], [0 ,1], seriestype=:straightline, linecolor=:gray, legend=false);
    Plots.annotate!(p, limits[1]+(2*δ), limits[2],
                text(string("R²=", round(R2, digits=2),
                            "\nRMSE=", round(RMSE, digits=2),
                            "\nRRMSE=", round(RRMSE, digits=2),
                            "\nMAE=", round(MAE, digits=2),
                            "\nMBE=", round(MBE, digits=2),
                            "\nRAE=", round(RAE, digits=2)),
                     :left, :top, 10)
                   );
    ### Output metrics and plot
    return(correlation_pearson, correlation_spearman, correlation_kendall,
           R2, R2_adj,
           MAE, MBE, RAE,
           MSE, RMSE, RRMSE, RMSLE,
           SMAE, SMBE, SRMSE,
           p)
end

function CV_MULTIVAR(nrep::Int64, nfold::Int64, syncx::String, maf::Float64, phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_col::Int=2, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], filter_genotype::Bool=true, transform_phenotype::Bool=true, standardise::Bool=false, model::Function=[OLS, GLMNET, MM][1], params=[["N<<P"], [0.5], ["RRBLUP", "ML", "GradientDescent", true, "N<<P", "XTX"]][1], save_plots::Bool=false, save_predictions::Bool=false, save_summary_plot::Bool=false, out::String="")::String
    # ####### TEST ########
    # nrep = 3
    # nfold = 10
    # syncx = "../test/test_Lr.syncx"
    # maf = 0.01
    # phenotype = "../test/test_Lr.csv"
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # filter_genotype = true
    # transform_phenotype = true
    # standardise = false
    # save_plots = true
    # save_predictions = true
    # save_summary_plot = true
    # ### OLS
    # model = OLS
    # params = ["N<<P"]
    # out = ""
    # CV_MULTIVAR(nrep, nfold, syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, filter_genotype, transform_phenotype, standardise, model, params, save_plots, save_predictions, save_summary_plot, out)
    # ### MM
    # model = MM
    # MM_model = ["GBLUP", "RRBLUP"][2]
    # MM_method = ["ML", "REML"][1]
    # MM_inner_optimizer = ["GradientDescent", "LBFGS", "BFGS", "SimulatedAnnealing", "NelderMead"][1]
    # MM_optim_trace = false
    # FE_method = ["CANONICAL", "N<<P"][2]
    # GBLUP_K = ["XTX", "COR"][2]
    # params = [MM_model, MM_method, MM_inner_optimizer, MM_optim_trace, FE_method, GBLUP_K]
    # out = ""
    # CV_MULTIVAR(nrep, nfold, syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, filter_genotype, transform_phenotype, standardise, model, params, save_plots, save_predictions, save_summary_plot, out)
    #####################
    ### Output tab-delimeted file including the (1) chromosome name, (2) position, (3) allele, (4) allele frequency, (5) allele effect, an (6) p-value
    if out==""
        if string(model) == "MM"
            out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-", string(model), "_", params[2], "_CV.tsv")
        else    
            out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-", string(model), "_CV.tsv")
        end
    end
    if isfile(out)
        out_basename = join(split(out, ".")[1:(end-1)], ".")
        out_extension = split(out, ".")[end]
        out = string(out_basename, "-", Dates.now(Dates.UTC), ".", out_extension)
    end
    ### Load genotype data
    χ = LOAD(syncx, true)
    ### Load phenotype data
    ϕ = LOAD(phenotype, delimiter, header, id_col, [phenotype_col], missing_strings)
    ### Remove missing observations
    idx = .!ismissing.(ϕ.phe[:,1])
    χ.cou = χ.cou[:, idx]
    ϕ = Phenotype(ϕ.iid[idx], ϕ.tid, ϕ.phe[idx, 1:1])
    ### Filter genotype data
    if filter_genotype
        X, vec_chr, vec_pos, vec_ale = FILTER(χ, maf)
    else
        X = Float64.(χ.cou')
        vec_chr = repeat(χ.chr, inner=7)
        vec_pos = repeat(χ.pos, inner=7)
        vec_ale = repeat(["A", "T", "C", "G", "INS", "DEL", "N"], outer=length(χ.chr))
    end
    n, p = size(X)
    ### Transform phenotype data
    y = Float64.(ϕ.phe[:,1])
    if transform_phenotype
        y, vec_fun, vec_max = TRANSFORM(y)
    else
        vec_fun = []
        vec_max = []
    end
    ### Standardise X and y so we can get rid of the intercept, i.e. for one less parameter to estimate
    if standardise
        μX = mean(X)
        σX = std(X)
        μy = mean(y)
        σy = std(y)
        X = (X .- μX) ./ σX
        y = (y .- μy) ./ σy
    end
    ### Save predictions if we waht to save the final plot
    if save_summary_plot
        save_predictions = true
    end
    ### Check is we have the same number of individuals in the genotype and phenotype data
    @assert n == length(y) "Genotype and phenotype data mismatch!"
    ### Divide the observations into nfold partitions
    if (n/nfold) < 1
        println("Sorry you dataset is too small. Exiting function now.")
        return("Error")
    end
    vec_fld = repeat(collect(1:nfold), inner=Int(floor(n/nfold)))
    if length(vec_fld) < n
        append!(vec_fld, repeat([nfold], n-length(vec_fld)))
    end
    ### Generate nrep random permutations of the partitioning of the observations
    mat_idx_rand = zeros(Int, n, nrep)
    for i in 1:nrep
        # i = 1
        mat_idx_rand[:, i] = sample(vec_fld, n, replace=false)
    end
    ### Cross-validate for nrep randomisations and k-fold cross-validation
    ### Parallel execution
    vec_vec_metrics = @sync @showprogress "Cross-validation: " @distributed (hcat) for i in 1:(nrep*nfold)
        # i = 1
        ### Determine the rep and fold numbers
        rep = Int(ceil(i/nfold))
        fold = mod(i, nfold)
        if fold == 0
            fold = nfold
        end
        ### Set the pre-randomised training and validation sets
        idx_training = mat_idx_rand[:, rep] .!= fold
        idx_validate = mat_idx_rand[:, rep] .== fold
        
        X_training = X[idx_training, :]
        X_validate = X[idx_validate, :]
        y_training = y[idx_training]
        y_validate = y[idx_validate]

        β̂ = model(X_training, y_training, params...)

        ŷ_validate = hcat(ones(length(y_validate)), X_validate) * β̂
        if transform_phenotype
            y_validate = UNTRANSFORM(y_validate, vec_fun, vec_max)
            ŷ_validate = UNTRANSFORM(ŷ_validate, vec_fun, vec_max)
        end
        if standardise
            y_validate = (y_validate .* σy) .+ μy
            ŷ_validate = (ŷ_validate .* σy) .+ μy
        end
        correlation_pearson, correlation_spearman, correlation_kendall, R2, R2_adj, MAE, MBE, RAE, MSE, RMSE, RRMSE, RMSLE, SMAE, SMBE, SRMSE, p = CV_METRICS(y_validate, ŷ_validate)
        ### Vectorise for hcat-ing
        [rep, fold, correlation_pearson, correlation_spearman, correlation_kendall, R2, R2_adj, MAE, MBE, RAE, MSE, RMSE, RRMSE, RMSLE, SMAE, SMBE, SRMSE, join(round.(y_validate, digits=2), ","), join(round.(ŷ_validate, digits=2), ",")]
    end
    ### Consolidate metrics
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
    vec_SMAE = []
    vec_SMBE = []
    vec_SRMSE = []
    if save_predictions
        vec_true = []
        vec_pred = []
    end
    t = size(vec_vec_metrics,2)
    for j in 1:t
        rep, fold, correlation_pearson, correlation_spearman, correlation_kendall, R2, R2_adj, MAE, MBE, RAE, MSE, RMSE, RRMSE, RMSLE, SMAE, SMBE, SRMSE, y_true, y_pred = vec_vec_metrics[:, j]
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
        append!(vec_SMAE, SMAE)
        append!(vec_SMBE, SMBE)
        append!(vec_SRMSE, SRMSE)
        if save_predictions
            push!(vec_true, y_true)
            push!(vec_pred, y_pred)
        end
    end
    ### Output
    if save_predictions
        out_pred = string(out, "-predictions.tsv")
        f = open(out_pred, "a")
        write(f, string(join(["rep", "fold", "true", "pred"], "\t"), "\n"))
        for i in 1:t
            y_true = parse.(Float64, split(vec_true[i], ","))
            y_pred = parse.(Float64, split(vec_pred[i], ","))
            for j in 1:length(y_true)
                write(f, string(join([vec_rep[i], vec_fold[i], y_true[j], y_pred[j]], "\t"), "\n"))
            end
        end
        close(f)
    end
    if save_summary_plot
        vec_y = []; vec_ŷ = []
        f = open(out_pred, "r")
        _ = readline(f)
        while !eof(f)
            line = split(readline(f), "\t")
            push!(vec_y, parse(Float64, line[3]))
            push!(vec_ŷ, parse(Float64, line[4]))
        end
        close(f)
        _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, p = CV_METRICS(Float64.(vec_y), Float64.(vec_ŷ))
        Plots.savefig(p, string(out, "-summary.svg"));
    end
    file_out = open(out, "a")
    line = string(join(["rep", "fold", "correlation_pearson", "correlation_spearman", "correlation_kendall", "R2", "R2_adj", "MAE", "MBE", "RAE", "MSE", "RMSE", "RRMSE", "RMSLE"], "\t"), "\n")
    write(file_out, line)
    for i in 1:length(vec_rep)
        line = string(join([Int(vec_rep[i]), Int(vec_fold[i])], "\t"), "\t", join([vec_correlation_pearson[i], vec_correlation_spearman[i], vec_correlation_kendall[i], vec_R2[i], vec_R2_adj[i], vec_MAE[i], vec_MBE[i], vec_RAE[i], vec_MSE[i], vec_RMSE[i], vec_RRMSE[i], vec_RMSLE[i]], "\t"), "\n")
        write(file_out, line)
    end
    close(file_out)
    return(out)
end


# ### Test the idea of variance estimation using GBLUP and estimate effects using RRBLUP equations:

# n = 100
# m = 1_000
# a = 10
# X = rand(n, m)
# b = zeros(m)
# b_loc = sort(sample(collect(1:m), a))
# b[b_loc] = rand(Distributions.Normal(0, 2.5), a)
# y = X * b + rand(Distributions.Normal(0, 0.5), n)

# ### Estimate variances and use the  σ2u and σ2e estimates to calculate the BLUEs and BLUPs
# function test(X, y)
#     n, p = size(X)
#     K = GENERATE_COVARIATE(X, n, "COR")
#     σ2u, σ2e = OPTIM_MM(ones(n),
#                         y,
#                         1.0*I,
#                         K,
#                         "N<<P",
#                         "ML")

#     ### Random effects variance-covariance matrix
#     D = (σ2u/p)*I
#     ### Error variance-covariance matrix (homoscedastic)
#     R = σ2e*I
#     ### Solve the mixed model equations
#     Z = X
#     X = ones(n)
#     β̂, μ̂ = MM(X, y, Z, D, R, "N<<P", true)
#     out = vcat(β̂[1], μ̂)
#     return(out)
# end

# @time grrblup = test(X, y);
# @time rrblup = MM(X, y, "RRBLUP");
# @time ablup = MM(X, y, "ABLUP");
# @time gblup = MM(X, y, "GBLUP");

# using UnicodePlots
# UnicodePlots.histogram(grrblup[2:end])
# UnicodePlots.histogram(rrblup[2:end])
# UnicodePlots.histogram(ablup[2:end])
# UnicodePlots.histogram(gblup[2:end])

# UnicodePlots.scatterplot(grrblup[2:end], rrblup[2:end])
# UnicodePlots.scatterplot(grrblup[2:end], ablup[2:end])
# UnicodePlots.scatterplot(grrblup[2:end], gblup[2:end])

# cor(grrblup[2:end], rrblup[2:end])
# cor(grrblup[2:end], ablup[2:end])
# cor(grrblup[2:end], gblup[2:end])

# UnicodePlots.histogram(y - hcat(ones(n), X) * grrblup)
# UnicodePlots.histogram(y - hcat(ones(n), X) * rrblup)
# UnicodePlots.histogram(y - hcat(ones(n), X) * ablup)
# UnicodePlots.histogram(y - hcat(ones(n), X) * gblup)

# sqrt(mean((y - hcat(ones(n), X) * grrblup).^2))
# sqrt(mean((y - hcat(ones(n), X) * rrblup).^2))
# sqrt(mean((y - hcat(ones(n), X) * ablup).^2))
# sqrt(mean((y - hcat(ones(n), X) * gblup).^2))