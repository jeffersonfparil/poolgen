### GENOMIC PREDICTION

# include("structs.jl")
# include("functions_io.jl")
# include("functions_filterTransform.jl")
# include("functions_linearModel.jl")
# using .structs: PileupLine, SyncxLine, LocusAlleleCounts, Window, PhenotypeLine, Phenotype, MinimisationError

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
    # correlation_pearson, correlation_spearman, correlation_kendall, R2, R2_adj, MAE, MBE, RAE, MSE, RMSE, RRMSE, RMSLE, p = CV_METRICS(y, ŷ)
    min_y = [minimum(y)<0 ? abs(minimum(y)) : 0][1] ### to get rid of negative values for 
    min_ŷ = [minimum(ŷ)<0 ? abs(minimum(ŷ)) : 0][1] ### to get rid of negative values for 
    max_y = [maximum(y)<0 ? abs(maximum(y)) : 0][1] ### to get rid of negative values for 
    max_ŷ = [maximum(ŷ)<0 ? abs(maximum(ŷ)) : 0][1] ### to get rid of negative values for 
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

function CV_MULTIVAR(nrep::Int64, nfold::Int64, syncx::String, maf::Float64, phenotype::String, delimiter::String, header::Bool=true, id_col::Int=1, phenotype_col::Int=2, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], filter_genotype::Bool=true, transform_phenotype::Bool=true, standardise::Bool=false, model::Function=OLS, params=["N<<P"], save_plots::Bool=false, save_predictions::Bool=false, save_summary_plot::Bool=false, out::String="")::String
    # # n = 5                 ### number of founders
    # # m = 10_000            ### number of loci
    # # l = 135_000_000       ### total genome length
    # # k = 5                 ### number of chromosomes
    # # ϵ = Int(1e+15)        ### some arbitrarily large number to signify the distance at which LD is nil
    # # a = 2                 ### number of alleles per locus
    # # vec_chr_lengths = [0] ### chromosome lengths
    # # vec_chr_names = [""]  ### chromosome names 
    # # dist_noLD = 500_000   ### distance at which LD is nil (related to ϵ)
    # # o = 1_000               ### total number of simulated individuals
    # # t = 10                ### number of random mating constant population size generation to simulate
    # # nQTL = 10             ### number of QTL to simulate
    # # heritability = 0.2    ### narrow(broad)-sense heritability as only additive effects are simulated
    # # LD_chr = ""           ### chromosome to calculate LD decay from
    # # LD_n_pairs = 10_000   ### number of randomly sampled pairs of loci to calculate LD
    # # plot_LD = false       ### simulate# plot simulated LD decay
    # # @time vec_chr, vec_pos, _X, _y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    # # npools = 50
    # # @time X, y = POOL(_X, _y, npools)
    # # ### Skew y
    # # using UnicodePlots
    # # y = (y .- minimum(y)) ./ (maximum(y)-minimum(y))
    # # UnicodePlots.histogram(Float64.(y))
    # # y = y.^7
    # # UnicodePlots.histogram(Float64.(y))
    # # @time syncx, phenotype = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X, y)
    # syncx = "../test/test_Lr.syncx"
    # phenotype = "../test/test_Lr.csv"
    # # ###########################
    # nfold = 10
    # nrep = 3
    # maf = 0.001
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # _model = ["GBLUP", "RRBLUP"][2]
    # _method = ["ML", "REML"][1]
    # _FE_method = ["CANONICAL", "N<<P"][2]
    # _inner_optimizer=["GradientDescent", "LBFGS", "BFGS", "SimulatedAnnealing", "NelderMead"][1]
    # _optim_trace = false
    # _GBLUP_K = ["", "XTX", "COR"][2]
    # # model = OLS; params = [_FE_method]
    # # model = GLMNET; params = [1.00]
    # model = MM; params = [_model, _method, _inner_optimizer, _optim_trace, _FE_method, _GBLUP_K]
    # filter_genotype = true
    # transform_phenotype = false
    # standardise = false
    # save_plots = false
    # save_predictions = false
    # save_summary_plot = false
    # out = ""
    # CV_MULTIVAR(nfold, nrep, syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, filter_genotype, transform_phenotype, standardise, model, params, save_plots, save_predictions, save_summary_plot, out)
    ### Output tab-delimeted file including the (1) chromosome name, (2) position, (3) allele, (4) allele frequency, (5) allele effect, an (6) p-value
    if out==""
        if string(model) == "LMM_MULTIVAR"
            out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-", string(model), "_", params[2], "_CV.tsv")
        else    
            out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-", string(model), "_CV.tsv")
        end
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
        y, vec_fun, vec_par = TRANSFORM(y, true, true, true)
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
    # BLAS.set_num_threads(1)
    # vec_rep = []
    # vec_fold = []
    # vec_correlation_pearson = []
    # vec_correlation_spearman = []
    # vec_correlation_kendall = []
    # vec_R2 = []
    # vec_R2_adj = []
    # vec_MAE = []
    # vec_MBE = []
    # vec_RAE = []
    # vec_MSE = []
    # vec_RMSE = []
    # vec_RRMSE = []
    # vec_RMSLE = []
    # vec_SMAE = []
    # vec_SMBE = []
    # vec_SRMSE = []
    # vec_y_true = []
    # vec_y_pred = []
    # @showprogress "Cross-validating: " for i in 1:(nrep*nfold)
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
        # append!(vec_rep, rep)
        # append!(vec_fold, fold)
        # append!(vec_correlation_pearson, correlation_pearson)
        # append!(vec_correlation_spearman, correlation_spearman)
        # append!(vec_correlation_kendall, correlation_kendall)
        # append!(vec_R2, R2)
        # append!(vec_R2_adj, R2_adj)
        # append!(vec_MAE, MAE)
        # append!(vec_MBE, MBE)
        # append!(vec_RAE, RAE)
        # append!(vec_MSE, MSE)
        # append!(vec_RMSE, RMSE)
        # append!(vec_RRMSE, RRMSE)
        # append!(vec_RMSLE, RMSLE)
        # append!(vec_SMAE, SMAE)
        # append!(vec_SMBE, SMBE)
        # append!(vec_SRMSE, SRMSE)
        # append!(vec_y_true, y_validate)
        # append!(vec_y_pred, ŷ_validate)
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
