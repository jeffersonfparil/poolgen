### IMPUTATION

####### TEST ########
# include("structs.jl")
# using .structs: PileupLine, SyncxLine, LocusAlleleCounts, Window, PhenotypeLine, Phenotype, MinimisationError
# include("functions_io.jl")
# using ProgressMeter, Distributions
# include("functions_filterTransform.jl")
# using LinearAlgebra, GLMNet, Optim
# using MultivariateStats
# include("functions_linearModel.jl")
#####################

function SIMULATE_MISSING!(window::Window, rate::Float64=0.01)::Window
    ####### TEST ########
    # syncx = "../test/test.syncx"
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
    # close(file)
    # rate = 0.01
    #####################
    p, n = size(window.cou)
    vec_idx = collect(1:p)
    s = Int(ceil(rate*p))
    idx_alleles = sample(vec_idx[mean(window.cou, dims=2)[:,1] .> 0], s)
    idx_pools = sample(repeat([1, 3], inner=10), s)
    window.cou[idx_alleles, idx_pools] .= missing
    return(window)
end

function SIMULATE_MISSING(syncx::String, init::Int, term::Int, window_size::Int=100, rate::Float64=0.01, out::String="")::String
    ####### TEST ########
    # syncx = "../test/test.syncx"
    # file = open(syncx, "r")
    # seekend(file)
    # pos = position(file)
    # threads = 4
    # vec_positions = Int.(round.(collect(0:(position(file)/threads):position(file))))
    # close(file)
    # init = vec_positions[1]
    # term = vec_positions[2]
    # window_size = 100
    # rate = 0.01
    # out = ""
    #####################
    ### Output syncx with simulated missing data
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-SIMULATED_MISSING.syncx")
    end
    out = AVOID_FILE_OVERWRITE(out)
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
    while (i < term) && (!eof(file))
        j = 0
        window = []
        while (j < window_size) && (i < term) && (!eof(file))
            j += 1
            locus = PARSE(SyncxLine(j, readline(file)))
            i = position(file)
            push!(window, locus)
        end
        window = PARSE(Array{LocusAlleleCounts}(window))
        SIMULATE_MISSING!(window, rate)
        SAVE(window, out)
    end
    close(file)
    return(out)
end

function IMPUTE!(window::Window, model::String=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance::Bool=true)::Window
    ####### TEST ########
    # syncx = "../test/test.syncx"
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
    # close(file)
    # SIMULATE_MISSING!(window)
    # model = ["Mean", "OLS", "RR", "LASSO", "GLMNET"][2]
    # distance = true
    #####################
    ### Find the indices of pools with missing data. These will be used independently and iteratively as our response variables
    p, n = size(window.cou)
    idx_pools = sum(ismissing.(window.cou), dims=1)[1,:] .> 0
    ### If we have at least one pool with no missing data, then we proceed with imputation
    if (sum(.!idx_pools) >= 1) && (sum(idx_pools) >= 1)
        ### Explanatory variables
        X = Int.(window.cou[:, .!idx_pools])
         ### Distance covariate (only add if the window is within a single chromosome)
         if (distance) && (length(unique(window.chr))==1) && (model != "Mean")
            m = length(window.pos)
            D = zeros(Int, m, m)
            @simd for i in 1:m
                for j in 1:m
                    D[i,j] = abs(window.pos[i] - window.pos[j])
                end
            end
            Z = MultivariateStats.projection(MultivariateStats.fit(PCA, repeat(D, inner=(7,1)); maxoutdim=1))
            X = hcat(X, (X .!= 0) .* Z) ### multiply Z by (X .!= 0) so we get rid of covariate effects when X_ij is zero
        end
        ### Impute using pools without missing data
        for j in collect(1:n)[idx_pools]
            # j = collect(1:n)[idx_pools][1]
            y = window.cou[:, j]
            idx_loci = ismissing.(y)
            y_train = Float64.(y[.!idx_loci])
            X_train = X[.!idx_loci, :]
            nf, pf = size(X_train)
            ### Train models
            if model == "Mean"
                β = append!([0.0], repeat([1/pf], pf))
            elseif model == "OLS"
                β = try
                    OLS(X_train, y_train, "MULTIALGORITHMIC")
                catch
                    try
                        OLS(X_train, y_train, "CANONICAL")
                    catch
                        missing
                    end
                end
            elseif (model == "RR") || (model == "LASSO") || (model == "GLMNET")
                model=="RR" ? alpha1=0 : model=="LASSO" ? alpha1=1 : alpha1=0.5
                β = try
                        GLMNET(X_train, y_train, alpha1)
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
    ####### TEST ########
    # ### Simulate missing
    # syncx_raw = "../test/test_Lr.syncx"
    # file = open(syncx_raw, "r")
    # seekend(file)
    # pos = position(file)
    # threads = 4
    # vec_positions = Int.(round.(collect(0:(position(file)/threads):position(file))))
    # close(file)
    # vec_out_filenames = []
    # for i in 2:length(vec_positions)
    #     # i = 2
    #     init = vec_positions[i-1]
    #     term = vec_positions[i]
    #     out = SIMULATE_MISSING(syncx_raw, init, term, 100, 0.01, string(join([syncx_raw, init, term], "-"), ".tmp"))
    #     push!(vec_out_filenames, out)
    # end
    # ### Merge simulated missing Pool-seq data
    # syncx = "../test/test_Lr-SIMULATED_MISSING.syncx"
    # MERGE(string.(vec_out_filenames), syncx)
    # file = open(syncx, "r")
    # seekend(file)
    # threads = 4
    # vec_positions = Int.(round.(collect(0:(position(file)/threads):position(file))))
    # close(file)
    # init = vec_positions[2]
    # term = vec_positions[3]
    # window_size = 100
    # model = "OLS"
    # distance = true
    # out = ""
    #####################
    ### Output syncx
    if out==""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-IMPUTED.syncx")
    end
    out = AVOID_FILE_OVERWRITE(out)
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
    while (i < term) && (!eof(file))
        if window == []
            while (j < window_size) && (i < term) && (!eof(file))
                j += 1
                locus = PARSE(SyncxLine(j, readline(file)))
                i = position(file)
                push!(window, locus)
            end
            window = PARSE(Array{LocusAlleleCounts}(window))
            IMPUTE!(window, model, distance)
            SAVE(EXTRACT(window, 1), out) ### Save the first locus
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
            SLIDE!(window, locus) ### Slide one locus at-a-time
            SAVE(EXTRACT(window, 1), out) ### Save the leading locus
            IMPUTE!(window, model, distance)
        end
    end
    if !eof(file)
        SAVE(EXTRACT(window, 2:window_size), out) ### Save the remaining trailing loci
    else
        SAVE(window, out)
    end
    close(file)
    return(out)
end
