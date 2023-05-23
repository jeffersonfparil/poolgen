### GP
library(glmnet)

PERF = function(y_test, y_hat, max_y) {
    n = length(y_test)
    if ((var(y_hat)==0) | (var(y_test)==0)) {
        cor = 0.0
    } else {
        cor = cor(y_test, y_hat)
    }
    e = y_test - y_hat
    mbe = ( sum(e)/n ) / max_y
    mae = ( sum(abs(e))/n ) / max_y
    mse = ( sum(e^2)/n ) / max_y
    rmse = ( sqrt(sum(e^2)/n) ) / max_y
    return(list(cor=cor, mbe=mbe, mae=mae, mse=mse, rmse=rmse, y_hat=y_hat, y_test=y_test))
}

OLS = function(x_train, x_test, y_train, y_test) {
    b_hat = t(x_train) %*% solve(x_train %*% t(x_train)) %*% y_train
    y_hat = (x_test %*% b_hat)[,1]
    y_max = max(c(y_test, y_train), na.rm=TRUE)
    # plot(abs(b_hat))
    # for (x in idx_b) {
    #     abline(v=x, col="red", lwd=2)
    # }
    return(PERF(y_test, y_hat, y_max))
}

LASSO = function(x_train, x_test, y_train, y_test) {
    mod_lasso = cv.glmnet(x=x_train, y=y_train, alpha=1.0)
    # lambda = mod_lasso$lambda.min
    # b_hat   = coef(mod_lasso, s="lambda.min")[,1]
    # b_hat   = b_hat[2:length(b_hat)]
    # plot(abs(b_hat))
    # for (x in idx_b) {
    #     abline(v=x, col="red", lwd=2)
    # }
    y_hat = predict(mod_lasso, newx=x_test, s="lambda.min")[,1]
    y_max = max(c(y_test, y_train), na.rm=TRUE)
    return(PERF(y_test, y_hat, y_max))
}

RIDGE = function(x_train, x_test, y_train, y_test) {
    mod_ridge = cv.glmnet(x=x_train, y=y_train, alpha=0.0)
    y_hat = predict(mod_ridge, newx=x_test, s="lambda.min")[,1]
    y_max = max(c(y_test, y_train), na.rm=TRUE)
    return(PERF(y_test, y_hat, y_max))
}

KFOLD_CV = function(x, y, r=5, k=10) {
    # k = 10
    ols = c()
    lasso = c()
    ridge = c()

    n = length(y)
    s = floor(n/k)
    idx = sample(c(1:n), n, replace=FALSE)

    pb = txtProgressBar(min=0, max=r*k, initial=0, style=3)
    for (rep in 1:r) {
        for (fold in 1:k) {
            # fold = 1
            i = (fold-1)*s + 1
            j = fold*s
            if (fold == k) {
                j = n
            }
            bool_test = c(1:n) %in% c(i:j)
            bool_train = !bool_test
            idx_train = idx[bool_train]
            idx_test = idx[bool_test]

            x_train = x[idx_train, ]
            x_test = x[idx_test, ]
            y_train = y[idx_train]
            y_test = y[idx_test]

            ols = c(ols, OLS(x_train, x_test, y_train, y_test))
            lasso = c(lasso, LASSO(x_train, x_test, y_train, y_test))
            ridge = c(ridge, RIDGE(x_train, x_test, y_train, y_test))
            setTxtProgressBar(pb, ((rep-1)*k)+fold)
        }
    }
    close(pb)
    return(list(ols=ols,
                lasso=lasso,
                ridge=ridge))
}


### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### Simulated data
### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
start_time = Sys.time()
n = 100
p = 50000
q = 10
maf = 1e-4
h2 = 0.75
X_sim = matrix(runif(n*p, min=maf, max=1-maf), nrow=n)
b = rep(0, p)
idx_b = sort(sample(c(1:p), q))
b[idx_b] = 1.0
xb = X_sim %*% b
v_xb = var(xb)
v_e = (v_xb/h2) - v_xb
e = rnorm(n, mean=0, sd=sqrt(v_e))
y = xb + e
y_sim = scale(y, center=T, scale=T)


kfold = c()
model = c()
ols = c()
lasso = c()
ridge = c()
plot_model = c()
plot_y_hat = c()
plot_y_test = c()

y = y_sim
x = X_sim

k = 10
r = 5
kfold_out = KFOLD_CV(x, y, r, k)

idx_steps = seq(from=1, to=(r*k*7), by=7)
for (mod in names(kfold_out)) {
    # mod = names(kfold_out)[1]
    p_y_hat = c()
    p_y_test = c()
    for (i_ in idx_steps) {
        # i_ = 1
        eval(parse(text=paste0(mod, " = c(", mod, ", kfold_out$", mod, "[i_:(i_+4)])")))
        p_y_hat = c(p_y_hat, unlist(eval(parse(text=paste0("kfold_out$", mod, "[i_+5]")))))
        p_y_test = c(p_y_test, unlist(eval(parse(text=paste0("kfold_out$", mod, "[i_+6]")))))
    }
    plot_model = c(plot_model, rep(mod, length(p_y_hat)))
    plot_y_hat = c(plot_y_hat, p_y_hat)
    plot_y_test = c(plot_y_test, p_y_test)
}

plot_df = data.frame(model=plot_model, y_hat=plot_y_hat, y_test=plot_y_test)
# par(mfrow=c(2,2))
# plot(x=plot_df$y_test, y=plot_df$y_hat, xlab="Observed", ylab="Predicted", pch=19, main="All models"); grid()
# legend("topright", legend=paste0("cor=", round(100*cor(x=plot_df$y_test, y=plot_df$y_hat),2), "%"))
vec_mod = unique(plot_df$model)
vec_cor = c()
vec_rmse = c()
for (mod in vec_mod) {    
    idx = plot_df$model == mod
    vec_cor = c(vec_cor, round(100*cor(x=plot_df$y_test[idx], y=plot_df$y_hat[idx]),2))
    vec_rmse = c(vec_rmse, sqrt(mean((plot_df$y_test[idx]-plot_df$y_hat[idx])^2)))
    # plot(x=plot_df$y_test[idx], y=plot_df$y_hat[idx], xlab="Observed", ylab="Predicted", pch=19, main=mod); grid()
    # legend("topright", legend=paste0("cor = ", round(100*cor(x=plot_df$y_test[idx], y=plot_df$y_hat[idx]),2), "%"))
}
out = data.frame(model=vec_mod, correlation=vec_cor, rmse=vec_rmse)
print(out)
end_time = Sys.time()
print(paste0("Time elapsed: ", end_time - start_time))


### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### Empirical data
### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

setwd("/data/weedomics/2.c_60_populations_genotyping/READ_LINKS/")
G = read.csv("Lolium_SEAU-1682289942.6635504-1682294227.0662844-allele_frequencies.csv", header=TRUE)
Y = read.csv("Lolium_SEAU.csv", header=TRUE)

p = nrow(G)
n = ncol(G) - 3

### Build the X matrix while removing non-polymorphic loci
X = rep(1, p)
pb = txtProgressBar(min=0, max=p, initial=0, style=3)
for (j in 1:p) {
    # j = 1
    xj = t(G[j, 4:(n+3)])
    if (sum(is.na(xj)) > 0) {
        next
    }
    if (var(xj) > 0.001) {
        X = cbind(X, xj)
    }
    setTxtProgressBar(pb, j)
}
close(pb)

vec_herbicides = colnames(Y)[3:ncol(Y)]

herbi = c()
kfold = c()
model = c()
ols = c()
lasso = c()
ridge = c()
plot_herbi = c()
plot_model = c()
plot_y_hat = c()
plot_y_test = c()

for (i in 1: length(vec_herbicides)) {
    # i = 1
    ### Prepare genotype matrix and phenotype vector
    y = Y[, c(FALSE, FALSE, vec_herbicides==vec_herbicides[i])]
    idx = y > -10
    y = y[idx]
    x = X[idx, ]
    # y = scale(y, center=TRUE, scale=TRUE)
    y = as.matrix(y)
    x = as.matrix(x)

    k = 10
    ### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    start_time = Sys.time()
    kfold_out = KFOLD_CV(x, y, k)
    end_time = Sys.time()
    print(paste0("Time elapsed: ", end_time - start_time))
    ### @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    m = length(kfold_out)
    herbi = c(herbi, rep(vec_herbicides[i], k*m))
    kfold = c(kfold, c(1:(k*m)))
    model = c(model, rep(names(kfold_out), each=k))

    idx_steps = seq(from=1, to=(k*7), by=7)
    for (mod in names(kfold_out)) {
        # mod = names(kfold_out)[1]
        p_y_hat = c()
        p_y_test = c()
        for (i_ in idx_steps) {
            # i_ = 1
            eval(parse(text=paste0(mod, " = c(", mod, ", kfold_out$", mod, "[i_:(i_+4)])")))
            p_y_hat = c(p_y_hat, unlist(eval(parse(text=paste0("kfold_out$", mod, "[i_+5]")))))
            p_y_test = c(p_y_test, unlist(eval(parse(text=paste0("kfold_out$", mod, "[i_+6]")))))
        }
        plot_herbi = c(plot_herbi, rep(vec_herbicides[i], length(p_y_hat)))
        plot_model = c(plot_model, rep(mod, length(p_y_hat)))
        plot_y_hat = c(plot_y_hat, p_y_hat)
        plot_y_test = c(plot_y_test, p_y_test)
    }
}

df_ols = as.data.frame(matrix(unlist(ols), ncol=5, byrow=TRUE))
df_lasso = as.data.frame(matrix(unlist(lasso), ncol=5, byrow=TRUE))
df_ridge = as.data.frame(matrix(unlist(ridge), ncol=5, byrow=TRUE))
df = rbind(df_ols, df_lasso, df_ridge)
colnames(df) = c("cor", "mbe", "mae", "mse", "rmse")
df = as.data.frame(cbind(herbi, kfold, model, df))

plot_df = data.frame(herbi=plot_herbi, model=plot_model, y_hat=plot_y_hat, y_test=plot_y_test)

for (herbi in vec_herbicides) {
    # herbi = vec_herbicides[1]
    idx = plot_df$herbi == herbi
    subdf = droplevels(plot_df[idx, ])
    svg(paste0(herbi, "-scatterplots-observedXpredicted.svg"), width=10, height=7)
    par(mfrow=c(2,2))
    plot(x=subdf$y_test, y=subdf$y_hat, xlab="Observed", ylab="Predicted", pch=19, main="All models"); grid()
    legend("topright", legend=paste0("cor=", round(100*cor(x=subdf$y_test, y=subdf$y_hat),2), "%"))
    for (mod in unique(subdf$model)) {    
        idx = subdf$model == mod
        plot(x=subdf$y_test[idx], y=subdf$y_hat[idx], xlab="Observed", ylab="Predicted", pch=19, main=mod); grid()
        legend("topright", legend=paste0("cor = ", round(100*cor(x=subdf$y_test[idx], y=subdf$y_hat[idx]),2), "%"))
    }
    dev.off()
}


# ######################################
# ### LEAVE-ONE-OUT-CROSS-VALIDATION ###
# ######################################

# library(rrBLUP)
# library(glmnet)
# library(BLR)
# library(foreach)
# library(doParallel)
# library(stringr)

# syncx = "/data-weedomics-2/poolgen/test/test_Lr.syncx"
# csv = "/data-weedomics-2/poolgen/test/test_Lr.csv"
# # syncx = "/data-weedomics-2/2.d_40_populations_model_validation/Imputed_geno-mean.syncx"
# # csv = "/data-weedomics-2/2.d_40_populations_model_validation/phenotype_data.csv"
# setwd(dirname(syncx))

# # nfold =  n
# # nfold =  ceiling(n/2)
# nfold =  10

# dat = read.table(syncx, header=FALSE, sep="\t")
# ### cbind is faster than rbind
# C = matrix(as.numeric(unlist(strsplit(dat[,3], ":"))), ncol=7, byrow=TRUE)
# X = matrix(t(C / rowSums(C)), ncol=1, byrow=FALSE)
# print("Loading syncx allele counts matrix: ")
# pb = txtProgressBar(4, ncol(dat), style=3)
# for (i in 4:ncol(dat)) {
#     setTxtProgressBar(pb, i)
#     C = matrix(as.numeric(unlist(strsplit(dat[,i], ":"))), ncol=7, byrow=TRUE)
#     F = matrix(t(C / rowSums(C)), ncol=1, byrow=FALSE)
#     X = cbind(X, F)
# }
# close(pb)
# X = t(X)
# AF_mu = apply(X, MARGIN=2, FUN=mean)
# AF_sd = apply(X, MARGIN=2, FUN=sd)
# idx = (AF_mu > 0) & (AF_sd > 0)
# X = X[, idx]



# Y = read.csv(csv, header=TRUE)
# p = Y[, ncol(Y)]
# # y = log10(y +1) ### Predictions are worse!

# n = nrow(X)
# vec_fold_id = rep(c(1:nfold), each=floor(n/nfold))
# if (length(vec_fold_id) < n) {
#     vec_fold_id = c(vec_fold_id, rep(nfold, times=n-length(vec_fold_id)))
# }


# RRBLUP = function(y, SNP, SNP_test, method="ML") {
#     # i = 1
#     # y = p[vec_fold_id!=i]
#     # SNP = X[vec_fold_id!=i, ]
#     # SNP_test = X[vec_fold_id==i, ]
#     # method = "ML"
#     mod = mixed.solve(y=y, Z=SNP, method=method)
#     y_pred = mod$beta[1] + (SNP_test %*% mod$u)
#     return(y_pred)
# }

# GLMNET = function(y, SNP, SNP_test, alpha=1.00) {
#     # i = 5
#     # y = p[vec_fold_id!=i]
#     # SNP = X[vec_fold_id!=i, ]
#     # SNP_test = X[vec_fold_id==i, ]
#     # alpha = 1.00
#     mod = cv.glmnet(y=y, x=SNP, alpha=alpha)
#     idx = tryCatch(
#         idx = which(mod$lambda == mod$lambda.min),
#         error = function(e){
#             idx = which(mod$cvm == min(mod$cvm, na.rm=TRUE))
#         }
#     )
#     if (length(idx) == 0){
#         y_pred = NA
#     } else {
#         y_pred = mod$glmnet.fit$a0[idx] + (SNP_test %*% mod$glmnet.fit$beta[, idx])
#     }
#     return(y_pred)
# }

# BAYESREG = function(y, SNP, SNP_test, shrinkage=c("BayesU", "BayesR", "BayesL")[1], nIter=1e4, burnIn=1e2, thin=1e1, id=1) {
#     # i = 5
#     # y = p[vec_fold_id!=i]
#     # SNP = X[vec_fold_id!=i, ]
#     # SNP_test = X[vec_fold_id==i, ]
#     # nIter=1e3
#     # burnIn=1e2
#     # thin=1e1
#     if (shrinkage == "BayesU") {
#         ### Unshruked estimates
#         prior = list(varE=list(df=1,S=1.0),
#                      varU=list(df=1,S=1.0), 
#                      varBR=list(df=1,S=1.0))
#         mod = BLR(y=y, XF=SNP, prior=prior, nIter=nIter, burnIn=burnIn, thin=thin, saveAt=paste0("temp_", id, "-"))
#         file.remove(paste0("temp_", id, "-bF.dat"))
#         file.remove(paste0("temp_", id, "-varE.dat"))
#         y_pred = mod$mu[1] + (SNP_test %*% mod$bF)
#     } else if (shrinkage == "BayesR") {
#         ### Ridge shrinkage of estimates
#         prior = list(varE=list(df=1,S=1.0),
#                      varU=list(df=1,S=1.0), 
#                      varBR=list(df=1,S=1.0))
#         mod = BLR(y=y, XR=SNP, prior=prior, nIter=nIter, burnIn=burnIn, thin=thin, saveAt=paste0("temp_", id, "-"))
#         file.remove(paste0("temp_", id, "-varBR.dat"))
#         file.remove(paste0("temp_", id, "-varE.dat"))
#         y_pred = mod$mu[1] + (SNP_test %*% mod$bR)
#     } else if (shrinkage == "BayesL") {
#         ### LASSO shrinkage of estimates
#         prior = list(varE=list(df=1,S=1.0),
#                      varU=list(df=1,S=1.0), 
#                      varBR=list(df=1,S=1.0), 
#                      lambda=list(shape=0.6,rate=1e-5,value=20,type='random'))
#         mod = BLR(y=y, XL=SNP, prior=prior, nIter=nIter, burnIn=burnIn, thin=thin, saveAt=paste0("temp_", id, "-"))
#         file.remove(paste0("temp_", id, "-lambda.dat"))
#         file.remove(paste0("temp_", id, "-varE.dat"))
#         y_pred = mod$mu[1] + (SNP_test %*% mod$bL)
#     }
#     return(y_pred)
# }


# params = expand.grid(c(1:nfold), c("RRBLUP", "RidgeR", "GLMNET", "LassoR", "BayesU", "BayesR", "BayesL"))
# colnames(params) = c("id", "model")
# N = nrow(params)
# ndigits = length(unlist(strsplit(as.character(N), "")))
# params = params[sample(c(1:N), N, replace=FALSE), ]
# # params = params[1:20,]; N = nrow(params)
# # params = params[params[,2]=="RRBLUP",]; N = nrow(params)


# print("Leave-one-out-cross-validation: ")
# ncores = detectCores()[1]
# clusters = makeCluster(ncores-1)
# registerDoParallel(clusters)
# OUT = foreach(i=1:N, .combine=rbind, .packages=c("rrBLUP", "glmnet", "BLR", "stringr")) %dopar% {
#     # i = 1
#     id = params$id[i]
#     model = params$model[i]
#     y = p[vec_fold_id!=id]
#     SNP = X[vec_fold_id!=id, ]
#     y_test = p[vec_fold_id==id]
#     SNP_test = X[vec_fold_id==id, ]
#     n = length(y_test)
#     if (model == "RRBLUP") {
#         y_pred = RRBLUP(y, SNP, SNP_test)
#     } else if (model == "RidgeR") {
#         y_pred = GLMNET(y, SNP, SNP_test, alpha=0.00)
#     } else if (model == "GLMNET") {
#         y_pred = GLMNET(y, SNP, SNP_test, alpha=0.50)
#     } else if (model == "LassoR") {
#         y_pred = GLMNET(y, SNP, SNP_test, alpha=1.00)
#     } else {
#         y_pred = BAYESREG(y, SNP, SNP_test, shrinkage=model, id=id)
#     }
#     out = data.frame(fold=rep(id, n), model=rep(model, n), y_true=y_test, y_pred=y_pred)
#     write.table(out, file=paste0(model, "-", str_pad(id, ndigits, pad="0"), ".csv.tmp"), sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE)
#     out
# }
# stopCluster(clusters)

# vec_syncx = unlist(strsplit(syncx, "\\."))
# fname_output = paste0(paste(vec_syncx[1:(length(vec_syncx)-1)], collapse="."), "-output.csv")
# write.table(OUT, file=fname_output, sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)

# # ####
# # # bash> cat *.tmp > R_runs_truncated_output_20221110.csv
# # # bash> R
# # dat = read.csv("test_Lr-output.csv", T)
# # dat$diff = dat$y_true - dat$y_pred
# # for (model in unique(dat$model)) {
# #     # model = unique(dat$model)[7]
# #     subdat = dat[dat$model == model, ]
# #     mod = lm(y_true ~ y_pred, data=subdat)
# #     summary_mod = summary(mod)
# #     R2_perc = round(summary_mod$r.squared*100, 2)
# #     ME = round(mean(subdat$diff, na.rm=TRUE), 2)
# #     MAE = round(mean(abs(subdat$diff), na.rm=TRUE), 2)
# #     RMSE = round(sqrt(mean(subdat$diff^2, na.rm=TRUE)), 2)
    
# #     svg(paste0(model, "-GP_LOOCV.svg"), width=10, height=6)
# #     layout(matrix(c(1,1,2,3), ncol=2))
# #     plot(subdat$y_pred, subdat$y_true, xlab="Predicted", ylab="Observed", main=model)
# #         abline(mod)
# #         grid()
# #         legend("topleft", legend=c(paste0("R2=", R2_perc, "%"), paste0("ME=", ME), paste0("MAE=", MAE), paste0("RMSE=", RMSE)), bg="white")
# #     hist(subdat$y_pred, xlab="Predicted", main="Predicted", col=rgb(0.8, 0.1, 0.1, alpha=0.5), bord=NA); grid()
# #     hist(subdat$y_true, xlab="Observed", main="Observed", col=rgb(0.8, 0.1, 0.1, alpha=0.5), bord=NA); grid()
# #     dev.off()
# # }
