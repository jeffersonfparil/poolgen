######################################
### LEAVE-ONE-OUT-CROSS-VALIDATION ###
######################################

library(rrBLUP)
library(glmnet)
library(BLR)
library(foreach)
library(doParallel)
library(stringr)

syncx = "/data-weedomics-2/poolgen/test/test_Lr.syncx"
csv = "/data-weedomics-2/poolgen/test/test_Lr.csv"
# syncx = "/data-weedomics-2/2.d_40_populations_model_validation/Imputed_geno-mean.syncx"
# csv = "/data-weedomics-2/2.d_40_populations_model_validation/phenotype_data.csv"
setwd(dirname(syncx))

# nfold =  n
# nfold =  ceiling(n/2)
nfold =  10

dat = read.table(syncx, header=FALSE, sep="\t")
### cbind is faster than rbind
C = matrix(as.numeric(unlist(strsplit(dat[,3], ":"))), ncol=7, byrow=TRUE)
X = matrix(t(C / rowSums(C)), ncol=1, byrow=FALSE)
print("Loading syncx allele counts matrix: ")
pb = txtProgressBar(4, ncol(dat), style=3)
for (i in 4:ncol(dat)) {
    setTxtProgressBar(pb, i)
    C = matrix(as.numeric(unlist(strsplit(dat[,i], ":"))), ncol=7, byrow=TRUE)
    F = matrix(t(C / rowSums(C)), ncol=1, byrow=FALSE)
    X = cbind(X, F)
}
close(pb)
X = t(X)
AF_mu = apply(X, MARGIN=2, FUN=mean)
AF_sd = apply(X, MARGIN=2, FUN=sd)
idx = (AF_mu > 0) & (AF_sd > 0)
X = X[, idx]



Y = read.csv(csv, header=TRUE)
p = Y[, ncol(Y)]
# y = log10(y +1) ### Predictions are worse!

n = nrow(X)
vec_fold_id = rep(c(1:nfold), each=floor(n/nfold))
if (length(vec_fold_id) < n) {
    vec_fold_id = c(vec_fold_id, rep(nfold, times=n-length(vec_fold_id)))
}


RRBLUP = function(y, SNP, SNP_test, method="ML") {
    # i = 1
    # y = p[vec_fold_id!=i]
    # SNP = X[vec_fold_id!=i, ]
    # SNP_test = X[vec_fold_id==i, ]
    # method = "ML"
    mod = mixed.solve(y=y, Z=SNP, method=method)
    y_pred = mod$beta[1] + (SNP_test %*% mod$u)
    return(y_pred)
}

GLMNET = function(y, SNP, SNP_test, alpha=1.00) {
    # i = 5
    # y = p[vec_fold_id!=i]
    # SNP = X[vec_fold_id!=i, ]
    # SNP_test = X[vec_fold_id==i, ]
    # alpha = 1.00
    mod = cv.glmnet(y=y, x=SNP, alpha=alpha)
    idx = tryCatch(
        idx = which(mod$lambda == mod$lambda.min),
        error = function(e){
            idx = which(mod$cvm == min(mod$cvm, na.rm=TRUE))
        }
    )
    if (length(idx) == 0){
        y_pred = NA
    } else {
        y_pred = mod$glmnet.fit$a0[idx] + (SNP_test %*% mod$glmnet.fit$beta[, idx])
    }
    return(y_pred)
}

BAYESREG = function(y, SNP, SNP_test, shrinkage=c("BayesU", "BayesR", "BayesL")[1], nIter=1e4, burnIn=1e2, thin=1e1, id=1) {
    # i = 5
    # y = p[vec_fold_id!=i]
    # SNP = X[vec_fold_id!=i, ]
    # SNP_test = X[vec_fold_id==i, ]
    # nIter=1e3
    # burnIn=1e2
    # thin=1e1
    if (shrinkage == "BayesU") {
        ### Unshruked estimates
        prior = list(varE=list(df=1,S=1.0),
                     varU=list(df=1,S=1.0), 
                     varBR=list(df=1,S=1.0))
        mod = BLR(y=y, XF=SNP, prior=prior, nIter=nIter, burnIn=burnIn, thin=thin, saveAt=paste0("temp_", id, "-"))
        file.remove(paste0("temp_", id, "-bF.dat"))
        file.remove(paste0("temp_", id, "-varE.dat"))
        y_pred = mod$mu[1] + (SNP_test %*% mod$bF)
    } else if (shrinkage == "BayesR") {
        ### Ridge shrinkage of estimates
        prior = list(varE=list(df=1,S=1.0),
                     varU=list(df=1,S=1.0), 
                     varBR=list(df=1,S=1.0))
        mod = BLR(y=y, XR=SNP, prior=prior, nIter=nIter, burnIn=burnIn, thin=thin, saveAt=paste0("temp_", id, "-"))
        file.remove(paste0("temp_", id, "-varBR.dat"))
        file.remove(paste0("temp_", id, "-varE.dat"))
        y_pred = mod$mu[1] + (SNP_test %*% mod$bR)
    } else if (shrinkage == "BayesL") {
        ### LASSO shrinkage of estimates
        prior = list(varE=list(df=1,S=1.0),
                     varU=list(df=1,S=1.0), 
                     varBR=list(df=1,S=1.0), 
                     lambda=list(shape=0.6,rate=1e-5,value=20,type='random'))
        mod = BLR(y=y, XL=SNP, prior=prior, nIter=nIter, burnIn=burnIn, thin=thin, saveAt=paste0("temp_", id, "-"))
        file.remove(paste0("temp_", id, "-lambda.dat"))
        file.remove(paste0("temp_", id, "-varE.dat"))
        y_pred = mod$mu[1] + (SNP_test %*% mod$bL)
    }
    return(y_pred)
}


params = expand.grid(c(1:nfold), c("RRBLUP", "RidgeR", "GLMNET", "LassoR", "BayesU", "BayesR", "BayesL"))
colnames(params) = c("id", "model")
N = nrow(params)
ndigits = length(unlist(strsplit(as.character(N), "")))
params = params[sample(c(1:N), N, replace=FALSE), ]
# params = params[1:20,]; N = nrow(params)
# params = params[params[,2]=="RRBLUP",]; N = nrow(params)


print("Leave-one-out-cross-validation: ")
ncores = detectCores()[1]
clusters = makeCluster(ncores-1)
registerDoParallel(clusters)
OUT = foreach(i=1:N, .combine=rbind, .packages=c("rrBLUP", "glmnet", "BLR", "stringr")) %dopar% {
    # i = 1
    id = params$id[i]
    model = params$model[i]
    y = p[vec_fold_id!=id]
    SNP = X[vec_fold_id!=id, ]
    y_test = p[vec_fold_id==id]
    SNP_test = X[vec_fold_id==id, ]
    n = length(y_test)
    if (model == "RRBLUP") {
        y_pred = RRBLUP(y, SNP, SNP_test)
    } else if (model == "RidgeR") {
        y_pred = GLMNET(y, SNP, SNP_test, alpha=0.00)
    } else if (model == "GLMNET") {
        y_pred = GLMNET(y, SNP, SNP_test, alpha=0.50)
    } else if (model == "LassoR") {
        y_pred = GLMNET(y, SNP, SNP_test, alpha=1.00)
    } else {
        y_pred = BAYESREG(y, SNP, SNP_test, shrinkage=model, id=id)
    }
    out = data.frame(fold=rep(id, n), model=rep(model, n), y_true=y_test, y_pred=y_pred)
    write.table(out, file=paste0(model, "-", str_pad(id, ndigits, pad="0"), ".csv.tmp"), sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE)
    out
}
stopCluster(clusters)

vec_syncx = unlist(strsplit(syncx, "\\."))
fname_output = paste0(paste(vec_syncx[1:(length(vec_syncx)-1)], collapse="."), "-output.csv")
write.table(OUT, file=fname_output, sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)

# ####
# # bash> cat *.tmp > R_runs_truncated_output_20221110.csv
# # bash> R
# dat = read.csv("test_Lr-output.csv", T)
# dat$diff = dat$y_true - dat$y_pred
# for (model in unique(dat$model)) {
#     # model = unique(dat$model)[7]
#     subdat = dat[dat$model == model, ]
#     mod = lm(y_true ~ y_pred, data=subdat)
#     summary_mod = summary(mod)
#     R2_perc = round(summary_mod$r.squared*100, 2)
#     ME = round(mean(subdat$diff, na.rm=TRUE), 2)
#     MAE = round(mean(abs(subdat$diff), na.rm=TRUE), 2)
#     RMSE = round(sqrt(mean(subdat$diff^2, na.rm=TRUE)), 2)
    
#     svg(paste0(model, "-GP_LOOCV.svg"), width=10, height=6)
#     layout(matrix(c(1,1,2,3), ncol=2))
#     plot(subdat$y_pred, subdat$y_true, xlab="Predicted", ylab="Observed", main=model)
#         abline(mod)
#         grid()
#         legend("topleft", legend=c(paste0("R2=", R2_perc, "%"), paste0("ME=", ME), paste0("MAE=", MAE), paste0("RMSE=", RMSE)), bg="white")
#     hist(subdat$y_pred, xlab="Predicted", main="Predicted", col=rgb(0.8, 0.1, 0.1, alpha=0.5), bord=NA); grid()
#     hist(subdat$y_true, xlab="Observed", main="Observed", col=rgb(0.8, 0.1, 0.1, alpha=0.5), bord=NA); grid()
#     dev.off()
# }
