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


Y = read.csv(csv, header=TRUE)
y = Y[, ncol(Y)]
# y = log10(y +1) ### Predictions are worse!

n = nrow(X)
nfold = n
vec_id = rep(c(1:nfold), each=floor(n/nfold))
if (length(vec_id) < n) {
    vec_id = c(vec_id, rep(10, times=n-length(vec_id)))
}


RRBLUP = function(y, SNP, SNP_test, method="ML") {
    # i = 5
    # y = y[vec_id!=i]
    # SNP = X[vec_id!=i, ]
    # SNP_test = X[vec_id==i, ]
    # method = "ML"
    mod = mixed.solve(y=y, Z=SNP, method=method)
    return(mod$beta[1] + (SNP_test %*% mod$u))
}

GLMNET = function(y, SNP, SNP_test, alpha=1.00) {
    # i = 5
    # y = y[vec_id!=i]
    # SNP = X[vec_id!=i, ]
    # SNP_test = X[vec_id==i, ]
    # alpha = 1.00
    mod = cv.glmnet(y=y, x=SNP, alpha=alpha)
    idx = which(mod$lambda == mod$lambda.min)    
    return(mod$glmnet.fit$a0[idx] + (SNP_test %*% mod$glmnet.fit$beta[, idx]))
}

BAYESREG = function(y, SNP, SNP_test, shrinkage=c("Bayes-unshrunk", "Bayes-ridge", "Bayes-lasso")[1], nIter=1e4, burnIn=1e2, thin=1e1, id=1) {
    # i = 5
    # y = y[vec_id!=i]
    # SNP = X[vec_id!=i, ]
    # SNP_test = X[vec_id==i, ]
    # nIter=1e3
    # burnIn=1e2
    # thin=1e1
    if (shrinkage == "Bayes-unshrunk") {
        ### Unshruked estimates
        prior = list(varE=list(df=1,S=1.0),
                     varU=list(df=1,S=1.0), 
                     varBR=list(df=1,S=1.0))
        mod = BLR(y=y, XF=SNP, prior=prior, nIter=nIter, burnIn=burnIn, thin=thin, saveAt=paste0("temp_", id, "-"))
        file.remove(paste0("temp_", id, "-bF.dat"))
        file.remove(paste0("temp_", id, "-varE.dat"))
        y_hat = mod$mu[1] + (SNP_test %*% mod$bF)
    } else if (shrinkage == "Bayes-ridge") {
        ### Ridge shrinkage of estimates
        prior = list(varE=list(df=1,S=1.0),
                     varU=list(df=1,S=1.0), 
                     varBR=list(df=1,S=1.0))
        mod = BLR(y=y, XR=SNP, prior=prior, nIter=nIter, burnIn=burnIn, thin=thin, saveAt=paste0("temp_", id, "-"))
        file.remove(paste0("temp_", id, "-varBR.dat"))
        file.remove(paste0("temp_", id, "-varE.dat"))
        y_hat = mod$mu[1] + (SNP_test %*% mod$bR)
    } else if (shrinkage == "Bayes-lasso") {
        ### LASSO shrinkage of estimates
        prior = list(varE=list(df=1,S=1.0),
                     varU=list(df=1,S=1.0), 
                     varBR=list(df=1,S=1.0), 
                     lambda=list(shape=0.6,rate=1e-5,value=20,type='random'))
        mod = BLR(y=y, XL=SNP, prior=prior, nIter=nIter, burnIn=burnIn, thin=thin, saveAt=paste0("temp_", id, "-"))
        file.remove(paste0("temp_", id, "-lambda.dat"))
        file.remove(paste0("temp_", id, "-varE.dat"))
        y_hat = mod$mu[1] + (SNP_test %*% mod$bL)
    }
    return(y_hat)
}


params = paste(y, vec_id, sep=",")
params = expand.grid(params, c("RRBLUP", "Ridge", "GLMNET", "Lasso", "Bayes-unshrunk", "Bayes-ridge", "Bayes-lasso"))
params = data.frame(matrix(as.numeric(unlist(strsplit(as.character(params[,1]), ","))), ncol=2, byrow=TRUE), model=as.character(params[,2]))
colnames(params) = c("y_true", "id", "model")
N = nrow(params)
ndigits = length(unlist(strsplit(as.character(N), "")))
params = params[sample(c(1:N), N, replace=FALSE), ]
# params = params[1:20,]; N = nrow(params)


# Progress combine function (from: https://gist.github.com/kvasilopoulos/d49499ea854541924a8a4cc43a77fed0)
f <- function(iterator){
  pb <- txtProgressBar(min = 1, max = iterator - 1, style = 3)
  count <- 0
  function(...) {
    count <<- count + length(list(...)) - 1
    setTxtProgressBar(pb, count)
    flush.console()
    rbind(...) # this can feed into .combine option of foreach
  }
}

print("Leave-one-out-cross-validation: ")
ncores = detectCores()[1]
clusters = makeCluster(ncores-1)
registerDoParallel(clusters)
OUT = foreach(i=1:N, .combine=f(N), .packages=c("rrBLUP", "glmnet", "BLR", "stringr")) %dopar% {
    # i = 1
    id = params$id[i]
    model = params$model[i]
    y = y[vec_id!=id]
    SNP = X[vec_id!=id, ]
    SNP_test = X[vec_id==id, ]
    if (model == "RRBLUP") {
        y_pred = RRBLUP(y, SNP, SNP_test)
    } else if (model == "Ridge") {
        y_pred = GLMNET(y, SNP, SNP_test, alpha=0.00)
    } else if (model == "GLMNET") {
        y_pred = GLMNET(y, SNP, SNP_test, alpha=0.50)
    } else if (model == "Lasso") {
        y_pred = GLMNET(y, SNP, SNP_test, alpha=1.00)
    } else {
        y_pred = BAYESREG(y, SNP, SNP_test, shrinkage=model, id=id)
    }
    write.table(data.frame(model, id, params$y_true[i], y_pred), file=paste0(model, "-", str_pad(id, ndigits, pad="0"), ".csv.tmp"), sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE)
    y_pred
}
stopCluster(clusters)

OUT = cbind(params, OUT)
colnames(OUT) = c(colnames(OUT)[1:3], "y_pred")
OUT$abs_diff = abs(OUT$y_true - OUT$y_pred)
OUT
