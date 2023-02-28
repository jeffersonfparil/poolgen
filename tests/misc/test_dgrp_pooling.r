args = commandArgs(trailingOnly=TRUE)
# args = c("DGRP_phenotype_lifetime_fecundity.csv", "1", "2", "Lifetime_fecundity", "5", "0.2", "0.2", "0.2", "0.2", "0.2")
fname = args[1]
line_col = as.numeric(args[2])
phen_col = as.numeric(args[3])
phen_name = args[4]
n_pools = as.numeric(args[5])
pool_sizes = as.numeric(args[6:(n_pools+5)])

root_name = paste(head(unlist(strsplit(basename(fname), '[.]')), -1), collapse='.')

dat = read.csv(fname, header=FALSE)
dat = data.frame(Line=dat[,line_col], Pheno=dat[,phen_col])
colnames(dat) = c("Line", phen_name)
dat = dat[order(dat[,2], decreasing=FALSE), ]

### REMOVE MISSING DGRP LINES MISSING IN THE GENOTYPE DATA
dat = dat[dat$Line != "RAL_591", ]
dat = dat[dat$Line != "RAL_272", ]
dat = dat[dat$Line != "RAL_398", ]
dat = dat[dat$Line != "RAL_556", ]
### Something went wrong with the indexing the bams
dat = dat[dat$Line != "RAL_379", ]
dat = dat[dat$Line != "RAL_852", ]
dat = dat[dat$Line != "RAL_476", ]
dat = dat[dat$Line != "RAL_313", ]
dat = dat[dat$Line != "RAL_589", ]
dat = dat[dat$Line != "RAL_378", ]

n = nrow(dat)
pool_sizes = cumsum(floor(pool_sizes * n))
if (pool_sizes[n_pools] < n) {
    idx = ceiling(n_pools/2)
    add = n - pool_sizes[n_pools]
    pool_sizes[idx:n_pools] = pool_sizes[idx:n_pools] + add
}

pool_sizes = c(0, pool_sizes)
vec_perc = head(tail(pool_sizes, -1), -1) / n
MIN = min(dat[,2])
MAX = max(dat[,2])
STD = sd(dat[,2])
vec_lines = list()
vec_means = c()
vec_sd = c()
vec_q = c()
for (i in 1:n_pools) {
    # i = 1
    start = pool_sizes[i] + 1
    end = pool_sizes[i+1]
    eval(parse(text=paste0("vec_lines$Pool_", i, " = dat[start:end, 1]")))
    y = dat[start:end, 2]
    vec_means = c(vec_means, mean(y))
    vec_sd = c(vec_sd, sd(y))
    vec_q = c(vec_q, max(y))
}
vec_q = head(vec_q, -1)
py = c(paste0('Pheno_name= "', phen_name, '";'),
       paste0('sig= ', STD, ';'),
       paste0('MIN= ', MIN, ';'),
       paste0('MAX= ', MAX, ';'),
       paste0('perc= [', paste(vec_perc, collapse=','), '];'),
       paste0('q= [', paste(vec_q, collapse=','), '];')
      )
csv = data.frame(Pool=c(1:n_pools), y=vec_means)
colnames(csv) = c("Pool", phen_name)

### List of bams to merge per pool
for (i in 1:length(vec_lines)) {
    # i = 1
    bam_list = paste0(vec_lines[[i]], '*.bam')
    write.table(bam_list, file=paste0(root_name, "_Pool_", i, "_bam_list.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
}

### Phenotype files, i.e. py for GWAlpha.py and csv for poolgen
write.table(py, file=paste0(root_name, "_pheno.py"), quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(csv, file=paste0(root_name, "_pheno.csv"), sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
