# poolgen

Quantitative and population genetics analyses using pool sequencing data

|**Build Status**|**License**|
|:---:|:---:|
| <a href="https://github.com/jeffersonfparil/poolgen/actions"><img src="https://github.com/jeffersonfparil/poolgen/actions/workflows/rust.yml/badge.svg"></a> | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |

## Quickstart

1. Install rustup from [https://www.rust-lang.org/tools/install](https://www.rust-lang.org/tools/install).

2. Download this repository

```shell
git clone https://github.com/jeffersonfparil/poolgen.git
```
3. Compile and run

```shell
cd poolgen/
cargo build --release

cargo run -- sync2syncx \
                -f ./tests/test-pileup2sync_default.sync \
                -o ./tests/test-pileup2sync_default_sync2syncx.syncx \
                --min-cov 2 \
                --n-threads 2
```

## File formats

### Pileup

Summarised or piled up base calls of aligned reads to a reference genome.

- *Column 1*:       name of chromosome, scaffold or contig
- *Column 2*:       locus position
- *Column 3*:       reference allele
- *Column 4*:       coverage, i.e. number of times the locus was sequenced
- *Column 5*:       read codes, i.e. "." ("," for reverse strand) reference allele; "A/T/C/G" ("a/t/c/g" for reverse strand) alternative alleles; "`\[+-][0-9]+[ACGTNacgtn]`" insertions and deletions; "^[" start of read including the mapping quality score; "$" end of read; and "*" deleted or missing locus.
- *Column 6*:       base qualities encoded as the `10 ^ -((ascii value of the character - 33) / 10)`
- *Columns 7 - 3n*: coverages, reads, and base qualities of *n* pools (3 columns per pool).

### Sync

[popoolation2's](https://academic.oup.com/bioinformatics/article/27/24/3435/306737) sync or synchronised pileup file format. A header line showing the names of each column including the names of the pools is prepended by '#'. Additional header line/s and comments prepended with '#' may be added anywhere within the file.

- *Header line/s*:  optional header line/s including the names of the pools, e.g. `# chr pos ref pool1 pool2 pool3 pool4 pool5`
- *Column 1*:       chromosome or scaffold name
- *Column 2*:       locus position 
- *Column 3*:       reference allele, e.g. A, T, C, G 
- *Column 4 to n*:  colon-delimited allele counts: A:T:C:G:DEL:N, where "DEL" is deletion, and "N" is unclassified (one column per pool).

### Phenotypes

1. A simple delimited file, e.g. "csv" and "tsv" with a column for the individual IDs, and at least one column for the phenotypic values. Header line/s and comments should be prepended by '#'.

2. GWAlpha-compatible text file (i.e. "py"):

- *Line 1*: phenotype name
- *Line 2*: standard deviation of the phenotype across pools or for the entire population
- *Line 3*: minimum phenotype value
- *Line 4*: maximum phenotype value
- *Line 5*: cummulative pool sizes percentiles (e.g. 0.2,0.4,0.6,0.8,1.0)
- *Line 6*: phenotype values corresponding to each percentile (e.g. 0.16,0.20,0.23,0.27,0.42)

## Details

### Ordinary least squares regression

The simplest regression model implemented is the ordinary least squares (OLS), where the allele effects are estimated as:

- if $n >= p$, then $\hat{\beta} = (X^{T}X)^{-1} X^{T} y$
- if $n < p$, then $\hat{\beta} = X^{T} (XX^{T})^{-1} y$

where: $n$ is the number of observations, $p$ is the number of predictors, $X$ is an $n \times p$ matrix consisting of a vector of ones and a vector or matrix of allele frequences, and $y$ is the vector of phenotype values.

### GWAlpha: genome-wide estimation of additive effects based on quantile distributions from pool-sequencing experiments

GWAlpha ([Fournier-Level, et al, 2016](https://doi.org/10.1093/bioinformatics/btw805)) iteratively estimates for each locus the effect of each allele on the phenotypic ranking of each pool. This allele effect is defined as $\hat{\alpha} = W {{(\hat{\mu}_{0} - \hat{\mu}_{1})} \over {\sigma_{y}}}$, where:

- let $a$ be the allele in question, and $q$ be the frequency of $a$, then

- $W = 2 \sqrt{q*(1-q)}$ is the penalisation for low allele frequency,

- $\mu_{0}$ is the mean of the beta distribution representing $a$ across pools,

- $\mu_{1}$ is the mean of the beta distribution representing $1-a$ (i.e. additive inverse of $a$) across pools,

- $\sigma_{y}$ is the standard deviation of the phenotype,

- $\beta(\Theta=\{\theta_{1}, \theta_{2}\})$ is used to model the distributions of $a$ and $1-a$ across pools, where:

  - $\Theta$ is estimated via maximum likelihood, i.e. $L(\Theta \mid Q) \sim \Pi^{k}_{i=1} \beta_{pdf}(q_{i} \mid \Theta)$ for the $i^{th}$ pool,

  - $Q = \{q_{1},...,q_{k}\}$ is the cumulative sum of allele frequencies across increasing-phenotypic-value-sorted pools where $k$ is the number of pools, and

  - $\beta_{pdf}(q_{i} \mid \Theta)$ is the probability density function for the $\beta$ distribution.

  - $q_{i} = \beta_{cdf}(y^{\prime}_{i},\Theta) - \beta_{cdf}(y^{\prime}_{i-1},\Theta)$, where $y^{\prime}_{i} \in Y^{\prime}$

  - $Y'$ is the inverse quantile-normalized into phenotype data such that $Y^{\prime} \in [0,1]$.

### Penalised linear regression

I have attempted to create a penalisation algorithm similar to elastic-net or glmnet package ([Friedman et al, 2010](https://doi.org/10.18637/jss.v033.i01)), where ridge regression is implemented if $\alpha = 0$, and lasso if $\alpha = 1$. My implementation is a simplistic contraction of the OLS estimates based on the scaled norms of multiple or univariate OLS, i.e. if the scaled norm is below a $\lambda$ threshold, followed by the expansion of OLS estimates above the threshold. The consequence of this is that while elastic-net uses coordinate descent - an iterative optimisation, my algorithm is parallelisable which means it can potentially outperform glmnet in terms of speed. However, comparisons in terms of accuracy or model fit performance is still to be determined.

### TODO

- [] Fst and Tajima's D estimation

- [] Improve genomic prediction cross-validation's use of PredictionPerformance fields, i.e. output the predictions and predictor distributions

- [] Simulation of genotype and phenotype data

- [] Imputation

Performs OLS or elastic-net regression to predict missing allele counts per window for each pool with at least one locus with missing data. This imputation method requires at least one pool without missing data across the window. It follows that to maximise the number of loci we can impute, we need to impose a maximum window size equal to the length of the sequencing read used to generate the data, e.g. 100 bp to 150 bp for Illumina reads.

For each pool with missing data we model the allele frequencies in the locus with some missing data as:
$$ y_{p} = X_{p}\beta + \epsilon. $$

We then estimate unbiased estimators, $\beta$ as:
$$ \hat{\beta} = f(X_{p}, y_{p{}}), $$

where $f$ can be OLS or elastic-net regression. Imputation is achieved by predicting the missing allele counts:
$$ \hat{y_{m}} = X_{m} \hat{\beta}, $$

where:

- $y_{p}$ is the vector of allele counts of one of the pools with missing data at the loci without missing data (length: mₚ non-missing loci × 7 alleles);
- $X_{p}$ is the matrix of allele counts of pools without missing data at the loci without missing data in the other pools (dimensions: mₚ non-missing loci × 7 alleles, nₚ pools without missing loci);
- $\hat{\beta}$ is the vector of estimates of the effects of each pool without missing data on the allele counts of one of the pools with missing data (length: nₚ pools without missing loci);
- $\hat{y_{m}}$ is the vector of imputed allele counts of one of the pools with missing data (length: mₘ missing loci × 7 alleles); and
- $X_{m}$ is the matrix of allele counts of pools without missing data at the loci with missing data in the other pools (dimensions: mₘ non-missing loci × 7 alleles, nₚ pools without missing loci).

Finally, the imputed allele counts are averaged across the windows sliding one locus at a time.
