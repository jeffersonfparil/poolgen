# poolgen
Quantitative and population genetics analyses using pool sequencing data

|**Build Status**|**License**|
|:---:|:---:|
| <a href="https://github.com/jeffersonfparil/poolgen/actions"><img src="https://github.com/jeffersonfparil/poolgen/actions/workflows/julia.yml/badge.svg"></a> | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |

## Quickstart

1. Download and install Julia. See [https://julialang.org/downloads/](https://julialang.org/downloads/).

2. Install the poolgen package. Open Julia and run the following:
```julia
using Pkg
Pkg.add(url="https://github.com/jeffersonfparil/poolgen.git")
```

3. Usage and documentation access. Open Julia and browse the documentation for all the currently available functions:
```julia
using poolgen
?poolgen
?poolgen.convert
?poolgen.pileup2syncx
?poolgen.filter
?poolgen.filter
?poolgen.impute
?poolgen.simulate
?poolgen.gwalpha
?poolgen.genomic_prediction
?poolgen.genomic_prediction_CV
```
## File formats

### Pileup
Summarised or piled up base calls of aligned reads to a reference genome.
- *Column 1*:       name of chromosome, scaffold or contig
- *Column 2*:       locus position
- *Column 3*:       reference allele
- *Column 4*:       coverage, i.e. number of times the locus was included in a read
- *Column 5*:       Read codes, i.e. "." ("," for reverse strand) reference allele; "A/T/C/G" ("a/t/c/g" for reverse strand) alternative alleles; "`\[+-][0-9]+[ACGTNacgtn]`" insertions and deletions; "^" start of read; "$" end of read; and "*" deleted or missing locus.
- *Column 6*:       base qualities encoded as the `10 ^ -((ascii value of the character - 33) / 10)`
- *Columns 7 - 3n*: coverages, reads, and base qualities of *n* pools (3 columns per pool).

### Syncx
Spiritual successor to [popoolation2's](https://academic.oup.com/bioinformatics/article/27/24/3435/306737) sync or synchronised pileup file format:
- *Column 1*:      chromosome or scaffold name
- *Column 2*:      locus position 
- *Column 3 to n*: colon-delimited allele counts: A:T:C:G:INS:DEL:N, where "INS" is insertion, "DEL" is deletion, and "N" is unclassified (one column per pool).

### Phenotypes

1. A simple delimited file, e.g. "csv" and "tsv" with a column for the individual IDs, and at least one column for the phenotypic values
2. Text with a ".py" extension for GWAlpha:
```python
	Pheno_name='Phenotype Name';
	sig=0.06724693662723039; # standard deviation
	MIN=0.0; # minimum phenotype value
	MAX=0.424591738712776; # maximum phenotype value
	perc=[0.2,0.4,0.6,0.8]; # cummulative pool sizes percentiles excluding the last pool
	q=[0.16,0.20,0.23,0.27,0.42]; # phenotype values corresponding to each percentile
```

## Details

### Ordinary least squares regression

The simplest regression model implemented is the ordinary least squares (OLS), where the allele effects are estimated as $\hat{\beta} = (X^{T}X)^{-1} X^{T} y$, where: $X$ consists of a vector of ones and a vector or matrix of allele frequences, and $y$ is the vector of phenotype values.

### GWAlpha: genome-wide estimation of additive effects based on quantile distributions from pool-sequencing experiments

GWAlpha ([Fournier-Level, et al, 2016](https://doi.org/10.1093/bioinformatics/btw805)) iteratively estimates for each locus the effect of each allele on the phenotypic ranking of each pool. This allele effect is defined as $\hat{\alpha} = W \frac{(\hat{\mu}_{0} - \hat{\mu}_{1})} {\sigma_{y}}$, where:
- let $a$ be the allele in question, and $q$ be the frequency of $a$, then
- $W = 2 \sqrt{q*(1-q)}$ is the penalisation for low allele frequency,
- $\mu_{0}$ is the mean of the beta distribution representing $a$ across pools,
- $\mu_{1}$ is the mean of the beta distribution representing $1-a$ (i.e. additive inverse of $a$) across pools,
- $\sigma_{y}$ is the standard deviation of the phenotype,
- $\beta(\Theta=\{\theta_{1}, \theta_{2}\})$ is used to model the distributions of $a$ and $1-a$ across pools, where:
  + $\Theta$ is estimated via maximum likelihood, i.e. $L(\Theta \mid Q) \propto \Pi^{k}_{i=1} \beta_{pdf}(q_{i} \mid \Theta)$ for the $i^{th}$ pool,
  + $Q = \{q_{1},...,q_{k}\}$ is the cumulative sum of allele frequencies across increasing-phenotypic-value-sorted pools where $k$ is the number of pools, and
  + $\beta_{pdf}(q_{i} \mid \Theta)$ is the probability density function for the $\beta$ distribution.
  + $q_{i} = \beta_{cdf}(y'_{i},\Theta) - \beta_{cdf}(y'_{i-1},\Theta)$, where $y'_{i} ∈ Y'$
 + $Y'$ is the inverse quantile-normalized into phenotype data such that $Y' \in [0,1]$.

### Linear mixed models

The linear mixed model is defined by [Henderson in 1975](https://doi.org/10.2307/2529430) is $y = X\beta + Zu + \epsilon$, where:
- $y$ isthe vector of phenotype values;
- $X\beta$ are the fixed effects, where in:
  + **GBLUP** (genomic best linear unbiased prediction): $X$ consists of a column of ones and allele frequencies matrix or vector, and $\beta$ is the vector of **intercept and allele effects**, and in
  + **ABLUP** (additive genetic effects best linear unbiased prediction) and **RRBLUP** (ridge regression best linear unbiased prediction): $X$ consists of a column of ones and non-genomic covariates, if available e.g. kinship/pool relationship matrix, testing environments, time, and blocks, and $\beta$ is the vector of **intercept and covariate effects**, if available;
- $Zu$ are the random effects, where in:
  + **GBLUP**: $Z$ is an identity matrix, and $u$ are the effects of each pool, where:
    - $u \sim N(0, G=\sigma^{2}_{g}K)$, and
    - $\sigma^{2}_{g}K$ relationship matrix between pools; and in
  + **ABLUP** and **RRBLUP**: $Z$ is the allele frequencies matrix, and $u$ is the vector of allele effects, where:
    - $u \sim N(0, G=\sigma^{2}_{a}I)$, and 
    - $\sigma^{2}_{a}I$ is the additive variance-covarince matrix showing identicaly and independently distributed allele effects;
- $\epsilon$ is the identically and independently distributed residual effects vector, i.e. $\epsilon \sim N(0, R=\sigma^{2}I)$
- variance components ($\sigma^{2}$, and $\sigma^{2}_{g}$ in GBLUP, or $\sigma^{2}_{a}$ in ABLUP and RRBLUP) are estimated via maximum likelihood (ML) or restricted maximum likelihood (REML);
- let $n$ be the number of pools, $m$ be the number if fixed effects to estimate, $p$ be the number of random effects to estmate, and $V = (Z G Z^{T}) + R$;
- fixed effects are estimated by solving:
  + $\hat{\beta} = (X^{T} V^{-1} X)^{-1} (X^{T} V^{-1} y)$, if $n > m$, while
  + $\hat{\beta} = (X^{T} V^{-1}) (X X^{T} V^{-1})^{-1} y$, if $n << m$ which is the default as most pool sequencing experiments have this dimensionality problem; and finally
- random effects are estimated by solving:
  + **GBLUP** and **ABLUP**: $\hat{\mu} = (G Z^{T} V^{-1}) (y - X\hat{\beta})$, or
  + **RRBLUP**: $\hat{\mu} = (Z^{T}Z + \lambda I)^{-1} (Z^{T}y)$, where $\lambda = \sigma^{2}/\sigma^{2}_{a}$.

### Elastic-net regression using glmnet

Elastic net performs penalised maximum likelihood and is implemented on the glmnet package ([Friedman et al, 2010](https://doi.org/10.18637/jss.v033.i01)), where ridge regression is implemented if $\alpha = 0$, and lasso if $\alpha = 1$. For details see: Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2010. “Regularization Paths for Generalized Linear Models via Coordinate Descent.” Journal of Statistical Software, Articles 33 (1): 1–22. [https://doi.org/10.18637/jss.v033.i01](https://doi.org/10.18637/jss.v033.i01).

### Imputation details

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

## Acknowledgements

- [Benjamin Camm](https://adaptive-evolution.biosciences.unimelb.edu.au/people/bencamm.html)

![ben](misc/BJC_pic.jpg)

- [Adaptive Evolution lab](https://adaptive-evolution.biosciences.unimelb.edu.au)
- [Grains Research & Development Corporation](https://grdc.com.au/)

## To do list
- [X] Types: pileup line
- [X] Types: syncx line
- [X] Types: phenotype line
- [X] Types: locus
- [X] Types: window
- [X] Types: phenotype
- [X] pileup I/O
- [X] syncx I/O
- [X] pileup filtering: pileup to pileup and pileup to syncx
- [X] syncx filtering
- [X] imputation pileup to syncx
- [X] simple genotype and phenotype simulation framework (with LD)
- [X] OLS
- [X] Elastic-net
- [X] LMM: gBLUP
- [X] LMM: rrBLUP
- [X] GP cross-validation and plotting
- [X] GWAlpha
- [X] iterative OLS
- [X] iterative elastic-net
- [X] iterative LMM: gBLUP
- [X] iterative LMM: rrBLUP
- [X] Empirical p-values via bootstrapping
- [X] Empirical p-values via maximum likelihood
- [X] GWAS plotting
- [ ] Machine-learning: random forest
- [ ] Machine-learning: multilayer perceptron