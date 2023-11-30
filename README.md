# poolgen

Quantitative and population genetics analyses using pool sequencing data (i.e. SNP data  where each sample is a pool or group of individuals, a population or a single polyploid individual).

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
    ./target/release/poolgen -h
    ```

4. Detailed documentation

    ```shell
    cargo doc --open
    ```

## File formats

## Note

**Header line/s and comments should be prepended by '#'.**

### Pileup

Summarised or piled up base calls of aligned reads to a reference genome.

- *Column 1*:       name of chromosome, scaffold or contig
- *Column 2*:       locus position
- *Column 3*:       reference allele
- *Column 4*:       coverage, i.e. number of times the locus was sequenced
- *Column 5*:       read codes, i.e. "." ("," for reverse strand) reference allele; "A/T/C/G" ("a/t/c/g" for reverse strand) alternative alleles; "`\[+-][0-9]+[ACGTNacgtn]`" insertions and deletions; "^[" start of read including the mapping quality score; "$" end of read; and "*" deleted or missing locus.
- *Column 6*:       base qualities encoded as the `10 ^ -((ascii value of the character - 33) / 10)`
- *Columns 7 - 3n*: coverages, reads, and base qualities of *n* pools (3 columns per pool).

### Variant call format (vcf)

Canonical variant calling or genotype data format for individual samples. This should include the `AD` field (allele depth), and genotype calls are not required since allele frequencies from allele depth will be used. The input vcf file can be generated with bcftools mpileup like: `bcftools mpileup -a AD...`. The [`vcf2sync`](#vcf2sync) utility is expected to work with vcf versions 4.2 and 4.3. See [VCFv4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf) and [VCFv4.3](https://samtools.github.io/hts-specs/VCFv4.3.pdf) for details in the format specifications.

### Sync

An extension of [popoolation2's](https://academic.oup.com/bioinformatics/article/27/24/3435/306737) sync or synchronised pileup file format, which includes a header line prepended by '#' showing the names of each column including the names of each pool. Additional header line/s and comments prepended with '#' may be added anywhere within the file.

- *Header line/s*:  optional header line/s including the names of the pools, e.g. `# chr pos ref pool1 pool2 pool3 pool4 pool5`
- *Column 1*:       chromosome or scaffold name
- *Column 2*:       locus position 
- *Column 3*:       reference allele, e.g. A, T, C, G 
- *Column/s 4 to n*:  colon-delimited allele counts: A:T:C:G:DEL:N, where "DEL" refers to insertion/deletion, and "N" is unclassified. A pool or population or polyploid individual is represented by a single column of this colon-delimited allele counts.

### Phenotypes

1. A simple delimited file, e.g. "csv" and "tsv" with a column for the individual IDs, and at least one column for the phenotypic values. Header line/s and comments should be prepended by '#'.

2. GWAlpha-compatible text file (i.e. "py"):

- *Line 1*: phenotype name
- *Line 2*: standard deviation of the phenotype across pools or for the entire population
- *Line 3*: minimum phenotype value
- *Line 4*: maximum phenotype value
- *Line 5*: cummulative pool sizes percentiles (e.g. 0.2,0.4,0.6,0.8,1.0)
- *Line 6*: phenotype values corresponding to each percentile (e.g. 0.16,0.20,0.23,0.27,0.42)

## Utilities

### pileup2sync

Convert pileup from `samtools mpileup` into a [synchronised pileup format](#Sync). Pileup from alignments can be generated similar to below:

```shell
samtools mpileup \
    -b /list/of/samtools/-/indexed/bam/or/cram/files.txt \
    -l /list/of/SNPs/in/tab/-/delimited/format/or/bed/-/like.txt \
    -d 100000 \
    -q 30 \
    -Q 30 \
    -f /reference/genome.fna \
    -o /output/file.pileup
```

### vcf2sync

Convert the most widely used genotype data format, [variant call format (`*.vcf`)](https://samtools.github.io/hts-specs/VCFv4.3.pdf) into a [synchronised pileup format](#Sync), making use of allele depths to estimate allele frequencies and omitting genotype classes information including genotype likelihoods. This utility should be compatible with vcf versions 4.2 and 4.3.

### sync2csv

Convert [synchronised pileup format](#Sync) into a matrix ($n$ pools x $p$ alleles across loci) and write into a comma-delimited (csv) file.

<!-- ### impute (redacted for now 2023-11-10)

Impute allele frequencies set to missing according to another minimum depth parameter, i.e. `--min-depth-set-to-missing`. Two imputation algorithms are currently available (a third one is in the works):

1. computationally efficient mean value imputation, and
2. adaptive linkage disequilibrium-based k-nearest neighbour imputation (an extension of [LD-kNNi](https://doi.org/10.1534/g3.115.021667)). -->

### fisher_exact_test

Perform Fisher's exact test per locus.

### chisq_test

Perform Chi-square test per locus.

### pearson_corr

Calculate correlations between allele frequencies per locus and phenotype data.

### ols_iter

Perform ordinary linear least squares regression between allele frequencies and phenotypes per locus, independently.

### ols_iter_with_kinship

Perform ordinary linear least squares regression between allele frequencies and phenotypes using a kinship matrix ($XX' \over p$) as a covariate per locus, independently.

### mle_iter

Perform linear regression between allele frequencies and phenotypes using maximum likelihood estimation per locus, independently.

### mle_iter_with_kinship

Perform linear regression between allele frequencies and phenotypes using maximum likelihood estimation a kinship matrix ($XX' \over p$) as a covariate per locus, independently.

### gwalpha

Perform parametric genomewide association study using pool sequencing data, i.e. pool-GWAS. Refer to [Fournier-Level, et al, 2017](https://academic.oup.com/bioinformatics/article/33/8/1246/2729762) for more details.

### ridge_iter

Perform ridge regression between allele frequencies and phenotypes per locus, independently.

### genomic_prediction_cross_validation

Perform genomic prediction cross-validation using various models including ordinary least squares (OLS), ridge regression (RR), least absolute shrinkage and selection operator (LASSO), and elastic-net ([glmnet](https://glmnet.stanford.edu/articles/glmnet.html)).

### fst

Estimate pairwise genetic differentiation between pools using unbiased estimates of heterozygosity ($\pi$ or $\theta_{\pi}=4N_{e}\mu$ - similar to [Korunes & Samuk 2019](https://doi.org/10.1111/1755-0998.13326) which assumes biallelic loci). Mean genomewide estimates, and per sliding (overlapping/non-overlapping) window estimates are generated.

### heterozygosity

Estimates per sliding window heterozygosities within populations using the unbiased method discussed [above](#fst).

### watterson_estimator

Estimates of Watterson's estimator of $\theta$ which is the expected heterozygosity given the number of polymorphic loci per sliding window (overlappining/non-overlapping).

### tajima_d

Computes [Tajima's D](https://en.wikipedia.org/wiki/Tajima%27s_D) per sliding (overlapping/non-overlapping) window.

### gudmc

Genomewide unbiased determination of the modes of convergent evolution. Per population, significant troughs (selective sweeps) and peaks (balancing selection) are detected and the widths of which are measured. Per population pair, significant deviations from mean genomewide Fst within the identified significant Tajima's D peaks and troughs are also identified. Narrow Tajima's D troughs/peaks imply *de novo* mutation as the source of the genomic segment under selection, while wider ones imply standing genetic variation as the source. At the loci under selection, pairwise Fst which are significantly lower than genomewide Fst imply migration of causal variants between populations, significantly higher implies independent evolution within each population, and non-significantly deviated pairwise Fst implies a shared source of the variants under selection.


## Details

### Ordinary least squares regression

The simplest regression model implemented is the ordinary least squares (OLS), where the allele effects are estimated as:

- if $n >= p$, then $\hat{\beta} = (X^{T}X)^{-1} X^{T} y$
- if $n < p$, then $\hat{\beta} = X^{T} (XX^{T})^{-1} y$

where: $n$ is the number of observations, $p$ is the number of predictors, $X$ is an $n \times p$ matrix consisting of a vector of ones and a vector or matrix of allele frequences, and $y$ is the vector of phenotype values. The inverses are calculated as the [Moore-Penrose pseudoinverse via singular value decomposition](https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Singular_value_decomposition_(SVD)), where the tolerance value uses the machine epsilon (`f64::EPSILON`).

### GWAlpha: genomewide estimation of additive effects based on quantile distributions from pool-sequencing experiments

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

- [X] Fst and Tajima's D estimation
- [X] Imputation
- [X] Canonical variant call format (vcf) file to sync file conversion
- [ ] Simulation of genotype and phenotype data
- [ ] Improve genomic prediction cross-validation's use of PredictionPerformance fields, i.e. output the predictions and predictor distributions
- [ ] Additional imputation algorithm. See details below:

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
