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

## Checklist
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
- [ ] Empirical p-values via bootstrapping
- [ ] Empirical p-values via maximum likelihood
- [ ] GWAS plotting
- [ ] Machine-learning: random forest
- [ ] Machine-learning: multilayer perceptron

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
	sig=0.06724693662723039;		  # standard deviation
	MIN=0.0;				              # minimum phenotype value
	MAX=0.424591738712776;			  # maximum phenotype value
	perc=[0.2,0.4,0.6,0.8];			  # cummulative pool sizes percentiles excluding the last pool
	q=[0.16,0.20,0.23,0.27,0.42];	# phenotype values corresponding to each percentile
```

## Details

### GWAlpha details (needs updating)
The GWAlpha model is defined as **α = W(μₐₗₗₑₗₑ-μₐₗₜₑᵣₙₐₜᵢᵥₑ)/σᵧ**, where:
- **μ** is the mean of the beta distribution, **Beta(θ)** where **θ={θ₁,θ₂}**
- **θ** is estimated via maximum likelihood **L(θ|Q) ∝ πᵢ₌₁₋ₖf(qᵢ|θ)**
- **Q = {q₁,...,qₖ}** is the cumulative sum of allele frequencies across increasing-phenotypic-value-sorted pools where **k** is the number of pools
- **E(allele|θ) = Beta_cdf(yᵢ',θ) - Beta_cdf(yᵢ₋₁',θ)**, where **yᵢ' ∈ Y'**
- **Y'** is the inverse quantile-normalized into phenotype data such that **Y' ∈ [0,1]**
- **W = 2√{E(allele)*(1-E(allele))}** is the penalization for low allele frequency
Cite: Fournier-Level A, Robin C, Balding DJ (2016). [GWAlpha: Genome-Wide estimation of additive effects (Alpha) based on trait quantile distribution from pool-sequencing experiments.](https://doi.org/10.1093/bioinformatics/btw805)

### Linear mixed models (needs updating)
The linear mixed model is defined as **y = Xb + Zu + e**, where:
- **X** [n,p] is the centered matrix of allele frequencies
- **Z** [n,n] is the square symmetric matrix of relatedness
- **y** [n,1] is the centered vector of phenotypic values
- no intercept is explicitly fitted but implicitly set at the mean phenotypic value as a consequence of centering **y**
- **u ~ N(0, σ²uI)**
- **e ~ N(0, σ²eI)**
- **y ~ N(0, V)**
- **V = (Z (σ²uI) Z') + (σ²eI)**
- variance components (**σ²e**, **σ²u**) are estimated via maximum likelihood (ML) or restricted maximum likelihood (REML)
- fixed effects (**b**) are estimated by solving: **(X' * V⁻¹) * (X * X' * V⁻¹)⁻¹ * y**
- random effects (**u**) are estimated by solving: **(σ²uI) * Z' * V⁻¹ * (y - (X*b))**

Empirical p-values were calculated by modelling the additive allelic effects (α) using a normal distribution with mean and variance estimated using maximum likelihood.

### Elastic-net details (needs updating)
- ridge regression at α = 0.0
- LASSO regression at α = 0.0
- See: Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2010. “Regularization Paths for Generalized Linear Models via Coordinate Descent.” Journal of Statistical Software, Articles 33 (1): 1–22. https://doi.org/10.18637/jss.v033.i01.

### Imputation details
Performs simple linear regression to predict missing allele counts per window for each pool with at least one locus with missing data. This imputation method requires at least one pool without missing data across the window. It follows that to maximise the number of loci we can impute, we need to impose a maximum window size equal to the length of the sequencing read used to generate the data, e.g. 100 bp to 150 bp for Illumina reads.

- For each pool with missing data we estimate β̂ as:
```
          yₚ = Xₚβ
        → β̂ = inverse(XₚᵀXₚ) (Xₚᵀyₚ).
```

- For each pool with missing data, imputation is achieved by predicting the missing allele counts:
```
          ŷₘ = XₘB̂.
```
- Where:
    + **yₚ** is the vector of allele counts of one of the pools with missing data at the loci without missing data (length: mₚ non-missing loci × 7 alleles);
    + **Xₚ** is the matrix of allele counts of pools without missing data at the loci without missing data in the other pools (dimensions: mₚ non-missing loci × 7 alleles, nₚ pools without missing loci);
    + **β̂** is the vector of estimates of the effects of each pool without missing data on the allele counts of one of the pools with missing data (length: nₚ pools without missing loci);
    + **inverse()** is the Moore-Penrose pseudoinverse if the automatic Julia solver fails;
    + **ŷₘ** is the vector of imputed allele counts of one of the pools with missing data (length: mₘ missing loci × 7 alleles); and
    + **Xₘ** is the matrix of allele counts of pools without missing data at the loci with missing data in the other pools (dimensions: mₘ non-missing loci × 7 alleles, nₚ pools without missing loci).

- The imputed allele counts are averaged across the windows sliding one locus at a time.
