# poolgen
Quantitative and population genetics analyses using pool sequencing data

|**Build Status**|**License**|
|:---:|:---:|
| <a href="https://github.com/jeffersonfparil/poolgen/actions"><img src="https://github.com/jeffersonfparil/poolgen/actions/workflows/julia.yml/badge.svg"></a> | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |

## High-level plan
- Phase 1: code in Julia
- Phase 2: try to implement in Rust

## Plan
- [X] Types: (1/3) pileup line
- [X] Types: (2/3) locus
- [X] Types: (3/3) window
- [X] pileup I/O
- [X] syncx I/O
- [X] pileup filtering: pileup to pileup and pileup to syncx
- [X] syncx filtering
- [X] imputation pileup-to-syncx
- [ ] iterative genome-wide association analysis (per locus)
- [ ] genomic prediction modeling and cross-validation

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

## Imputation details
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
