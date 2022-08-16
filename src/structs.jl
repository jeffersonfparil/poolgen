### Naming convention:
### (1) variable names: snake_case
### (2) structure names: CamelCase
### (3) function names: SCREAMING
### (4) function names user-exposed: snake_case

module structs

struct PileupLine
    index::Int      ### line number
    line::String    ### a line of the pileup file
end

struct SyncxLine
    index::Int      ### line number (not critical)
    line::String    ### a line of the syncx file
end

struct LocusAlleleCounts
    chr::String         ### chromosome
    pos::Int            ### position
    ref::Char           ### reference allele
    dep::Vector{Int}    ### vector of depths per pool
    A::Vector{Int}      ### counts per pool of the A allele
    T::Vector{Int}      ### counts per pool of the T allele
    C::Vector{Int}      ### counts per pool of the C allele
    G::Vector{Int}      ### counts per pool of the G allele
    I::Vector{Int}      ### counts per pool of the I allele (insertion)
    D::Vector{Int}      ### counts per pool of the D allele (deletion)
    N::Vector{Int}      ### counts per pool of the N allele (missing)
end

mutable struct Window
    ### Note that this Window struct does not really even need the mutable keyword since its matrix and vector componenets are mutable it seems
    chr::Vector{String} ### vector of chromosome names
    pos::Vector{Int}    ### vector of positions
    ref::Vector{Char}   ### vector of reference alleles
    cou::Matrix{Any}    ### n (window size*number of alleles) rows x p (number of pools) columns
    imp::Matrix{Int}    ### number of times a locus has been imputed (corresponds to the elements of cou)
end

struct PhenotypeLine
    index::Int          ### line number (not critical)
    line::String        ### a line of delimited phenotype text file
    delimiter::String   ### string (single or multiple characters) delimeter
    id_col::Int         ### column number of ID names corresponsing to the columns in the genotype file
    y_cols::Vector{Int} ### vector of column numbers corresponding to the phenotype columns needed
end

struct Phenotype
    ids::Vector{String}         ### names of individuals (n individuals)
    trait_names::Vector{String} ### names of traits (m traits)
    Y::Matrix{Any}              ### n x m matrix of phenotype values
end

end