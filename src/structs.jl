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

end