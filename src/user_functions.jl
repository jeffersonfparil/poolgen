### Naming convention:
### (1) variable names: snake_case
### (2) structure names: CamelCase
### (3) function names: SCREAMING
### (4) function names user-exposed: snake_case

module user_functions

using Distributed
using ProgressMeter

include("functions.jl")
using .functions: PileupLine, SyncxLine, LocusAlleleCounts, Window
using .functions: PARSE, SAVE, SPLIT, MERGE, CONVERT, PILEUP2SYNCX, FILTER, IMPUTE, GWALPHA
using .functions: GWAS_PLOT_MANHATTAN, ESTIMATE_LOD
using .functions: SIMULATE, POOL, EXPORT_SIMULATED_DATA
using .functions: OLS_MULTIVAR, ELA_MULTIVAR, LMM_MULTIVAR
using .functions: CV_MULTIVAR

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
###########
### I/O ###
###########
"""
# ___________________________________
# Convert synchronised pileup formats from sync to syncx and vice-versa

`convert(syncx_or_sync::String; out::String="")::String`

# Inputs
1. syncx_or_sync [String]: synchronised pileup file, i.e. sync or syncx
2. out [String]: output filename (defaults to the basename of `syncx_or_sync` with the extension name converted from sync to syncx and vice-versa)

# Output
1. [String]: filename of the output

# Examples
1. Single-core conversion
```julia
using poolgen
n=5; m=10_000; l=135_000_000; k=5; ϵ=Int(1e+15); a=2; vec_chr_lengths=[0]; vec_chr_names=[""]; dist_noLD=500_000; o=100; t=10; nQTL=10; heritability=0.5; LD_chr=""; LD_n_pairs=10_000; plot_LD=false; npools=5
map, bim, ped, fam, syncx, csv = poolgen.simulate(n=n, m=m, l=l, k=k, ϵ=ϵ, a=a, vec_chr_lengths=vec_chr_lengths, vec_chr_names=vec_chr_names, dist_noLD=dist_noLD, o=o, t=t, nQTL=nQTL, heritability=heritability, npools=npools, LD_chr=LD_chr, LD_n_pairs=LD_n_pairs, plot_LD=plot_LD)
sync_out = poolgen.convert(syncx, out="test.sync")
sync_out = poolgen.convert(syncx, out="test.syncx")
```

2. Multi-core conversion
```julia
using Distributed
Distributed.addprocs(length(Sys.cpu_info())-1)
@everywhere using poolgen
n=5; m=10_000; l=135_000_000; k=5; ϵ=Int(1e+15); a=2; vec_chr_lengths=[0]; vec_chr_names=[""]; dist_noLD=500_000; o=100; t=10; nQTL=10; heritability=0.5; LD_chr=""; LD_n_pairs=10_000; plot_LD=false; npools=5
map, bim, ped, fam, syncx, csv = poolgen.simulate(n=n, m=m, l=l, k=k, ϵ=ϵ, a=a, vec_chr_lengths=vec_chr_lengths, vec_chr_names=vec_chr_names, dist_noLD=dist_noLD, o=o, t=t, nQTL=nQTL, heritability=heritability, npools=npools, LD_chr=LD_chr, LD_n_pairs=LD_n_pairs, plot_LD=plot_LD)
sync_out = poolgen.convert(syncx, out="test.sync")
sync_out = poolgen.convert(syncx, out="test.syncx")
```
"""
function convert(syncx_or_sync::String; out::String="")::String
    # using Distributed
    # Distributed.addprocs(length(Sys.cpu_info())-1)
    # @everywhere using ProgressMeter
    # @everywhere include("functions.jl")
    # @everywhere using .functions: PileupLine, LocusAlleleCounts, Window
    # @everywhere using .functions: SPLIT, MERGE, CONVERT
    # syncx_or_sync = "/home/jeffersonfparil/Documents/poolgen/test/test_3.syncx"
    # # syncx_or_sync = "/home/jeffersonfparil/Documents/poolgen/test/test_3-A.syncx"
    # out = ""
    ### Define output file if not specified
    if out == ""
        extension = [split(syncx_or_sync, ".")[end][end] == 'x' ? ".sync" : ".syncx"][1]
        out = string(join(split(syncx_or_sync, ".")[1:(end-1)], "."), extension)
    end
    ### Find file positions for parallel processing
    threads = length(Distributed.workers())
    positions_init, positions_term = SPLIT(threads, syncx_or_sync)
    digit = length(string(length(positions_init)))
    ### Convert into syncx or sync format
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        id = lpad(i, (digit+1)-length(string(i)), "0")
        tmp = string(syncx_or_sync, "-CONVERT2SYNCXOR2SYNC-", id, ".tmp")
        filename = CONVERT(syncx_or_sync, init, term, tmp)
        [filename]
    end
    ### Sort the output files from parallel processing and merge
    MERGE(filenames_out, out)
    return(out)
end

"""
# ___________________________________
# Convert pileup into the extended synchronised pileup format

`pileup2syncx(pileup::String; out::String="")::String`

# Inputs
1. pileup [String]: pileup file
2. out [String]: output filename (defaults to the basename of `pileup` with the extension name converted to syncx)

# Output
1. [String]: filename of the output

# Examples
1. Single-core parsing
```julia
using poolgen
pileup = "test.pileup"
f = open(pileup, "a")
write(f, "seq1\t272\tT\t24\t,......,,.,.,...,,,.,..^+.\t<<<+;<<<<<<<<<<<=<;<;7<&\t24\t,......,,.,.,...,,,.,..^+.\t<<<+;<<<<<<<<<<<=<;<;7<&\n")
write(f, "seq1\t273\tT\t23\t,.....,,.,.,...,,,.,..A\t<<<;<<<<<<<<<3<=<<<;<<+\t23\t,.....,,.,.,...,,,.,..A\t<<<;<<<<<<<<<3<=<<<;<<+\n")
write(f, "seq1\t274\tT\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\n")
write(f, "seq1\t275\tA\t23\t,....,,.,.,...,,,.,...^l.\t<+;9*<<<<<<<<<=<<:;<<<<\t23\t,....,,.,.,...,,,.,...^l.\t<+;9*<<<<<<<<<=<<:;<<<<\n")
write(f, "seq2\t276\tG\t22\t...T,,.,.,...,,,.,....\t33;+<<7=7<<7<&<<1;<<6<\t22\t...T,,.,.,...,,,.,....\t33;+<<7=7<<7<&<<1;<<6<\n")
write(f, "seq2\t277\tT\t22\t....,,.,.,.C.,,,.,..G.\t+7<;<<<<<<<&<=<<:;<<&<\t22\t....,,.,.,.C.,,,.,..G.\t+7<;<<<<<<<&<=<<:;<<&<\n")
write(f, "seq2\t278\tG\t23\t....,,.,.,...,,,.,....^k.\t%38*<<;<7<<7<=<<<;<<<<<\t23\t....,,.,.,...,,,.,....^k.\t%38*<<;<7<<7<=<<<;<<<<<\n")
write(f, "seq2\t279\tC\t23\tA..T,,.,.,...,,,.,.....\t75&<<<<<<<<<=<<<9<<:<<<\t23\tA..T,,.,.,...,,,.,.....\t75&<<<<<<<<<=<<<9<<:<<<\n")
close(f)
@time poolgen.pileup2syncx(pileup)
```

2. Multi-core parsing
```julia
using Distributed
Distributed.addprocs(length(Sys.cpu_info())-1)
@everywhere using poolgen
pileup = "test.pileup"
f = open(pileup, "a")
for i in 1:10_000
    write(f, "seq1\t272\tT\t24\t,......,,.,.,...,,,.,..^+.\t<<<+;<<<<<<<<<<<=<;<;7<&\t24\t,......,,.,.,...,,,.,..^+.\t<<<+;<<<<<<<<<<<=<;<;7<&\n")
    write(f, "seq1\t273\tT\t23\t,.....,,.,.,...,,,.,..A\t<<<;<<<<<<<<<3<=<<<;<<+\t23\t,.....,,.,.,...,,,.,..A\t<<<;<<<<<<<<<3<=<<<;<<+\n")
    write(f, "seq1\t274\tT\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\n")
end
close(f)
@time poolgen.pileup2syncx(pileup)
```
"""
function pileup2syncx(pileup::String; out::String="")::String
    ### Define output file if not specified
    if out == ""
        out = string(join(split(pileup, ".")[1:(end-1)], "."), ".syncx")
    end
    ### Find file positions for parallel processing
    threads = length(Distributed.workers())
    positions_init, positions_term = SPLIT(threads, pileup)
    digit = length(string(length(positions_init)))
    ### Parse to convert from pileup to syncx format
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        id = lpad(i, (digit+1)-length(string(i)), "0")
        tmp = string(pileup, "-PILEUP2SYNCX-", id, ".tmp")
        filename = PILEUP2SYNCX(pileup, init, term, tmp)
        [filename]
    end
    ### Sort the output files from parallel processing and merge
    MERGE(filenames_out, out)
    return(out)
end

"""
# ___________________________________
# Filter pileup

`filter(pileup::String, outype::String; maximum_missing_fraction::Float64=0.10, alpha1::Float64=0.05, maf::Float64=0.01, alpha2::Float64=0.50, minimum_coverage::Int64=5, out::String="")::String`

# Inputs
1. pileup [String]: pileup file
2. outype [String]: type of the output file, i.e. "pileup" or "syncx"
3. maximum_missing_fraction [Float64]: maximum fraction of individuals with missing loci (dafault=0.10)
4. alpha1 [Float64]: percentile allele frequency to implement `maf` (default=0.05)
5. maf [Float64]: minimum allele frequency (default=0.01)
6. alpha2 [Float64]: percentile coverage to implement `minimum_coverage` (default=0.50)
7. minimum_coverage [Int64]: minmum sequence coverage (default=5)
8. out [String]: output filename

# Output
1. [String]: filename of the output

# Examples
1. Single-core filtering
```julia
using poolgen
pileup = "test.pileup"
f = open(pileup, "a")
write(f, "seq1\t272\tT\t24\t,......,,.,.,...,,,.,..^+.\t<<<+;<<<<<<<<<<<=<;<;7<&\t24\t,......,,.,.,...,,,.,..^+.\t<<<+;<<<<<<<<<<<=<;<;7<&\n")
write(f, "seq1\t273\tT\t23\t,.....,,.,.,...,,,.,..A\t<<<;<<<<<<<<<3<=<<<;<<+\t23\t,.....,,.,.,...,,,.,..A\t<<<;<<<<<<<<<3<=<<<;<<+\n")
write(f, "seq1\t274\tT\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\n")
write(f, "seq1\t275\tA\t23\t,....,,.,.,...,,,.,...^l.\t<+;9*<<<<<<<<<=<<:;<<<<\t23\t,....,,.,.,...,,,.,...^l.\t<+;9*<<<<<<<<<=<<:;<<<<\n")
write(f, "seq2\t276\tG\t22\t...T,,.,.,...,,,.,....\t33;+<<7=7<<7<&<<1;<<6<\t22\t...T,,.,.,...,,,.,....\t33;+<<7=7<<7<&<<1;<<6<\n")
write(f, "seq2\t277\tT\t22\t....,,.,.,.C.,,,.,..G.\t+7<;<<<<<<<&<=<<:;<<&<\t22\t....,,.,.,.C.,,,.,..G.\t+7<;<<<<<<<&<=<<:;<<&<\n")
write(f, "seq2\t278\tG\t23\t....,,.,.,...,,,.,....^k.\t%38*<<;<7<<7<=<<<;<<<<<\t23\t....,,.,.,...,,,.,....^k.\t%38*<<;<7<<7<=<<<;<<<<<\n")
write(f, "seq2\t279\tC\t23\tA..T,,.,.,...,,,.,.....\t75&<<<<<<<<<=<<<9<<:<<<\t23\tA..T,,.,.,...,,,.,.....\t75&<<<<<<<<<=<<<9<<:<<<\n")
close(f)
@time poolgen.filter(pileup, "pileup")
@time poolgen.filter(pileup, "syncx")
```

2. Multi-core filtering
```julia
using Distributed
Distributed.addprocs(length(Sys.cpu_info())-1)
@everywhere using poolgen
pileup = "test.pileup"
f = open(pileup, "a")
for i in 1:10_000
    write(f, "seq1\t272\tT\t24\t,......,,.,.,...,,,.,..^+.\t<<<+;<<<<<<<<<<<=<;<;7<&\t24\t,......,,.,.,...,,,.,..^+.\t<<<+;<<<<<<<<<<<=<;<;7<&\n")
    write(f, "seq1\t273\tT\t23\t,.....,,.,.,...,,,.,..A\t<<<;<<<<<<<<<3<=<<<;<<+\t23\t,.....,,.,.,...,,,.,..A\t<<<;<<<<<<<<<3<=<<<;<<+\n")
    write(f, "seq1\t274\tT\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\n")
end
close(f)
@time poolgen.filter(pileup, "pileup")
@time poolgen.filter(pileup, "syncx")
```
"""
function filter(pileup::String, outype::String; maximum_missing_fraction::Float64=0.10, alpha1::Float64=0.05, maf::Float64=0.01, alpha2::Float64=0.50, minimum_coverage::Int64=5, out::String="")::String
    # using Distributed
    # Distributed.addprocs(length(Sys.cpu_info())-1)
    # @everywhere using ProgressMeter
    # @everywhere include("functions.jl")
    # @everywhere using .functions: PileupLine, LocusAlleleCounts, Window
    # @everywhere using .functions: SPLIT, MERGE, FILTER
    # pileup = "/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup"    
    # # pileup = "/home/jeffersonfparil/Documents/poolgen/test/test_2.pileup"    
    # outype = ["pileup", "syncx"][1]
    # maximum_missing_fraction = 0.10
    # alpha1 = 0.05
    # maf = 0.001
    # alpha2 = 0.50
    # minimum_coverage = 1
    # out = ""
    ### Define output file if not specified
    if out == ""
        if outype == "pileup"
            out = string(join(split(pileup, ".")[1:(end-1)], "."), "-FILTERED-missing_", maximum_missing_fraction, "-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, ".pileup")
        elseif outype == "syncx"
            out = string(join(split(pileup, ".")[1:(end-1)], "."), "-FILTERED-missing_", maximum_missing_fraction, "-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, ".syncx")
        end
    end
    ### Find file positions for parallel processing
    threads = length(Distributed.workers())
    positions_init, positions_term = SPLIT(threads, pileup)
    digit = length(string(length(positions_init)))
    ### Filter loci by minimum allele frequency (alpha1 and maf) and minimum coverage (alpha2 and minimum_coverage)
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        id = lpad(i, (digit+1)-length(string(i)), "0")
        tmp = string(pileup, "-FILTERED-", id, ".tmp")
        filename = FILTER(pileup, outype, init, term, maximum_missing_fraction, alpha1, maf, alpha2, minimum_coverage, tmp)
        [filename]
    end
    ### Sort the output files from parallel processing and merge
    MERGE(filenames_out, out)
    return(out)
end

"""
# ___________________________________
# Filter syncx

`filter(syncx::String; maximum_missing_fraction::Float64=0.10, alpha1::Float64=0.05, maf::Float64=0.01, alpha2::Float64=0.50, minimum_coverage::Int64=5, out::String="")::String`

# Inputs
1. syncx [String]: extended synchronised pileup file
2. maximum_missing_fraction [Float64]: maximum fraction of individuals with missing loci (dafault=0.10)
3. alpha1 [Float64]: percentile allele frequency to implement `maf` (default=0.05)
4. maf [Float64]: minimum allele frequency (default=0.01)
5. alpha2 [Float64]: percentile coverage to implement `minimum_coverage` (default=0.50)
6. minimum_coverage [Int64]: minmum sequence coverage (default=5)
7. out [String]: output filename

# Output
1. [String]: filename of the output

# Examples
1. Single-core conversion
```julia
using poolgen
n=5; m=10_000; l=135_000_000; k=5; ϵ=Int(1e+15); a=2; vec_chr_lengths=[0]; vec_chr_names=[""]; dist_noLD=500_000; o=100; t=10; nQTL=10; heritability=0.5; LD_chr=""; LD_n_pairs=10_000; plot_LD=false; npools=5
map, bim, ped, fam, syncx, csv = poolgen.simulate(n=n, m=m, l=l, k=k, ϵ=ϵ, a=a, vec_chr_lengths=vec_chr_lengths, vec_chr_names=vec_chr_names, dist_noLD=dist_noLD, o=o, t=t, nQTL=nQTL, heritability=heritability, npools=npools, LD_chr=LD_chr, LD_n_pairs=LD_n_pairs, plot_LD=plot_LD)
@time poolgen.filter(syncx)
```

2. Multi-core conversion
```julia
using Distributed
Distributed.addprocs(length(Sys.cpu_info())-1)
@everywhere using poolgen
n=5; m=10_000; l=135_000_000; k=5; ϵ=Int(1e+15); a=2; vec_chr_lengths=[0]; vec_chr_names=[""]; dist_noLD=500_000; o=100; t=10; nQTL=10; heritability=0.5; LD_chr=""; LD_n_pairs=10_000; plot_LD=false; npools=5
map, bim, ped, fam, syncx, csv = poolgen.simulate(n=n, m=m, l=l, k=k, ϵ=ϵ, a=a, vec_chr_lengths=vec_chr_lengths, vec_chr_names=vec_chr_names, dist_noLD=dist_noLD, o=o, t=t, nQTL=nQTL, heritability=heritability, npools=npools, LD_chr=LD_chr, LD_n_pairs=LD_n_pairs, plot_LD=plot_LD)
@time poolgen.filter(syncx)
```
"""
function filter(syncx::String; maximum_missing_fraction::Float64=0.10, alpha1::Float64=0.05, maf::Float64=0.01, alpha2::Float64=0.50, minimum_coverage::Int64=5, out::String="")::String
    # using Distributed
    # Distributed.addprocs(length(Sys.cpu_info())-1)
    # @everywhere using ProgressMeter
    # @everywhere include("structs.jl")
    # @everywhere include("functions.jl")
    # @everywhere include("user_functions.jl")
    # @everywhere using .structs: PileupLine, LocusAlleleCounts, Window
    # @everywhere using .functions: SPLIT, MERGE, FILTER
    # @everywhere using .user_functions: pileup2syncx
    # # syncx = pileup2syncx("/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup")
    # syncx = pileup2syncx("/home/jeffersonfparil/Documents/poolgen/test/test_2.pileup")
    # maximum_missing_fraction = 0.90
    # alpha1 = 0.05
    # maf = 0.001
    # alpha2 = 0.50
    # minimum_coverage = 1
    # out = ""
    ### Define output file if not specified
    if out == ""
        out = string(join(split(syncx, ".")[1:(end-1)], "."), "-FILTERED-missing_", maximum_missing_fraction, "-alpha1_", alpha1, "-maf_", maf, "-alpha2_", alpha2, "-cov_", minimum_coverage, ".syncx")
    end
    ### Find file positions for parallel processing
    threads = length(Distributed.workers())
    positions_init, positions_term = SPLIT(threads, syncx)
    digit = length(string(length(positions_init)))
    ### Filter loci by minimum allele frequency (alpha1 and maf) and minimum coverage (alpha2 and minimum_coverage)
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        id = lpad(i, (digit+1)-length(string(i)), "0")
        tmp = string(syncx, "-FILTERED-", id, ".tmp")
        filename = FILTER(syncx, init, term, maximum_missing_fraction, alpha1, maf, alpha2, minimum_coverage, tmp)
        [filename]
    end
    ### Sort the output files from parallel processing and merge
    MERGE(filenames_out, out)
    return(out)
end

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################
### IMPUTATION ###
##################
"""
# ___________________________________
# Impute pileup or synchronised pileup formats

`impute(filename::String; window_size::Int=100, model::String=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance::Bool=true, out::String="")::String`

# Inputs
1. filename [String]: filename of the input pileup or syncx file
2. window_size [Int]: sliding window size (slides 1 position; default=100)
3. model [String]: linear regression model to use in imputation. Select from "Mean", "OLS", "RR", "LASSO", and "GLMNET" (default="OLS")
4. distance [Bool]: use distances as covariate (default=true)
5. out [String]: output filename

# Output
1. [String]: filename of the output

# Examples
1. Single-core conversion
```julia
using poolgen
pileup = "test.pileup"
f = open(pileup, "a")
for i in 1:10_000
    write(f, "seq1\t272\tT\t24\t,......,,.,.,...,,,.,..^+.\t<<<+;<<<<<<<<<<<=<;<;7<&\t24\t,......,,.,.,...,,,.,..^+.\t<<<+;<<<<<<<<<<<=<;<;7<&\n")
    write(f, "seq1\t273\tT\t23\t,.....,,.,.,...,,,.,..A\t<<<;<<<<<<<<<3<=<<<;<<+\t23\t,.....,,.,.,...,,,.,..A\t<<<;<<<<<<<<<3<=<<<;<<+\n")
    write(f, "seq1\t274\tT\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\n")
    write(f, "seq1\t275\tT\t0\t*\t*\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\n")
    write(f, "seq1\t276\tT\t0\t*\t*\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\n")
end
close(f)
syncx = poolgen.pileup2syncx(pileup)
@time poolgen.impute(pileup)
@time poolgen.impute(syncx)
```

2. Multi-core conversion
```julia
using Distributed
Distributed.addprocs(length(Sys.cpu_info())-1)
@everywhere using poolgen
pileup = "test.pileup"
f = open(pileup, "a")
for i in 1:10_000
    write(f, "seq1\t272\tT\t24\t,......,,.,.,...,,,.,..^+.\t<<<+;<<<<<<<<<<<=<;<;7<&\t24\t,......,,.,.,...,,,.,..^+.\t<<<+;<<<<<<<<<<<=<;<;7<&\n")
    write(f, "seq1\t273\tT\t23\t,.....,,.,.,...,,,.,..A\t<<<;<<<<<<<<<3<=<<<;<<+\t23\t,.....,,.,.,...,,,.,..A\t<<<;<<<<<<<<<3<=<<<;<<+\n")
    write(f, "seq1\t274\tT\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\n")
    write(f, "seq1\t275\tT\t0\t*\t*\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\n")
    write(f, "seq1\t276\tT\t0\t*\t*\t23\t,.....,,.,.,...,,,.,...\t7<7;<;<<<<<<<<<=<;<;<<6\n")
end
close(f)
syncx = poolgen.pileup2syncx(pileup)
@time poolgen.impute(pileup)
@time poolgen.impute(syncx)
```
"""
function impute(filename::String; window_size::Int=100, model::String=["Mean", "OLS", "RR", "LASSO", "GLMNET"][2], distance::Bool=true, out::String="")::String
    # using Distributed
    # Distributed.addprocs(length(Sys.cpu_info())-1)
    # @everywhere using ProgressMeter
    # @everywhere include("functions.jl")
    # @everywhere include("user_functions.jl")
    # @everywhere using .functions: PileupLine, SyncxLine, LocusAlleleCounts, Window
    # @everywhere using .functions: PARSE, SPLIT, MERGE, PILEUP2SYNCX, IMPUTE
    # @everywhere using .user_functions: pileup2syncx
    # filename = "/home/jeffersonfparil/Documents/poolgen/test/test_1.pileup"
    # # filename = "/home/jeffersonfparil/Documents/poolgen/test/test_2.pileup"
    # # filename = "/home/jeffersonfparil/Documents/poolgen/test/test_1.syncx"
    # # filename = "/home/jeffersonfparil/Documents/poolgen/test/test_2.syncx"
    # window_size = 20
    # # window_size = 100
    # model = "LASSO"
    # distance = true
    # out = ""
    ### Define output file if not specified
    if out == ""
        out = string(join(split(filename, '.')[1:(end-1)], '.'), "-IMPUTED-window_", window_size, "-model_", model, "-distance_", distance, ".syncx")
    end
    ### Define the full path to the input and output files since calling functions within @distributed loop will revert back to the root directory from where julia was executed from
    if dirname(filename) == ""
        filename = string(pwd(), "/", filename)
    end
    if dirname(out) == ""
        out = string(pwd(), "/", out)
    end
    ### Convert to syncx if the input file is in pileup format: output syncx - filename of syncx file
    syncx = try
        file = open(filename, "r")
        PARSE(SyncxLine(1, readline(file)))
        close(file)
        filename
    catch
        pileup2syncx(filename)
    end
    ### Find file positions for parallel processing
    threads = length(Distributed.workers())
    positions_init, positions_term = SPLIT(threads, syncx, window_size)
    idx = collect(1:length(positions_term))[positions_term .== maximum(positions_term)][1]
    positions_init = positions_init[1:idx]
    positions_term = positions_term[1:idx]
    digit = length(string(length(positions_init)))
    ### Impute
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        id = lpad(i, (digit+1)-length(string(i)), "0")
        tmp = string(syncx, "-IMPUTED-", id, ".syncx.tmp")
        filename = IMPUTE(syncx, init, term, window_size, model, distance, tmp)
        [filename]
    end
    ### Merge the chunks
    if length(filenames_out) > 1
        MERGE(filenames_out, window_size, out)
    else
        MERGE(filenames_out, out)
    end
    return(out)
end

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##################
### SIMULATION ###
##################
"""
# ___________________________________
# Simulation function

`simulate(;n::Int64, m::Int64, l::Int64, k::Int64, ϵ::Int64=Int(1e+15), a::Int64=2, vec_chr_lengths::Vector{Int64}=Int64.([0]), vec_chr_names::Vector{String}=[""], dist_noLD::Int64=10_000, o::Int64=1_000, t::Int64=10, nQTL::Int64=10, heritability::Float64=0.5, npools::Int64=5, LD_chr::String="", LD_n_pairs::Int64=10_000, plot_LD::Bool=true)::Tuple{String, String, String, String, String, String}`

# Inputs
1.  n [Int64]: number of heterozygote founders
2.  m [Int64]: number of loci
3.  l [Int64]: total genome size in base-pairs
4.  k [Int64]: number of chromosomes
5.  ϵ [Int64]: an arbitrarily large Int64 number to indicate no LD, i.e. the distance between the termini of 2 adjacent chromosomes is infinitely large LD-wise but we don't want to use Inf as it is not Int64 (default=Int(1e+15))
6.  a [Int64]: number of alleles per locus (default=2 or biallelic)
7.  vec_chr_lengths [Vector{Int64}]: lengths of chromosomes (default: equally/similarly sized)
8.  vec_chr_names [Vector{String}]: chromosome names (default: consecutive numbers, i.e. "1", "2", ..., string(k))
9.  dist_noLD [Int64]: distance in base-pairs a with LD decays to zero (default=10_000)
10. o [Int64]: constant number of individuals per simulated generation, which also correspond to the total number of output individuals simulated (default=1_000)
11. t [Int64]: number of random mating generations (default=10)
12. nQTL [Int64]: number of completely additive quantitative trait loci (default=10)
13. heritability [Float64]: broad-sense or narrow-sense heritability as only additive effects are simulated (default=0.5)
14. npools [Int64]: number of pools (default=5)
15. LD_chr [String]: chromosome to use in linkage disequillibrium (LD) estimation (default: first chromosome)
16. LD_n_pairs [Int64]: number of randomly sampled loci pairs to use in LD estimation (default=10_000)
17. plot_LD [Bool]: plot the LD decay (default=true)

# Outputs
1. [String]: filename of map file ([plink1.9 map format](https://www.cog-genomics.org/plink/1.9/formats#map))
2. [String]: filename of bim file ([plink1.9 bim format](https://www.cog-genomics.org/plink/1.9/formats#bim))
3. [String]: filename of ped file ([plink1.9 ped format](https://www.cog-genomics.org/plink/1.9/formats#ped))
4. [String]: filename of fam file ([plink1.9 fam format](https://www.cog-genomics.org/plink/1.9/formats#fam))
5. [String]: filename of pool genotype file (extended synchronised pileup "syncx" format)
6. [String]: filename of phenotype data (comma-separated file; with a header where column 1 refers to the pool IDs, and column 2 is the phenotype values)

# Example
```julia
using poolgen
n=5; m=10_000; l=135_000_000; k=5; ϵ=Int(1e+15); a=2; vec_chr_lengths=[0]; vec_chr_names=[""]; dist_noLD=500_000; o=100; t=10; nQTL=10; heritability=0.5; LD_chr=""; LD_n_pairs=10_000; plot_LD=false; npools=5
@time map, bim, ped, fam, syncx, csv = poolgen.simulate(n=n, m=m, l=l, k=k, ϵ=ϵ, a=a, vec_chr_lengths=vec_chr_lengths, vec_chr_names=vec_chr_names, dist_noLD=dist_noLD, o=o, t=t, nQTL=nQTL, heritability=heritability, npools=npools, LD_chr=LD_chr, LD_n_pairs=LD_n_pairs, plot_LD=plot_LD)
```
"""
function simulate(;n::Int64, m::Int64, l::Int64, k::Int64, ϵ::Int64=Int(1e+15), a::Int64=2, vec_chr_lengths::Vector{Int64}=Int64.([0]), vec_chr_names::Vector{String}=[""], dist_noLD::Int64=10_000, o::Int64=1_000, t::Int64=10, nQTL::Int64=10, heritability::Float64=0.5, npools::Int64=5, LD_chr::String="", LD_n_pairs::Int64=10_000, plot_LD::Bool=true)::Tuple{String, String, String, String, String, String}
    # n = 5                 ### number of founders
    # m = 10_000            ### number of loci
    # l = 135_000_000       ### total genome length
    # k = 5                 ### number of chromosomes
    # ϵ = Int(1e+15)        ### some arbitrarily large number to signify the distance at which LD is nil
    # a = 2                 ### number of alleles per locus
    # vec_chr_lengths = [0] ### chromosome lengths
    # vec_chr_names = [""]  ### chromosome names 
    # dist_noLD = 500_000   ### distance at which LD is nil (related to ϵ)
    # o = 100               ### total number of simulated individuals
    # t = 10                ### number of random mating constant population size generation to simulate
    # nQTL = 10             ### number of QTL to simulate
    # heritability = 0.5    ### narrow(broad)-sense heritability as only additive effects are simulated
    # LD_chr = ""           ### chromosome to calculate LD decay from
    # LD_n_pairs = 10_000   ### number of randomly sampled pairs of loci to calculate LD
    # plot_LD = true        ### simulated LD decay
    vec_chr, vec_pos, X, y, b = SIMULATE(n, m, l, k, ϵ, a, vec_chr_lengths, vec_chr_names, dist_noLD, o, t, nQTL, heritability, LD_chr, LD_n_pairs, plot_LD)
    G, p = POOL(X, y, npools)
    map, bim, ped, fam = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, X, y)
    syncx, csv = EXPORT_SIMULATED_DATA(vec_chr, vec_pos, G, p)
    return(map, bim, ped, fam, syncx, csv)
end

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#############################
### QUANTITATIVE GENETICS ###
#############################
"""
# ___________________________________
# Genome-wide additive allelic effects (alpha) association using Pool-seq data

`gwalpha(;syncx::String, py_phenotype::String, maf::Float64=0.001, penalty::Bool=true, out::String="")::String`

# Original implementation
Compatible with the original implementation published and reposited in:
- Fournier-Level A, Robin C, Balding DJ (2016). [GWAlpha: Genome-Wide estimation of additive effects (Alpha) based on trait quantile distribution from pool-sequencing experiments.](https://doi.org/10.1093/bioinformatics/btw805)
- https://github.com/aflevel/GWAlpha


# Inputs
1. syncx [String]: extended synchronised pileup file
2. py_phenotype [String]: Text with a ".py" extension:
```python
    Pheno_name='Phenotype Name';
    sig=0.06724693662723039; # standard deviation
    MIN=0.0; # minimum phenotype value
    MAX=0.424591738712776; # maximum phenotype value
    perc=[0.2,0.4,0.6,0.8];	 # cummulative pool sizes percentiles excluding the last pool
    q=[0.16,0.20,0.23,0.27,0.42];  # phenotype values corresponding to each percentile
```
3. maf [Float64]: minimum allele frequency (default=0.001)
4. penalty [Bool]: use the penalization for low allele frequency (default=true)
5. out [String]: output filename (default: `syncx` with the extension converted to `.tsv`)

# Output
1. [String]: filename of the output in tab-delimited format
- *Column 1*: chromosome names
- *Column 2*: position
- *Column 3*: allele name
- *Column 4*: 0 for skipped locus; and 1 for included locus
- *Column 5*: allele frequency
- *Column 6*: alpha or allele effect

# Example
```julia
using Distributed
Distributed.addprocs(length(Sys.cpu_info())-1)
@everywhere using poolgen
n=5; m=10_000; l=135_000_000; k=5; ϵ=Int(1e+15); a=2; vec_chr_lengths=[0]; vec_chr_names=[""]; dist_noLD=500_000; o=100; t=10; nQTL=10; heritability=0.5; LD_chr=""; LD_n_pairs=10_000; plot_LD=false; npools=5
map, bim, ped, fam, syncx, csv = poolgen.simulate(n=n, m=m, l=l, k=k, ϵ=ϵ, a=a, vec_chr_lengths=vec_chr_lengths, vec_chr_names=vec_chr_names, dist_noLD=dist_noLD, o=o, t=t, nQTL=nQTL, heritability=heritability, npools=npools, LD_chr=LD_chr, LD_n_pairs=LD_n_pairs, plot_LD=plot_LD)

write(f, "Pheno_name=\'Phenotype Name\';\n")
write(f, "sig=0.06724693662723039;\n")
write(f, "MIN=0.0;\n")
write(f, "MAX=0.424591738712776;\n")
write(f, "perc=[0.2,0.4,0.6,0.8];\n")
write(f, "q=[0.16,0.20,0.23,0.27,0.42];\n")
```
"""
function gwalpha(;syncx::String, py_phenotype::String, maf::Float64=0.001, penalty::Bool=true, out::String="")::String
    # using Distributed
    # Distributed.addprocs(length(Sys.cpu_info())-1)
    # @everywhere using ProgressMeter
    # @everywhere include("functions.jl")
    # @everywhere include("user_functions.jl")
    # @everywhere using .functions: PileupLine, SyncxLine, LocusAlleleCounts, Window
    # @everywhere using .functions: SPLIT, MERGE, GWALPHA
    # syncx = "/home/jeffersonfparil/Documents/poolgen/test/test_3.syncx"
    # py_phenotype = "/home/jeffersonfparil/Documents/poolgen/test/test_3-pheno-any-filename.py"
    # maf = 0.001
    # penalty = false
    # out = ""
    ### Define output file if not specified
    if out == ""
        out = string(join(split(syncx, '.')[1:(end-1)], '.'), "-GWAlpha-maf_", round(maf, digits=5), ".tsv")
    end
    ### Define the full path to the input and output files since calling functions within @distributed loop will revert back to the root directory from where julia was executed from
    if dirname(syncx) == ""
        syncx = string(pwd(), "/", syncx)
    end
    if dirname(out) == ""
        out = string(pwd(), "/", out)
    end
    ### Find file positions for parallel processing
    threads = length(Distributed.workers())
    positions_init, positions_term = SPLIT(threads, syncx)
    digit = length(string(length(positions_init)))
    ### Impute
    @time filenames_out = @sync @showprogress @distributed (append!) for i in eachindex(positions_init)
        init = positions_init[i]
        term = positions_term[i]
        id = lpad(i, (digit+1)-length(string(i)), "0")
        tmp = string(syncx, "-GWAlpha-", id, ".tsv.tmp")
        filename = GWALPHA(syncx, py_phenotype, init, term, maf, penalty, tmp)
        [filename]
    end
    ### Merge the chunks
    MERGE(filenames_out, out)
    return(out)
end

"""
# ___________________________________
# Fit genomic prediction models

`genomic_prediction(;syncx::String, maf::Float64, phenotype::String, model::String=["OLS", "ELASTIC", "LMM"][1], delimiter::String=",", header::Bool=true, id_col::Int=1, phenotype_col::Int=2, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], FE_method::String=["CANONICAL", "N<<P"][2], alpha::Float64=1.0, covariate::String=["", "XTX", "COR"][2], MM_model::String=["GBLUP", "RRBLUP"][1], MM_method::String=["ML", "REML"][1], inner_optimizer=["LBFGS", "BFGS", "SimulatedAnnealing", "GradientDescent", "NelderMead"][1], optim_trace::Bool=false, out::String="")::String`

# Inputs
1.  syncx [String]: extended synchronised pileup file
2.  phenotype [String]: phenotype data (comma-separated file; with a header where column 1 refers to the pool IDs, and column 2 is the phenotype values)
3.  model [String]: genomic prediction model. Choose from "OLS", "ELASTIC", and "LMM"] (default="OLS")
4.  maf [Float64]: minimum allele frequency (default=0.001)
5.  delimiter [String]: delimited of the `phenotype` (default=",")
6.  header [Bool]: header of the `phenotype` (default=true)
7.  id_col [Int]: column of the `phenotype` containing the pool IDs (default=1)
8.  phenotype_col [Int]: column of the `phenotype` containing the phenotype values (default=2)
9.  missing_strings [Vector{String}]: missing phenotype data encoding (default=["NA", "NAN", "NaN", "missing", ""])
10. FE_method [String]: fixed effect estimation method. Choose from "CANONICAL" (inverse(XᵀX)(Xᵀy)), and "N<<P" ((Xᵀ)inverse(XXᵀ)(y)) (default="N<<P")
11. alpha [Float64]: elastic-net penalty. Ranges from 0.0 (ridge) to 1.0 (lasso) (default=1.0)
12. covariate [String]: covarate to use in linear mixed models. Choose from "" (none), "XTX" (unscaled kinship matrix or squared Euclidean distance between indvidual), and "COR" (Pearson's correlation matrix) (default: "XTX")
13. MM_model [String]: linear mixed model. Chose from "GBLUP", and "RRBLUP"][1]
    - GBLUP: `y = Xβ + g + ϵ`,
        where `X` are the intercept and SNPs,
        `g = Zμ = μ` are the individual genotype effects,
            where `Z = I(nxn)`, and `(μ==g)~MVN(0, D),`
                where `D = σ2u * K`
                    where `K ≈ (X'X)/n`
        and `ϵ~MVN(0, R)`
            where `R = σ2e * I`
    - RR-BLUP: `y = Xβ + Zμ + ϵ`,
        where `X` are the intercept and covariates, if any,
        `Z` are the SNPs,
        μ~MVN(0, D),
            where D = σ2u * I
        and ϵ~MVN(0, R)
            where R = σ2e * I
14. MM_method [String]: linear mixed model parameter estimation method. Choose from "ML" (maximum likelihood), and "REML" (restricted maximum likelihood) (default="ML")
15. inner_optimizer [String]: linear mixed model parameter estimation inner optimiser. Choose from "LBFGS", "BFGS", "SimulatedAnnealing", "GradientDescent", and "NelderMead" (default="LBFGS")
16. optim_trace [Bool]: print out optimisation progress (default=false)
17. out [String]: output filename (default: `syncx` with the extension converted to `.tsv`)

# Output
1. [String]: filename of the output in tab-delimited format
- *Column 1*: chromosome names
- *Column 2*: position
- *Column 3*: allele name
- *Column 4*: allele frequency
- *Column 5*: allele effect
- *Column 6*: "NA" (used as p-value column in GWAS)

# Example
```julia
using Distributed
Distributed.addprocs(length(Sys.cpu_info())-1)
@everywhere using poolgen
n=5; m=10_000; l=135_000_000; k=5; ϵ=Int(1e+15); a=2; vec_chr_lengths=[0]; vec_chr_names=[""]; dist_noLD=500_000; o=100; t=10; nQTL=10; heritability=0.5; LD_chr=""; LD_n_pairs=10_000; plot_LD=false; npools=5
map, bim, ped, fam, syncx, csv = poolgen.simulate(n=n, m=m, l=l, k=k, ϵ=ϵ, a=a, vec_chr_lengths=vec_chr_lengths, vec_chr_names=vec_chr_names, dist_noLD=dist_noLD, o=o, t=t, nQTL=nQTL, heritability=heritability, npools=npools, LD_chr=LD_chr, LD_n_pairs=LD_n_pairs, plot_LD=plot_LD)
@time poolgen.genomic_prediction(syncx=syncx, phenotype=csv)
@time poolgen.genomic_prediction(syncx=syncx, phenotype=csv, model="ELASTIC", alpha=0.50)
@time poolgen.genomic_prediction(syncx=syncx, phenotype=csv, model="LMM", MM_model="GBLUP", MM_method="ML")
```
"""
function genomic_prediction(;syncx::String, phenotype::String, model::String=["OLS", "ELASTIC", "LMM"][1], maf::Float64=0.001, delimiter::String=",", header::Bool=true, id_col::Int=1, phenotype_col::Int=2, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], FE_method::String=["CANONICAL", "N<<P"][2], alpha::Float64=1.0, covariate::String=["", "XTX", "COR"][2], MM_model::String=["GBLUP", "RRBLUP"][1], MM_method::String=["ML", "REML"][1], inner_optimizer::String=["LBFGS", "BFGS", "SimulatedAnnealing", "GradientDescent", "NelderMead"][1], optim_trace::Bool=false, out::String="")::String
    # model = ["OLS", "ELASTIC", "LMM"][1]
    # syncx = "../test/test_5.syncx"
    # maf = 0.001
    # phenotype = "../test/test_5.csv"
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # FE_method = ["CANONICAL", "N<<P"][2]
    # alpha = 1.0
    # covariate = ["", "XTX", "COR"][2]
    # MM_model = ["GBLUP", "RRBLUP"][1]
    # MM_method = ["ML", "REML"][1]
    # inner_optimizer=["LBFGS", "BFGS", "SimulatedAnnealing", "GradientDescent", "NelderMead"][1]
    # optim_trace = false
    # out = ""
    # out = genomic_prediction(model=model, syncx=syncx, maf=maf, phenotype=phenotype, delimiter=delimiter, header=header, id_col=id_col, phenotype_col=phenotype_col, missing_strings=missing_strings, FE_method=FE_method, alpha=alpha, covariate=covariate, MM_model=MM_model, MM_method=MM_method, inner_optimizer=inner_optimizer, optim_trace=optim_trace, out=out)
    ### Fit
    if model == "OLS"
        out = OLS_MULTIVAR(syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, FE_method, out)
    elseif model == "ELASTIC"
        out = ELA_MULTIVAR(syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, alpha, out)
    elseif model == "LMM"
        out = LMM_MULTIVAR(syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings, covariate, MM_model, MM_method, FE_method, inner_optimizer, optim_trace, out)
    else
        println(string("Sorry the genomic prodection model: ", model, " is not implemented."))
        println("Please select from: ")
        println("\t- 'OLS': ordinary least squares")
        println("\t- 'ELASTIC': elastic-net penalised regression regression")
        println("\t- 'LMM': linear mixed model.")
    end
    return(out)
end

"""
# ___________________________________
# Genomic prediction cross-validation

`genomic_prediction_CV(;nrep::Int64, nfold::Int64, syncx::String, maf::Float64, phenotype::String, model::String=["OLS", "ELASTIC", "LMM"][1], delimiter::String=",", header::Bool=true, id_col::Int=1, phenotype_col::Int=2, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], FE_method::String=["CANONICAL", "N<<P"][2], alpha::Float64=1.0, covariate::String=["", "XTX", "COR"][2], MM_model::String=["GBLUP", "RRBLUP"][1], MM_method::String=["ML", "REML"][1], inner_optimizer::String=["LBFGS", "BFGS", "SimulatedAnnealing", "GradientDescent", "NelderMead"][1], optim_trace::Bool=false, save_plots::Bool=false, save_predictions::Bool=false, out::String="")::String`

# Inputs
1.  nrep [Int64]: number of randomisation replication to perform k-fold cross-validation on
2.  nfold [Int64]: number of sets to divide the pools into for k-fold cross-validation. If this value results to less than 5 pools per set, the a lower value will be selected.
3.  syncx [String]: extended synchronised pileup file
4.  phenotype [String]: phenotype data (comma-separated file; with a header where column 1 refers to the pool IDs, and column 2 is the phenotype values)
5.  model [String]: genomic prediction model. Choose from "OLS", "ELASTIC", and "LMM"] (default="OLS")
6.  maf [Float64]: minimum allele frequency (default=0.001)
7.  delimiter [String]: delimited of the `phenotype` (default=",")
8.  header [Bool]: header of the `phenotype` (default=true)
9.  id_col [Int]: column of the `phenotype` containing the pool IDs (default=1)
10. phenotype_col [Int]: column of the `phenotype` containing the phenotype values (default=2)
11. missing_strings [Vector{String}]: missing phenotype data encoding (default=["NA", "NAN", "NaN", "missing", ""])
12. FE_method [String]: fixed effect estimation method. Choose from "CANONICAL" (inverse(XᵀX)(Xᵀy)), and "N<<P" ((Xᵀ)inverse(XXᵀ)(y)) (default="N<<P")
13. alpha [Float64]: elastic-net penalty. Ranges from 0.0 (ridge) to 1.0 (lasso) (default=1.0)
14. covariate [String]: covarate to use in linear mixed models. Choose from "" (none), "XTX" (unscaled kinship matrix or squared Euclidean distance between indvidual), and "COR" (Pearson's correlation matrix) (default: "XTX")
15. MM_model [String]: linear mixed model. Chose from "GBLUP", and "RRBLUP"][1]
    - GBLUP: `y = Xβ + g + ϵ`,
        where `X` are the intercept and SNPs,
        `g = Zμ = μ` are the individual genotype effects,
            where `Z = I(nxn)`, and `(μ==g)~MVN(0, D),`
                where `D = σ2u * K`
                    where `K ≈ (X'X)/n`
        and `ϵ~MVN(0, R)`
            where `R = σ2e * I`
    - RR-BLUP: `y = Xβ + Zμ + ϵ`,
        where `X` are the intercept and covariates, if any,
        `Z` are the SNPs,
        μ~MVN(0, D),
            where D = σ2u * I
        and ϵ~MVN(0, R)
            where R = σ2e * I
16. MM_method [String]: linear mixed model parameter estimation method. Choose from "ML" (maximum likelihood), and "REML" (restricted maximum likelihood) (default="ML")
17. inner_optimizer [String]: linear mixed model parameter estimation inner optimiser. Choose from "LBFGS", "BFGS", "SimulatedAnnealing", "GradientDescent", and "NelderMead" (default="LBFGS")
18. optim_trace [Bool]: print out optimisation progress (default=false)
19. save_plots [Bool]: save scatter plots of true vs. predicted phenotypes (default=false)
20. save_predictions [Bool]: save true and predicted phenotypes across replications and folds (default=false)
21. out [String]: output filename (default: `syncx` with the extension converted to `-MULTIVAR_CV.tsv`)

# Output
1. [String]: filename of the genomic prediction cross-validation accuracy metrics in tab-delimited format (header: "rep", "fold", "correlation_pearson", "correlation_spearman", "correlation_kendall", "R2", "R2_adj", "MAE", "MBE", "RAE", "MSE", "RMSE", "RRMSE", "RMSLE")
- *Column 1*:  replication number
- *Column 2*:  fold number
- *Column 3*:  coefficient of determination (model = predicted ~ 1 + true)
- *Column 4*:  adjusted coefficient of determination
- *Column 5*:  mean absolute error (MAE = mean(abs.(y .- ŷ)))
- *Column 6*:  mean bias error (MBE = mean(y .- ŷ))
- *Column 7*:  relative absolute error (RAE = sum(abs.(y .- ŷ)) / sum(abs.(y .- mean(y))))
- *Column 8*:  mean square error (MSE = mean((y .- ŷ).^2))
- *Column 9*:  root mean square error (RMSE = sqrt(MSE))
- *Column 10*: relative root mean square error (RRMSE = sqrt(MSE / sum(ŷ.^2)))
- *Column 11*: root mean squared logarithmic error (RMSLE = sqrt(mean((log.(y .+ min_y .+ 1) .- log10.(ŷ .+ min_ŷ .+ 1)).^2)))
2. [String; Optional output, i.e. if save_predictions==true]: filename of the true and predicted phenotype values across replications and folds (header: "rep", "fold", "true", "pred")
- *Column 1*:  replication number
- *Column 2*:  fold number
- *Column 3*:  true phenotype value
- *Column 4*:  predicted phenotype value

# Example
```julia
using Distributed
Distributed.addprocs(length(Sys.cpu_info())-1)
@everywhere using poolgen
n=5; m=10_000; l=135_000_000; k=5; ϵ=Int(1e+15); a=2; vec_chr_lengths=[0]; vec_chr_names=[""]; dist_noLD=500_000; o=100; t=10; nQTL=10; heritability=0.5; LD_chr=""; LD_n_pairs=10_000; plot_LD=false; npools=5
map, bim, ped, fam, syncx, csv = poolgen.simulate(n=n, m=m, l=l, k=k, ϵ=ϵ, a=a, vec_chr_lengths=vec_chr_lengths, vec_chr_names=vec_chr_names, dist_noLD=dist_noLD, o=o, t=t, nQTL=nQTL, heritability=heritability, npools=npools, LD_chr=LD_chr, LD_n_pairs=LD_n_pairs, plot_LD=plot_LD)
@time poolgen.genomic_prediction(nrep=3, nfold=10, syncx=syncx, phenotype=csv)
@time poolgen.genomic_prediction(nrep=3, nfold=10, syncx=syncx, phenotype=csv, model="ELASTIC", alpha=0.50)
@time poolgen.genomic_prediction(nrep=3, nfold=10, syncx=syncx, phenotype=csv, model="LMM", MM_model="GBLUP", MM_method="ML")
```
"""
function genomic_prediction_CV(;nfold::Int64, nrep::Int64, syncx::String, phenotype::String, model::String=["OLS", "ELASTIC", "LMM"][1], maf::Float64=0.001, delimiter::String=",", header::Bool=true, id_col::Int=1, phenotype_col::Int=2, missing_strings::Vector{String}=["NA", "NAN", "NaN", "missing", ""], FE_method::String=["CANONICAL", "N<<P"][2], alpha::Float64=1.0, covariate::String=["", "XTX", "COR"][2], MM_model::String=["GBLUP", "RRBLUP"][1], MM_method::String=["ML", "REML"][1], inner_optimizer::String=["LBFGS", "BFGS", "SimulatedAnnealing", "GradientDescent", "NelderMead"][1], optim_trace::Bool=false, save_plots::Bool=false, save_predictions::Bool=false, out::String="")::String
    # nfold = 10
    # nrep = 3
    # model = ["OLS", "ELASTIC", "LMM"][1]
    # syncx = "../test/test_5.syncx"
    # maf = 0.001
    # phenotype = "../test/test_5.csv"
    # delimiter = ","
    # header = true
    # id_col = 1
    # phenotype_col = 2
    # missing_strings = ["NA", "NAN", "NaN", "missing", ""]
    # FE_method = ["CANONICAL", "N<<P"][2]
    # alpha = 1.0
    # covariate = ["", "XTX", "COR"][2]
    # MM_model = ["GBLUP", "RRBLUP"][1]
    # MM_method = ["ML", "REML"][1]
    # inner_optimizer=["LBFGS", "BFGS", "SimulatedAnnealing", "GradientDescent", "NelderMead"][1]
    # optim_trace = false
    # save_plots = false
    # save_predictions = false
    # out = ""
    # out = poolgen.genomic_prediction_CV(nfold=nfold, nrep=nrep, model=model, syncx=syncx, maf=maf, phenotype=phenotype, delimiter=delimiter, header=header, id_col=id_col, phenotype_col=phenotype_col, missing_strings=missing_strings, FE_method=FE_method, alpha=alpha, covariate=covariate, MM_model=MM_model, MM_method=MM_method, inner_optimizer=inner_optimizer, optim_trace=optim_trace, out=out)
    ### Fit
    if model == "OLS"
        params=[FE_method]
        out = CV_MULTIVAR(nfold, nrep, syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings,
                          OLS_MULTIVAR,
                          params,
                          save_plots,
                          save_predictions,
                          out)
    elseif model == "ELASTIC"
        params=[alpha]
        out = CV_MULTIVAR(nfold, nrep, syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings,
                          ELA_MULTIVAR,
                          params,
                          save_plots,
                          save_predictions,
                          out)
    elseif model == "LMM"
        params=[covariate, MM_model, MM_method, FE_method, inner_optimizer, optim_trace]
        out = CV_MULTIVAR(nfold, nrep, syncx, maf, phenotype, delimiter, header, id_col, phenotype_col, missing_strings,
                          LMM_MULTIVAR,
                          params,
                          save_plots,
                          save_predictions,
                          out)
    else
        println(string("Sorry the genomic prodection model: ", model, " is not implemented."))
        println("Please select from: ")
        println("\t- 'OLS': ordinary least squares")
        println("\t- 'ELASTIC': elastic-net penalised regression regression")
        println("\t- 'LMM': linear mixed model.")
    end
    return(out)
end


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
###########################
### POPULATION GENETICS ###
###########################



end