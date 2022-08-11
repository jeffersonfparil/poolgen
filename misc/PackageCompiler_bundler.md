# Create the standalone binary file or bundle here for ease of distribution and use to the general public

Using [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl), be sure you've downloaded julia from the [official Julia website](https://julialang.org/downloads/) and not through a package manager because there is a high chance that generating the standaone bundle won't work.

## Instantiate a new package:

```{julia}
cd("Documents/")
run(`git clone https://github.com/jeffersonfparil/PoPoolImpute.jl`)
using Pkg
Pkg.generate("PoPoolImpute")
cp("PoPoolImpute.jl/src/PoPoolImpute.jl", "PoPoolImpute/src/PoPoolImpute.jl", force=true)
cp("PoPoolImpute.jl/src/functions.jl", "PoPoolImpute/src/functions.jl")
cp("PoPoolImpute.jl/Project.toml", "PoPoolImpute/Project.toml", force=true)
```

## Open `PoPoolImpute/src/PoPoolImpute.jl` and add the following function after loading the libraries bu before the function definition and markdown documentation:

```{julia snippet}
function julia_main()::Cint
    try
        pileup_with_missing = ARGS[1]
        window_size = try
            parse(Int, ARGS[2])
        catch
            100
        end
        model = try
            ARGS[3]
        catch
            "RR"
        end
        distance = try
            parse(Bool, ARGS[4])
        catch
            false
        end
        syncx_imputed = try
            ARGS[5]
        catch
            ""
        end
        threads = try
            parse(Int, ARGS[6])
        catch
            2
        end
        lines_per_chunk = try
            parse(Int, ARGS[7])
        catch
            10_000
        end
        fname_output = impute(pileup_with_missing,
                              window_size=window_size,
                              model=model,
                              distance=distance,
                              syncx_imputed=syncx_imputed,
                              threads=threads,
                              lines_per_chunk=lines_per_chunk)
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end
```

## Build the standalone bundle:

```{julia}
cd("Documents/")
using PackageCompiler
create_app("PoPoolImpute", "PoPoolImpute-BUILD", force=true)
```

## Test the build:

### Simulate missing data:

```{julia}
using Pkg
cd("Documents/PoPoolImpute.jl/test")
run(`tar -xvf test.pileup.tar.xz`)
Pkg.add(url="https://github.com/jeffersonfparil/PoPoolImpute.jl.git")
using PoPoolImpute

pileup_without_missing = "test.pileup"
pileup_with_missing = PoPoolImpute.functions.SIMULATESPARSITY(pileup_without_missing,
                                                              read_length=10,
                                                              missing_loci_fraction=0.50,
                                                              missing_pools_fraction=0.25)
```

### Impute:

```{sh}
cd ~/Documents/PoPoolImpute-BUILD/bin
time \
./PoPoolImpute \
    ~/Documents/PoPoolImpute.jl/test/test-SIMULATED_MISSING.pileup \
    20 \
    RR \
    false \
    ~/test.syncx \
    2 \
    45
cd -
```

### Cross-validation:

```{julia}
using PoPoolImpute
using UnicodePlots

pileup_without_missing = "Documents/PoPoolImpute.jl/test/test.pileup"
pileup_with_missing = "Documents/PoPoolImpute.jl/test/test-SIMULATED_MISSING.pileup"
syncx_imputed = "test.syncx"

syncx_without_missing = PoPoolImpute.functions.PILEUP2SYNCX(pileup_without_missing)
syncx_with_missing = PoPoolImpute.functions.PILEUP2SYNCX(pileup_with_missing)
csv_accuracy = PoPoolImpute.functions.CROSSVALIDATE(syncx_without_missing,
                                                    syncx_with_missing,
                                                    syncx_imputed)
### plot
file = open(csv_accuracy)
expected = []
imputed = []
expected_freq = []
imputed_freq = []
fraction_of_missing_imputed = 0
while !eof(file)
    line = split(readline(file), ',')
    if line[2] != "" 
        push!(expected, parse(Float64, line[1]))
        push!(imputed, parse(Float64, line[2]))
        push!(expected_freq, parse(Float64, line[3]))
        push!(imputed_freq, parse(Float64, line[4]))
    else
        fraction_of_missing_imputed = parse(Float64, line[1])
    end
end
close(file)
plot1 = UnicodePlots.scatterplot(Int.(expected), Int.(imputed),
                                    title="Counts",
                                    grid=true, color=:white, canvas=BlockCanvas)
plot2 = UnicodePlots.scatterplot(Float64.(expected_freq), Float64.(imputed_freq),
                                    title="Frequencies",
                                    grid=true, color=:white, canvas=BlockCanvas)
@show plot1
@show plot2
@show RMSE_count = sqrt(sum((expected .- imputed).^2)/length(expected))
@show RMSE_freqs = sqrt(sum((expected_freq .- imputed_freq).^2)/length(expected_freq))
@show fraction_of_missing_imputed
```