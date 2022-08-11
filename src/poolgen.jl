module poolgen

include("user_functions.jl")
using .user_functions: pileup2syncx, filter, impute

### Documentation
"""
# ____________________________________________________________________
# poolgen
# Usage
`poolgen.impute():::String`

# Inputs
1. 
# Output
syncx_imputed:
Syncx format (after popoolation2's sync or synchronised pileup file format):
- Column 1:   chromosome or scaffold name
- Column 2:   locus position repeated 7 times corresponding to alleles "A", "T", "C", "G", "INS", "DEL", "N", where "INS" is insertion, "DEL" is deletion, and "N" is unclassified
- Column 3-n: allele counts one column for each pool or population

# Examples
```
# Single-threaded execution
using poolgen
poolgen.impute("test.pileup")

# Multi-threaded execution
using Distributed
int_thread_count = length(Sys.cpu_info())-1
Distributed.addprocs(int_thread_count)
@everywhere using poolgen
poolgen.impute("test.pileup", window_size=20, threads=2, lines_per_chunk=30)
```

# Details

# Author
- Jeff Paril (jeffersonparil@gmail.com; https://orcid.org/0000-0002-5693-4123)
...
"""

end
