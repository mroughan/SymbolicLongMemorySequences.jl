using Pkg

Pkg.activate(; temp = true)
Pkg.add(url = "https://github.com/mroughan/IncCSV.jl")
Pkg.develop(path = dirname(@__DIR__))
Pkg.add("Aqua")

using Aqua
using SymbolicLongMemorySequences

Aqua.test_all(SymbolicLongMemorySequences; stale_deps = (ignore = [:StableRNGs],))

