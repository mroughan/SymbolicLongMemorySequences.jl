using Pkg

Pkg.activate(; temp = true)
Pkg.add(url = "https://github.com/mroughan/IncCSV.jl")
Pkg.develop(path = dirname(@__DIR__))
Pkg.add("JET")

using JET
using SymbolicLongMemorySequences

JET.test_package(SymbolicLongMemorySequences; target_modules = (SymbolicLongMemorySequences,))

