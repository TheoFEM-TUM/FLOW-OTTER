using Pkg

Pkg.add([
    "ArgParse",
    "DelimitedFiles",
    "Distributions",
    "HDF5",
    "KrylovKit",
    "LinearAlgebra",
    "MPI",
    "MPIPreferences",
    "OhMyThreads",
    "PhysicalConstants",
    "Printf",
    "Random",
    "SparseArrays",
    "SpecialFunctions",
    "Statistics",
    "Unitful",
    "YAML",
])


Pkg.add([
    "Vampires",
    "Hamster",
    "conductivity",
])