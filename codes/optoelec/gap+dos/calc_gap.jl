using SparseArrays, LinearAlgebra, KrylovKit, DelimitedFiles

include("read_H.jl")

TB_path = ARGS[1]
t = parse(Int, ARGS[2])
guess_E_v = parse(Float64, ARGS[3]) 
guess_E_c = parse(Float64, ARGS[4])
output_path = ARGS[5]
hamiltonian_style = ARGS[6]


H = get_sparse_H(TB_path, t, hamiltonian_style)

E_c = real(eigsolve(H, 1, EigSorter(λ->abs(guess_E_c-λ), rev=false), ishermitian=true)[1][1])
E_v = real(eigsolve(H, 1, EigSorter(λ->abs(guess_E_v-λ), rev=false), ishermitian=true)[1][1])

gap = abs(real(E_c - E_v))

open(joinpath(output_path, "gap_$(t)_KPM.txt"), "w") do file
    write(file, "# gap     VBM     CBM \n")
    write(file, "$gap $(real(E_v)) $(real(E_c)) \n")
end
