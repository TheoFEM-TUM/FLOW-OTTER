using SparseArrays, LinearAlgebra, DelimitedFiles

include("read_hamiltonian.jl")


TB_path = ARGS[1]
output_path = ARGS[2]
t = parse(Int, ARGS[3])
hamiltonian_style = ARGS[4]  


H = get_dense_H(TB_path, t, hamiltonian_style)

E, _ = eigen(H)

open(joinpath(output_path, "EV_$(t).txt"), "w") do file
    write(file, "# eigenvalues of hamiltonian \n")
    for i in axes(E, 1)
        write(file, "$(E[i])" * "\n")
    end
end
