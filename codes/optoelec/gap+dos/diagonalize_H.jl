using SparseArrays, LinearAlgebra, DelimitedFiles
using MPI

include("read_H.jl")

TB_path = ARGS[1]
output_path = ARGS[2]
t = parse(Int, ARGS[3])
hamiltonian_style = ARGS[4]  

MPI.Init()
H = get_dense_H(TB_path, t, hamiltonian_style)
MPI.Finalize()

E, _ = eigen(H)

eps = 10^(-7)

if all(abs.(imag.(E)) .< eps)

    open(joinpath(output_path, "EV_$(t).txt"), "w") do file
        write(file, "# eigenvalues of hamiltonian \n")
        for i in axes(E, 1)
            write(file, "$(real(E[i]))" * "\n")
        end
    end

else

    open(joinpath(output_path, "EV_$(t).txt"), "w") do file
        write(file, "# eigenvalues of hamiltonian \n")
        for i in axes(E, 1)
            write(file, "$(E[i])" * "\n")
        end
    end

    error("Error: Eigenvalues are not purely real!")
end
