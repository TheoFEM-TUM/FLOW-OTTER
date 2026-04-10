using LinearAlgebra, SparseArrays
using MPI
using DelimitedFiles
using Plots
include("../gap+dos/read_H.jl")

function _broadening(dE::Float64, sigma::Float64, kind::Symbol)::Float64
    if kind === :gaussian
        return exp(-0.5 * (dE / sigma)^2) / (sigma * sqrt(2π))
    elseif kind === :lorentzian
        return (sigma / π) / (dE^2 + sigma^2)
    else
        throw(ArgumentError("Unknown broadening kind: $kind. Use :gaussian or :lorentzian."))
    end
end

"""
    compute_cohp_pair(H, i_orb, j_orb; n_E, sigma, broadening)

Compute the pairwise COHP between orbitals `i_orb` and `j_orb`.

For each eigenstate k with eigenvector c and eigenvalue ε_k, the symmetrised
orbital-pair contribution is:

    COHP_ij(E) = Σ_k [ H[i,j]·c[i,k]·c[j,k] + H[j,i]·c[j,k]·c[i,k] ] · δ(E - ε_k)

For a Hermitian (real symmetric) H this equals 2·H[i,j]·c[i,k]·c[j,k], but
the symmetrised form is kept explicitly for clarity and generality.

Returns: E_grid, cohp, energies, vecs
"""
function compute_cohp_pair(H::AbstractMatrix,
                           i_orb::Int,
                           j_orb::Int;
                           n_E::Int       = 1000,
                           sigma::Float64 = 0.05,
                           broadening::Symbol = :gaussian)

    n_orb = size(H, 1)
    @assert 1 <= i_orb <= n_orb "i_orb=$i_orb out of range [1, $n_orb]"
    @assert 1 <= j_orb <= n_orb "j_orb=$j_orb out of range [1, $n_orb]"

    # Diagonalise — H is gamma-only so dense eigen is fine
    energies, vecs = eigen(Hermitian(Matrix(H)))   # real eigenvalues guaranteed

    E_grid = range(minimum(energies) * 1.05, maximum(energies) * 1.05, length=n_E)
    cohp   = zeros(Float64, n_E)

    H_ij = H[i_orb, j_orb]   # scalar, reused for every eigenstate
    H_ji = H[j_orb, i_orb]

    @inbounds for k in 1:n_orb
        ε_k  = energies[k]
        ci   = vecs[i_orb, k]
        cj   = vecs[j_orb, k]

        # Symmetrised pairwise COHP weight for this eigenstate
        cohp_k = real(H_ij * conj(ci) * cj + H_ji * conj(cj) * ci)

        for iE in 1:n_E
            cohp[iE] += cohp_k * _broadening(E_grid[iE] - ε_k, sigma, broadening)
        end
    end

    return collect(E_grid), cohp, energies, vecs
end

# ---------------------------------------------------------------------------
# Script arguments
# ---------------------------------------------------------------------------
# ARGS[1]  path to Hamiltonian file
# ARGS[2]  output directory (trailing slash recommended)
# ARGS[3]  t   (integer, passed to get_dense_H)
# ARGS[4]  hamiltonian_style
# ARGS[5]  i_orb  (1-based orbital index)
# ARGS[6]  j_orb  (1-based orbital index)

#length(ARGS) >= 6 || error("Usage: julia cohp_pair.jl <H_path> <output_path> <t> <style> <i_orb> <j_orb>")

H_path            = ARGS[1]
output_path       = ARGS[2]
t                 = parse(Int, ARGS[3])
hamiltonian_style = ARGS[4]
i_orb             = parse(Int, ARGS[5])
j_orb             = parse(Int, ARGS[6])

# ---------------------------------------------------------------------------
# Build Hamiltonian
# ---------------------------------------------------------------------------
MPI.Init()
H = get_dense_H(H_path, t, hamiltonian_style)
MPI.Finalize()

# ---------------------------------------------------------------------------
# Compute pairwise COHP
# ---------------------------------------------------------------------------
@info "Computing COHP for orbital pair ($i_orb, $j_orb) ..."
E_grid, cohp, energies, vecs = compute_cohp_pair(H, i_orb, j_orb)

# ---------------------------------------------------------------------------
# Write data
# ---------------------------------------------------------------------------

mkpath(output_path)

pair_tag  = "$(i_orb)_$(j_orb)"
data_file = joinpath(output_path, "cohp_$(pair_tag).txt")

open(data_file, "w") do io
    println(io, "# Orbital pair: i=$i_orb  j=$j_orb")
    println(io, "# Energy   COHP")
    writedlm(io, hcat(E_grid, cohp))
end
@info "COHP data written to $data_file"

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
plot_file = joinpath(output_path, "cohp_$(pair_tag).pdf")

p = plot(
    -cohp, E_grid;                          # conventional: –COHP on x-axis
    xlabel  = "–COHP (arb. units)",
    ylabel  = "Energy",
    title   = "COHP  (orbitals $i_orb – $j_orb)",
    legend  = false,
    lw      = 1.5,
    color   = :steelblue,
)
hline!(p, [0.0]; linestyle = :dot, color = :gray, lw = 1, label = "E = 0")
vline!(p, [0.0]; linestyle = :dash, color = :black, lw = 0.8, label = "COHP = 0")

savefig(p, plot_file)
@info "COHP plot saved to $plot_file"