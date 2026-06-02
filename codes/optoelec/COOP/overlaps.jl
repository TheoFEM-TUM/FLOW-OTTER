begin
    snormalize(vec) = norm(vec) ≈ 0 ? vec : normalize(vec)
    function get_axis(bond, orb)
        if occursin('s', orb)
            return bond
        else
            return ifelse(occursin("px", orb), [1, 0, 0], ifelse(occursin("py", orb), [0, 1, 0], [0, 0, 1]))
        end
    end

    function get_overlap_label(bond, is_onsite_R, i, j, orb_el1, orb_el2)
        el1, orb1 = split(orb_el1, "-"); el2, orb2 = split(orb_el2, "-")
        orb_label = ifelse(is_onsite_R && i == j, "onsite", "")

        if is_onsite_R && i == j
            if occursin('s', orb1) && occursin('s', orb2)
                return "$el1-$(el2)_s-onsite"
            elseif occursin('p', orb1) || occursin('p', orb2)
                return "$el1-$(el2)_p-onsite"
            end
        elseif norm(bond) ≈ 0 && is_onsite_R
            return "$el1-$(el2)_offdiag-onsite"
        else
            axis_i = get_axis(bond, orb1)
            axis_j = get_axis(bond, orb2)

            orbs = sort([orb1[1], orb2[1]])

            dot1 = abs(dot(snormalize(bond), snormalize(axis_i)))
            dot2 = abs(dot(snormalize(bond), snormalize(axis_j)))
            #println("")
            #@show dot1
            #@show dot2
            mi = dot1 > 0.6 ? "σ" : "π"
            mj = dot2 > 0.6 ? "σ" : "π"
            if mi == "π" && mj == "π" && el1 == "Br" && el2 == "Pb"
                #println("---")
                #@show norm(bond)
                #@show dot1
                #@show dot2
            end

            ms = sort([mi, mj])
            return replace("$el1-$(el2)_$(orbs[1])$(orbs[2])$(ms[1])$(ms[2])", "↑"=>"", "↓"=>"")
        end
    end

    function find_hoppings(inds, path, poscar_path)

        basis_orbs = Dict(
            "Pb"   => ["s↑", "s↓", "px↑", "px↓", "py↑", "py↓", "pz↑", "pz↓"],
            "Br"   => ["px↑", "px↓", "py↑", "py↓", "pz↑", "pz↓"],
            "Cs"   => ["s↑", "s↓"],
        )

        poscar = Hamster.read_poscar(poscar_path)

        atom_types = Hamster.number_to_element.(poscar.atom_types)
        atom_inds = collect(1:length(atom_types))
        filter!(i-> atom_types[i] ∈ keys(basis_orbs), atom_inds)

        atom_types = atom_types[atom_inds]
        rs_atom = poscar.rs_atom[:, atom_inds]
        
        basis = mapreduce(vcat, atom_types) do type
            type .* "-" .* basis_orbs[type]
        end



        overlaps = Dict{String, Vector{Float64}}()
        for (n, ind) in enumerate(inds)
            rs = h5read(joinpath(path, "structures.h5"), "configs")
            positions = zeros(3, length(basis)); k = 0
            for (i, type) in enumerate(atom_types), _ in eachindex(basis_orbs[type])
                k += 1
                positions[:, k] = rs_atom[:, i]
            end
            Hr, vecs = read_ham(ind, filename=joinpath(path, "ham.h5"), space="r")
            for R in eachindex(Hr)
                is_onsite_R = vecs[1, R] == 0.0 && vecs[2, R] == 0.0 && vecs[3, R] == 0.0
                @views for (i, j, Hij) in zip(findnz(Hr[R])...)
                    ri = Hamster.frac_to_cart(positions[:, i], poscar.lattice)
                    rj = Hamster.frac_to_cart(positions[:, j] .- vecs[:, R], poscar.lattice)
                    bond = rj .- ri

                    if norm(bond) > 7
                        for R_ in axes(vecs, 2)
                            trial_bond = Hamster.frac_to_cart(positions[:, j] .- vecs[:, R_] .- positions[:, i], poscar.lattice)
                            if norm(trial_bond) < norm(bond)
                                bond = trial_bond
                            end
                        end
                    end

                    overlap = get_overlap_label(bond, is_onsite_R, i, j, sort(["$(basis[i])", "$(basis[j])"])...)
                    #if overlap == "Br-Pb_ppππ" && real(Hij) < 0.0
                    #    println("---")
                    #    @show vecs[:, R]
                    #    @show norm(bond)
                    #    @show dot(snormalize(bond), [1, 0, 0])
                    #    @show dot(snormalize(bond), [0, 1, 0])
                    #    @show dot(snormalize(bond), [0, 0, 1])
                    #end
                    
                    if haskey(overlaps, overlap)
                        push!(overlaps[overlap], real(Hij))
                    else
                        overlaps[overlap] = Float64[real(Hij)]
                    end
                end
            end
            println(" $n / $(length(inds))")
        end

        return overlaps
    end

    Ts = [300]
    for (t, T) in enumerate(Ts)
        path = joinpath(@__DIR__, "04_Hamster_large_scale_lammps", "4x4x4", "$(T)K")
        poscar_path = joinpath(path, "POSCAR")
        ham_path = joinpath(path, "ham.h5")


        inds = h5read(joinpath(path, "hamster_out.h5"), "config_inds")
        overlaps = find_hoppings(inds, path, poscar_path)

        path_2 = joinpath(@__DIR__, "04_Hamster_large_scale_lammps", "4x4x4_noortho", "$(T)K")
        overlaps_tb = find_hoppings(inds, path_2, joinpath(path_2, "POSCAR"))

        for (overlap, values) in overlaps
            fig = plot(framestyle=:box, guidefontsize=labelsize, tickfontsize=ticksize, legendfontsize=labelsize)
            
            hist = fit(Histogram, values, nbins=100)
            xs = hist.edges[1]
            ys = hist.weights ./ sum(hist.weights)
            plot!(xs[2:end], ys, label="Batched optim", lw=2, color=colors[t])

            hist = fit(Histogram, overlaps_tb[overlap], nbins=100)
            xs_tb = hist.edges[1]
            ys_tb = hist.weights ./ sum(hist.weights)
            plot!(xs_tb[2:end], ys_tb, label="No batches", lw=2, color=colors[t+1])
            
            xlabel!("Energy (eV)")
            ylabel!("Occurence (-)")
            title!("$overlap")
            savefig("00_figures/cspbbr3_2x2x2_hr_histogram_$overlap.pdf")
        end
    end
end
