# Author: Ivan Bioli (https://github.com/IvanBioli)
begin
    import Pkg
    Pkg.activate("elementifinitiunipv_pkg")
    using Revise

    # Load the necessary files
    includet("../../modules/Meshing.jl")
    includet("../../modules/Quadrature_adv.jl")
    includet("../../modules/Assembly.jl")

    import Meshes
    using Plots
    using LaTeXStrings
    import GridapGmsh
end

###################### PROBLEMA CON SOLUZIONE LISCIA ######################
begin
# Define data
β = [1, 0]
β_fun(x) = β
f(x, ϵ) = 2*ϵ*π^2  * sin(π  * x[1]) * sin(π * x[2]) + β[1] * π * cos(π  * x[1]) * sin(π * x[2])+ β[2] * π * sin(π  * x[1]) * cos(π * x[2])
u_exact(x) = sin(π  * x[1]) * sin(π * x[2])
∇u_exact(x) = [π * cos(π  * x[1]) * sin(π * x[2]), π * sin(π  * x[1]) * cos(π * x[2])]
g = u_exact
end;

ϵ_vec = [1e-4, 1e-8, 1e-12]
for ϵ = ϵ_vec
    ϵ_fun(x) = ϵ
    fϵ(x) = f(x, ϵ)

    msh_sizes = 2 .^ range(-1, -6, 6)
    SG_L2err = similar(msh_sizes)
    SG_H1err = similar(msh_sizes)
    NCAD_L2err = similar(msh_sizes)
    NCAD_H1err = similar(msh_sizes)
    NCSD_L2err = similar(msh_sizes)
    NCSD_H1err = similar(msh_sizes)
    SUPG_L2err = similar(msh_sizes)
    SUPG_H1err = similar(msh_sizes)
    for (i, h) in enumerate(msh_sizes)
        # Build the mesh with mesh-size h    
        out_file = mesh_square(h)
        T, p = get_nodes_connectivity(out_file)
        msh = Mesh(T, p)
        # Get Dirichlet dofs
        bnd_tags, bnd_coords = get_boundary_nodes(out_file)
        set_dirichletdofs!(msh, bnd_tags)
        initialize_assembly!(msh)
        
        # Standard Galerkin
        local_assembler_sg(Ke, fe, msh, cell_index) = transport_assemble_local!(Ke, fe, msh, cell_index, fϵ, ϵ_fun, β_fun)
        A, b = assemble_global(msh, local_assembler_sg)
        A_cond, b_cond, uh_SG = impose_dirichlet(A, b, g, msh)
        uh_SG[get_freedofs(msh)] = A_cond \ b_cond

        # NCAD
        local_assembler_ncad(Ke, fe, msh, cell_index) = transport_assemble_local!(Ke, fe, msh, cell_index, fϵ, ϵ_fun, β_fun, stab="NCAD")
        A, b = assemble_global(msh, local_assembler_ncad)
        A_cond, b_cond, uh_NCAD = impose_dirichlet(A, b, g, msh)
        uh_NCAD[get_freedofs(msh)] = A_cond \ b_cond

        # NCSD
        local_assembler_ncsd(Ke, fe, msh, cell_index) = transport_assemble_local!(Ke, fe, msh, cell_index, fϵ, ϵ_fun, β_fun, stab="NCSD")
        A, b = assemble_global(msh, local_assembler_ncsd)
        A_cond, b_cond, uh_NCSD = impose_dirichlet(A, b, g, msh)
        uh_NCSD[get_freedofs(msh)] = A_cond \ b_cond

        # SUPG
        local_assembler_supg(Ke, fe, msh, cell_index) = transport_assemble_local!(Ke, fe, msh, cell_index, fϵ, ϵ_fun, β_fun, stab="SUPG")
        A, b = assemble_global(msh, local_assembler_supg)
        A_cond, b_cond, uh_SUPG = impose_dirichlet(A, b, g, msh)
        uh_SUPG[get_freedofs(msh)] = A_cond \ b_cond

        # Compute errors
        Ih_uexact = dropdims(mapslices(u_exact, msh.p; dims=1); dims=1) # Interpolation of the exact solution
        SG_L2err[i] = L2error(u_exact, uh_SG, msh, Q2_ref)
        SG_H1err[i] = sqrt(H1semierror(∇u_exact, uh_SG, msh, Q2_ref)^2 + SG_L2err[i]^2)
        NCAD_L2err[i] = L2error(u_exact, uh_NCAD, msh, Q2_ref)
        NCAD_H1err[i] = sqrt(H1semierror(∇u_exact, uh_NCAD, msh, Q2_ref)^2 + NCAD_L2err[i]^2)
        NCSD_L2err[i] = L2error(u_exact, uh_NCSD, msh, Q2_ref)
        NCSD_H1err[i] = sqrt(H1semierror(∇u_exact, uh_NCSD, msh, Q2_ref)^2 + NCSD_L2err[i]^2)
        SUPG_L2err[i] = L2error(u_exact, uh_SUPG, msh, Q2_ref)
        SUPG_H1err[i] = sqrt(H1semierror(∇u_exact, uh_SUPG, msh, Q2_ref)^2 + SUPG_L2err[i]^2)
    end

    plt = plot(msh_sizes, SG_L2err, xaxis=:log, yaxis=:log, label="SG", marker=:diamond)
    plot!(msh_sizes, NCAD_L2err, xaxis=:log, yaxis=:log, label="NCAD", marker=:rect)
    plot!(msh_sizes, NCSD_L2err, xaxis=:log, yaxis=:log, label="NCSD", marker=:circ)
    plot!(msh_sizes, SUPG_L2err, xaxis=:log, yaxis=:log, label="SUPG", marker=:ltriangle)
    plot!(msh_sizes, NCAD_L2err[1] .* (msh_sizes / msh_sizes[1]) / 2, xaxis=:log, yaxis=:log, linestyle=:dash, label=L"$\mathcal{O}(h)$")
    plot!(msh_sizes, SUPG_L2err[1] .* ((msh_sizes/msh_sizes[1]) .^ 2) / 2, xaxis=:log, yaxis=:log, linestyle=:dash, label=L"$\mathcal{O}(h^2)$")

    title!(L"$\varepsilon = %$ϵ $")
    xlabel!(L"$h$")
    ylabel!(L"${\|\|u-u_h \|\|}_{L^2(\Omega)}$")
    plot!(legend=:bottomright)
    savefig(plt, "./figures/ex07_3_eps$ϵ.pdf") # Save plot
    plot(plt) # Show plot
end