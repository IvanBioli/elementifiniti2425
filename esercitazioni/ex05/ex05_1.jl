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
end

######################## HOMOGENEOUS DIRICHLET BCs  ########################
begin
# Define data
f(x) = 1.0
u_exact(x) = 0.25 - 0.25 * (x[1]^2 + x[2]^2)
∇u_exact(x) = [-0.5 * x[1], -0.5 * x[2]]
g = u_exact
end

begin
msh_sizes = 2 .^ range(-1, -6, 6)
Linf_error = similar(msh_sizes)
L2Q0_error = similar(msh_sizes)
L2Q2_error = similar(msh_sizes)
H1Q0_error = similar(msh_sizes)
H1Q2_error = similar(msh_sizes)
for (i, h) in enumerate(msh_sizes)
    # Build the mesh with mesh-size h    
    out_file = mesh_circle(h)
    T, p = get_nodes_connectivity(out_file)
    msh = Mesh(T, p)
    # Get Dirichlet dofs
    bnd_tags, bnd_coords = get_boundary_nodes(out_file)
    set_dirichletdofs!(msh, bnd_tags)

    # Assemble the system
    initialize_assembly!(msh)
    local_assembler(Ke, fe, msh, cell_index) = poisson_assemble_local!(Ke, fe, msh, cell_index, f)
    A, b = assemble_global(msh, local_assembler)
    # Condense system imposing Dirichlet BCs
    A_cond, b_cond, uh = impose_dirichlet(A, b, g, msh)
    # Solve linear system
    uh[get_freedofs(msh)] = A_cond \ b_cond

    # Compute Linfity error
    Ih_uexact = dropdims(mapslices(u_exact, msh.p; dims=1); dims=1) # Interpolation of the exact solution
    Linf_error[i] = norm(Ih_uexact - uh, Inf)
    L2Q0_error[i] = L2error(u_exact, uh, msh, Q0_ref)
    L2Q2_error[i] = L2error(u_exact, uh, msh, Q2_ref)
    H1Q0_error[i] = sqrt(H1semierror(∇u_exact, uh, msh, Q0_ref)^2 + L2Q0_error[i]^2)
    H1Q2_error[i] = sqrt(H1semierror(∇u_exact, uh, msh, Q2_ref)^2 + L2Q2_error[i]^2)
end
end

# Plotting
begin
    plt = plot(msh_sizes, Linf_error, xaxis=:log, yaxis=:log, label=L"${\|\|u-u_h \|\|}_{L^\infty}$", marker=:circle)
    plot!(msh_sizes, L2Q0_error, xaxis=:log, yaxis=:log, label=L"${\|\|u-u_h \|\|}_{L^2}$", marker=:rect)
    plot!(msh_sizes, H1Q0_error, xaxis=:log, yaxis=:log, label=L"${\|\|u-u_h \|\|}_{H^1}$", marker=:utriangle)
    plot!(msh_sizes, Linf_error[1] .* msh_sizes .^ 2, xaxis=:log, yaxis=:log, linestyle=:dash, label=L"$\mathcal{O}(h^2)$")

    title!("Formula dei baricentri")
    xlabel!(L"$h$")
    ylabel!("Errore")
    plot!(legend=:bottomright)
    savefig(plt, "./figures/ex05_1_1.pdf") # Save plot
    plot(plt) # Show plot
end

begin
    plt = plot(msh_sizes, Linf_error, xaxis=:log, yaxis=:log, label=L"${\|\|u-u_h \|\|}_{L^\infty}$", marker=:circle)
    plot!(msh_sizes, L2Q2_error, xaxis=:log, yaxis=:log, label=L"${\|\|u-u_h \|\|}_{L^2}$", marker=:rect)
    plot!(msh_sizes, H1Q2_error, xaxis=:log, yaxis=:log, label=L"${\|\|u-u_h \|\|}_{H^1}$", marker=:utriangle)
    plot!(msh_sizes, Linf_error[1] .* msh_sizes .^ 2, xaxis=:log, yaxis=:log, linestyle=:dash, label=L"$\mathcal{O}(h^2)$")
    plot!(msh_sizes, H1Q2_error[1] .* msh_sizes, xaxis=:log, yaxis=:log, linestyle=:dash, label=L"$\mathcal{O}(h)$")

    title!("Formula dei punti medi")
    xlabel!(L"$h$")
    ylabel!("Errore")
    plot!(legend=:bottomright)
    savefig(plt, "./figures/ex05_1_2.pdf") # Save plot
    plot(plt) # Show plot
end

################### PLOT THE SOLUTION ###################
begin
    # Get the mesh
    h = 0.1
    out_file = mesh_circle(h)
    T, p = get_nodes_connectivity(out_file)
    msh = Mesh(T, p)
    bnd_tags, bnd_coords = get_boundary_nodes(out_file)
    set_dirichletdofs!(msh, bnd_tags)

    # Define the data
    f(x) = 1.0
    u_exact(x) = 0.25 - 0.25 * (x[1]^2 + x[2]^2)
    g = u_exact

    # Assemble
    initialize_assembly!(msh)
    local_assembler(Ke, fe, msh, cell_index) = poisson_assemble_local!(Ke, fe, msh, cell_index, f)
    A, b = assemble_global(msh, local_assembler)
    A_cond, b_cond, uh = impose_dirichlet(A, b, g, msh)
    uh[get_freedofs(msh)] = A_cond \ b_cond

    # Plot the solution
    plt = plot_flat(msh, uh; plot_msh = true)
    savefig(plt, "./figures/ex05_1_3_solution_flat.pdf") # Save plot
    plt = plot_surf(msh, uh; plot_msh = true)
    PlotlyJS.savefig(plt, "./figures/ex05_1_3_solution_surf.pdf") # Save plot
end