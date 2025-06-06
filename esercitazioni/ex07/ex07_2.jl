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

################ PROBLEMA DEL LAYER CARATTERISTICO INTERNO ###############
begin
# Define data
ϵ = 1e-4
ϵ_fun(x) = ϵ
β = [2, 3]
β_fun(x) = β
f(x) = 1
g(x) = (3 * x[1] + x[2] - 1 > 0) ? 0 : 1 

h = 0.025
# Build the mesh with mesh-size h    
out_file = mesh_square(h)
T, p = get_nodes_connectivity(out_file)
msh = Mesh(T, p)
# Get Dirichlet dofs
bnd_tags, bnd_coords = get_boundary_nodes(out_file)
set_dirichletdofs!(msh, bnd_tags)
# Assemble the system
initialize_assembly!(msh)
end

##### STANDARD GALERKIN #####
begin
local_assembler(Ke, fe, msh, cell_index) = transport_assemble_local!(Ke, fe, msh, cell_index, f, ϵ_fun, β_fun)
A, b = assemble_global(msh, local_assembler)
# Condense system imposing Dirichlet BCs
A_cond, b_cond, uh = impose_dirichlet(A, b, g, msh)
# Solve linear system
uh[get_freedofs(msh)] = A_cond \ b_cond

# Plot the solution
plt = plot_flat(msh, uh; plot_msh = true)
savefig(plt, "./figures/ex07_2SG_solution_flat.pdf") # Save plot

plt = plot_surf(msh, uh; plot_msh = false)
PlotlyJS.savefig(plt, "./figures/ex07_2SG_solution_surf.pdf")
display(plt) # Save plot
end

##### NCAD STABILIZATION #####
begin
local_assembler(Ke, fe, msh, cell_index) = transport_assemble_local!(Ke, fe, msh, cell_index, f, ϵ_fun, β_fun; stab = "NCAD")
A, b = assemble_global(msh, local_assembler)
# Condense system imposing Dirichlet BCs
A_cond, b_cond, uh = impose_dirichlet(A, b, g, msh)
# Solve linear system
uh[get_freedofs(msh)] = A_cond \ b_cond

# Plot the solution
plt = plot_flat(msh, uh; plot_msh = true)
savefig(plt, "./figures/ex07_2NCAD_solution_flat.pdf") # Save plot

plt = plot_surf(msh, uh; plot_msh = false)
PlotlyJS.savefig(plt, "./figures/ex07_2NCAD_solution_surf.pdf")
display(plt) # Save plot
end

##### NCSD STABILIZATION #####
begin
local_assembler(Ke, fe, msh, cell_index) = transport_assemble_local!(Ke, fe, msh, cell_index, f, ϵ_fun, β_fun; stab = "NCSD")
A, b = assemble_global(msh, local_assembler)
# Condense system imposing Dirichlet BCs
A_cond, b_cond, uh = impose_dirichlet(A, b, g, msh)
# Solve linear system
uh[get_freedofs(msh)] = A_cond \ b_cond

# Plot the solution
plt = plot_flat(msh, uh; plot_msh = true)
savefig(plt, "./figures/ex07_2NCSD_solution_flat.pdf") # Save plot

plt = plot_surf(msh, uh; plot_msh = false)
PlotlyJS.savefig(plt, "./figures/ex07_2NCSD_solution_surf.pdf")
display(plt) # Save plot
end

##### SUPG STABILIZATION #####
begin
local_assembler(Ke, fe, msh, cell_index) = transport_assemble_local!(Ke, fe, msh, cell_index, f, ϵ_fun, β_fun; stab = "SUPG")
A, b = assemble_global(msh, local_assembler)
# Condense system imposing Dirichlet BCs
A_cond, b_cond, uh = impose_dirichlet(A, b, g, msh)
# Solve linear system
uh[get_freedofs(msh)] = A_cond \ b_cond

# Plot the solution
plt = plot_flat(msh, uh; plot_msh = true)
savefig(plt, "./figures/ex07_2SUPG_solution_flat.pdf") # Save plot

plt = plot_surf(msh, uh; plot_msh = false)
PlotlyJS.savefig(plt, "./figures/ex07_2SUPG_solution_surf.pdf")
display(plt) # Save plot
end