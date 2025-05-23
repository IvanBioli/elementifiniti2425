# Author: Ivan Bioli (https://github.com/IvanBioli)

import Pkg
Pkg.activate("elementifinitiunipv_pkg")
using Revise

# Load the necessary files
includet("../../modules/Meshing.jl")
includet("../../modules/Quadrature_adv.jl")

import Meshes
using Plots
using LaTeXStrings

######################## u(x,y) = exp(x^2+y^2)  ########################
u(x) = exp(x[1]^2 + x[2]^2);
u_str = L"e^{x^2+y^2}"
I_exact = π * (exp(1) - 1)

begin
    mesh_sizes = 10 .^ (range(0, -2, 10))
    Q0_approx = similar(mesh_sizes)
    Q1_approx = similar(mesh_sizes)
    Q2_approx = similar(mesh_sizes)
    for (i, h) in enumerate(mesh_sizes)
        # Build the mesh with mesh-size h    
        out_file = mesh_circle(h)
        T, p = get_nodes_connectivity(out_file)
        mesh = Mesh(T, p)
        # Compute the integral
        Q0_approx[i] = Quadrature(u, mesh, Q0_ref)
        Q1_approx[i] = Quadrature(u, mesh, Q1_ref)
        Q2_approx[i] = Quadrature(u, mesh, Q2_ref)
    end
end

# Plotting
begin
    plt = plot(mesh_sizes, abs.(Q0_approx .- I_exact), xaxis=:log, yaxis=:log, label=L"$Q_0$", marker=:circle)
    plot!(mesh_sizes, abs.(Q1_approx .- I_exact), xaxis=:log, yaxis=:log, label=L"$Q_1$", marker=:rect)
    plot!(mesh_sizes, abs.(Q2_approx .- I_exact), xaxis=:log, yaxis=:log, label=L"$Q_2$", marker=:diamond)
    plot!(mesh_sizes, mesh_sizes .^ 2, xaxis=:log, yaxis=:log, linestyle=:dash, label=L"$\mathcal{O}(h^2)$")
    plot!(plt, legend=:bottomright)
    title!("Andamento dell'errore con $u_str sul cerchio")
    xlabel!(L"$h$")
    ylabel!("Errore")
    savefig(plt, "./figures/ex03_2_1.pdf") # Save plot
    plot(plt) # Show plot
end

######################## u(x,y) = 1  ########################
u(x) = 1.0;
u_str = L"1"
I_exact = π

begin
    mesh_sizes = 10 .^ (range(0, -2, 10))
    Q0_approx = similar(mesh_sizes)
    Q1_approx = similar(mesh_sizes)
    Q2_approx = similar(mesh_sizes)
    for (i, h) in enumerate(mesh_sizes)
        # Build the mesh with mesh-size h    
        out_file = mesh_circle(h)
        T, p = get_nodes_connectivity(out_file)
        mesh = Mesh(T, p)
        # Compute the integral
        Q0_approx[i] = Quadrature(u, mesh, Q0_ref)
        Q1_approx[i] = Quadrature(u, mesh, Q1_ref)
        Q2_approx[i] = Quadrature(u, mesh, Q2_ref)
    end
end

# Plotting
begin
    plt = plot(mesh_sizes, abs.(Q0_approx .- I_exact), xaxis=:log, yaxis=:log, label=L"$Q_0$", marker=:circle)
    plot!(mesh_sizes, abs.(Q1_approx .- I_exact), xaxis=:log, yaxis=:log, label=L"$Q_1$", marker=:rect)
    plot!(mesh_sizes, abs.(Q2_approx .- I_exact), xaxis=:log, yaxis=:log, label=L"$Q_2$", marker=:diamond)
    plot!(mesh_sizes, mesh_sizes .^ 2, xaxis=:log, yaxis=:log, linestyle=:dash, label=L"$\mathcal{O}(h^2)$")
    plot!(plt, legend=:bottomright)
    title!("Andamento dell'errore con $u_str sul cerchio")
    xlabel!(L"$h$")
    ylabel!("Errore")
    savefig(plt, "./figures/ex03_2_2.pdf") # Save plot
    plot(plt) # Show plot
end
