# Author: Ivan Bioli (https://github.com/IvanBioli)
# Inspired by code written by Jochen Hinz (https://github.com/JochenHinz) for MATH-451 @ EPFL

"""
    struct TriQuad

A structure representing a triangular quadrature rule.

# Fields
- `name::String`: The name of the quadrature rule.
- `order::Integer`: The order of the quadrature rule.
- `points::Matrix`: The quadrature points.
- `weights::Array`: The quadrature weights.
"""
struct TriQuad
    name::String
    order::Integer
    points::Matrix
    weights::Array
end

Q0_ref = TriQuad(
    "Q0",
    2,
    reshape([
            1 / 3;
            1 / 3
        ], (2, 1)),
    [1 / 2]
)

Q1_ref = TriQuad(
    "Q1",
    2,
    [
        0.0 1.0 0.0;
        0.0 0.0 1.0
    ],
    [
        1/6 1/6 1/6
    ]
)

Q2_ref = TriQuad(
    "Q2",
    3,
    [
        0.5 0.5 0.0;
        0.0 0.5 0.5
    ],
    [
        1/6 1/6 1/6
    ]
)

"""
    Quadrature(u, mesh::Mesh, ref_quad::TriQuad)

Perform numerical integration of a function over a mesh using a given quadrature rule.

# Arguments
- `u`: The function to be integrated.
- `mesh::Mesh`: The mesh object.
- `ref_quad::TriQuad`: The reference quadrature rule.

# Returns
- `I_approx::Float64`: The approximate integral of the function over the mesh.
"""
function Quadrature(u, mesh::Mesh, ref_quad::TriQuad)
    # Compute matrices for pushforward of reference element
    Bk, ak = get_Bk!(mesh)
    # Compute the absolute value of the determinant
    detBk = get_detBk!(mesh)
    # Get quadrature points and weights on the reference element
    points_refelem, weights_refelem = ref_quad.points, ref_quad.weights
    points_elem = zeros(Float64, size(points_refelem))
    u_evals = zeros(Float64, size(weights_refelem))

    # Loop across all elements
    n_tri = size(mesh.T, 2)
    I_approx::Float64 = 0
    for i = 1:n_tri
        points_elem = Bk[:, :, i] * points_refelem .+ ak[:, i] # Points in the current element
        u_evals = eval_u(u, points_elem, mesh, i, ref_quad)
        I_approx += sum(u_evals .* weights_refelem) * detBk[i]
    end
    return I_approx
end

# Evaluation of a function
"""
    eval_u(u::Function, points_elem::Matrix, mesh::Mesh, tri_idx::Integer, quadrule::TriQuad)

Evaluate a function at given points within an element.

# Arguments
- `u::Function`: The function to be evaluated.
- `points_elem::Matrix`: The points at which to evaluate the function.
- `mesh::Mesh`: The mesh object (ignored).
- `tri_idx::Integer`: The index of the current element (ignored).
- `quadrule::TriQuad`: The quadrature rule (ignored).

# Returns
- `u_evals::Matrix`: The evaluated function values at the given points.
"""
function eval_u(u::Function, points_elem::Matrix, mesh::Mesh, tri_idx::Integer, quadrule::TriQuad)
    return mapslices(u, points_elem, dims=1)
end

"""
    eval_u(uh::Vector, points_elem::Matrix, mesh::Mesh, tri_idx::Integer, quadrule::TriQuad)

Evaluate a linear finite element solution at given quadrature points within an element.

# Arguments
- `uh::Vector`: The finite element solution vector.
- `points_elem::Matrix`: The points at which to evaluate the solution (ignored).
- `mesh::Mesh`: The mesh object.
- `tri_idx::Integer`: The index of the current element.
- `quadrule::TriQuad`: The quadrature rule.

# Returns
- `uh_evals::Matrix`: The evaluated solution values at the given points.
"""
function eval_u(uh::Vector, points_elem::Matrix, mesh::Mesh, tri_idx::Integer, quadrule::TriQuad)
    ######################
    ### COMPLETARE QUI ###
    ######################
end

"""
    L2error(u::Function, uh::Vector, mesh::Mesh, ref_quad::TriQuad)

Compute the L2 error between a function and a finite element solution over a mesh.

# Arguments
- `u::Function`: The exact solution function.
- `uh::Vector`: The finite element solution vector.
- `mesh::Mesh`: The mesh object.
- `ref_quad::TriQuad`: The reference quadrature rule.

# Returns
- `L2_error::Float64`: The L2 error between the exact solution and the finite element solution.
"""
function L2error(u::Function, uh::Vector, mesh::Mesh, ref_quad::TriQuad)
    ######################
    ### COMPLETARE QUI ###
    ######################
end

"""
    H1semierror(∇u::Function, uh::Vector, mesh::Mesh, ref_quad::TriQuad)

Compute the H1 semi-norm error between the gradient of a function and a finite element solution over a mesh.

# Arguments
- `∇u::Function`: The gradient of the exact solution function.
- `uh::Vector`: The finite element solution vector.
- `mesh::Mesh`: The mesh object.
- `ref_quad::TriQuad`: The reference quadrature rule.

# Returns
- `H1_semi_error::Float64`: The H1 semi-norm error between the gradient of the exact solution and the finite element solution.
"""
function H1semierror(∇u::Function, uh::Vector, mesh::Mesh, ref_quad::TriQuad)
     ######################
    ### COMPLETARE QUI ###
    ######################
end