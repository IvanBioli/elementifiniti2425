# Author: Ivan Bioli (https://github.com/IvanBioli)
# Inspired by code written by Jochen Hinz (https://github.com/JochenHinz) for MATH-451 @ EPFL

using Memoize
using SparseArrays

"""
    initialize_assembly_mixed!(mesh::Mesh)

Initializes the assembly process for a mixed finite element method on the given `mesh`.

This function performs the following steps:
- Computes the element transformation matrices by calling `get_Bk!`.
- Calculates the determinants of the transformation matrices by calling `get_detBk!`.

# Arguments
- `mesh::Mesh`: The mesh data structure on which the assembly is to be initialized.

# Side Effects
Modifies the internal state of `mesh` by updating its transformation matrices and their determinants.

# See Also
- [`get_Bk!`](@ref)
- [`get_detBk!`](@ref)
"""
function initialize_assembly_mixed!(mesh::Mesh)
    get_Bk!(mesh)
    get_detBk!(mesh)
end

########################### GLOBAL ASSEMBLER ########################### 
"""
    assemble_global_mixed(mesh::Mesh, local_assembler!)

Assembles the global mixed finite element system matrices and right-hand side vector for a given mesh.

# Arguments
- `mesh::Mesh`: The mesh data structure containing information about elements, edges, and connectivity.
- `local_assembler!`: A function that computes the local element matrices and force vector for a given cell. It should have the signature `local_assembler!(Ae, Be, fe, mesh, cell_index)`.

# Returns
- `K`: The global system matrix assembled as a block matrix, where the upper-left block corresponds to the velocity matrix, the upper-right and lower-left blocks correspond to the coupling (pressure) matrices, and the lower-right block is a zero matrix.
- `b`: The global right-hand side vector, concatenating zeros for velocity DOFs and the assembled force vector for pressure DOFs.
"""
function assemble_global_mixed(mesh::Mesh, local_assembler!)
    ###########################################################################
    ############################ ADD YOUR CODE HERE ###########################
    ###########################################################################
end


"""
    shapef_2D_RT0FE(quadrule::TriQuad)

Compute the Raviart-Thomas RT0 vector-valued shape functions on the reference triangle at the given quadrature points.

# Arguments
- `quadrule::TriQuad`: A quadrature rule object for triangles, containing the quadrature points as a 2×n matrix.

# Returns
- `shapef::Array{Float64,3}`: A 3D array of size (2, 3, n), where:
    - The first dimension (2) corresponds to the vector components (x and y).
    - The second dimension (3) corresponds to the three RT0 basis functions.
    - The third dimension (n) corresponds to the number of quadrature points.

# Details
The RT0 basis functions on the reference triangle are:
- `f1(x, y) = [x; y - 1]`
- `f2(x, y) = [x; y]`
- `f3(x, y) = [x - 1; y]`

The function evaluates these basis functions at each quadrature point and returns them in a single array for efficient use in finite element assembly.

# Memoization
The function is memoized to cache results for repeated calls with the same quadrature rule.
"""
#@memoize 
function shapef_2D_RT0FE(quadrule::TriQuad)
    ###########################################################################
    ############################ ADD YOUR CODE HERE ###########################
    ###########################################################################
end


"""
    divshapef_2D_RT0FE(quadrule::TriQuad)

Compute the divergence of the Raviart-Thomas RT0 vector-valued basis functions on a 2D triangle at the given quadrature points.

# Arguments
- `quadrule::TriQuad`: A quadrature rule object containing the quadrature points for integration over a triangle.

# Returns
- `divshapef::Array{Int,3}`: A `(1, 3, n)` array, where `n` is the number of quadrature points. Each entry contains the divergence of the three RT0 basis functions, which is constant and equal to 2 for each function.

# Notes
- The RT0 basis functions on the reference triangle are:
    - `f₁(x, y) = [x; y - 1]`
    - `f₂(x, y) = [x; y]`
    - `f₃(x, y) = [x - 1; y]`
- The divergence of each basis function is constant and equal to 2.
- The result is repeated for each quadrature point.
"""
#@memoize 
function divshapef_2D_RT0FE(quadrule::TriQuad)
    ###########################################################################
    ############################ ADD YOUR CODE HERE ###########################
    ###########################################################################
end

"""
    darcy_assemble_local_mixed!(Ae::Matrix, Be::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f, μ)

Assembles the local matrices and vector for the mixed finite element formulation of the Darcy problem on a single cell.

# Arguments
- `Ae::Matrix`: Local stiffness matrix (to be filled in-place).
- `Be::Matrix`: Local divergence matrix (to be filled in-place).
- `fe::Vector`: Local right-hand side vector (to be filled in-place).
- `mesh::Mesh`: Mesh data structure containing geometry and connectivity.
- `cell_index::Integer`: Index of the current cell in the mesh.
- `f`: Function representing the source term, evaluated as `f(x)` at a point `x`.
- `μ`: Function representing the permeability coefficient, evaluated as `μ(x)` at a point `x`.

# Description
This function computes the local contributions to the global system for the mixed finite element discretization of the Darcy problem. It uses RT0 (Raviart-Thomas of lowest order) basis functions and Piola transformation for mapping reference shape functions to the physical element. The function performs numerical integration using a quadrature rule, and assembles the local matrices and vector by looping over quadrature points and basis functions.

# Notes
- The function assumes that the shape functions are oriented consistently using `elems2orientation`.
- The local matrices and vector are reset to zero at the beginning of the function.
- The function modifies `Ae`, `Be`, and `fe` in-place and returns them for convenience.

# Returns
- `(Ae, Be, fe)`: The assembled local stiffness matrix, divergence matrix, and right-hand side vector.
"""
########################### DARCY PROBLEM ###########################
function darcy_assemble_local_mixed!(Ae::Matrix, Be::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f, μ)
    ###########################################################################
    ############################ ADD YOUR CODE HERE ###########################
    ###########################################################################
end


########################## DEFINE FUNCTIONS TO COMPUTE ERROR ##########################
"""
    L2error_mixed_p(p::Function, ph::Vector, mesh::Mesh, ref_quad::TriQuad) -> Float64

Compute the L² error between an exact function `p` and its discrete approximation `ph` over a given mesh using quadrature.

# Arguments
- `p::Function`: The exact function to be evaluated.
- `ph::Vector`: Vector of discrete approximations (one per element).
- `mesh::Mesh`: The mesh structure containing element connectivity and geometry.
- `ref_quad::TriQuad`: Quadrature rule on the reference triangle, providing points and weights.

# Returns
- `Float64`: The computed L² error over the mesh.
"""
function L2error_mixed_p(p::Function, ph::Vector, mesh::Mesh, ref_quad::TriQuad)
    ###########################################################################
    ############################ ADD YOUR CODE HERE ###########################
    ###########################################################################
end

"""
    H1diverror_mixed_u(u::Function, divu::Function, uh::Vector, mesh::Mesh, ref_quad::TriQuad) -> Float64

Compute the H(div) error norm between the exact solution `(u, divu)` and the finite element solution `uh` 
for a mixed finite element method using Raviart-Thomas (RT0) elements on a triangular mesh.

# Arguments
- `u::Function`: The exact vector-valued solution function, accepting a point and returning a 2D vector.
- `divu::Function`: The exact divergence of the solution, accepting a point and returning a scalar.
- `uh::Vector`: The vector of degrees of freedom for the finite element solution (RT0 coefficients).
- `mesh::Mesh`: The mesh data structure, containing element connectivity, geometry, and orientation.
- `ref_quad::TriQuad`: Quadrature rule on the reference triangle, providing points and weights.

# Returns
- `Float64`: The combined H(div) error norm, i.e., `sqrt(∫|u - uh|^2) + sqrt(∫|divu - div(uh)|^2)` over the domain.
"""
function H1diverror_mixed_u(u::Function, divu::Function, uh::Vector, mesh::Mesh, ref_quad::TriQuad)
    ###########################################################################
    ############################ ADD YOUR CODE HERE ###########################
    ###########################################################################
end