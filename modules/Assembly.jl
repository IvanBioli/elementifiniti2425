# Author: Ivan Bioli (https://github.com/IvanBioli)
# Inspired by code written by Jochen Hinz (https://github.com/JochenHinz) for MATH-451 @ EPFL

using Memoize
using SparseArrays

"""
    initialize_assembly!(mesh::Mesh)

Initialize the assembly process for the given mesh by computing the necessary geometric quantities.

# Arguments
- `mesh::Mesh`: The mesh object for which the assembly is initialized.
"""
function initialize_assembly!(mesh::Mesh)
    get_Bk!(mesh)
    get_detBk!(mesh)
    get_invBk!(mesh)
end

########################### GLOBAL ASSEMBLER ########################### 
"""
    assemble_global(mesh::Mesh, local_assembler!)

Assemble the global stiffness matrix and force vector for the given mesh using the provided local assembler function.

# Arguments
- `mesh::Mesh`: The mesh object.
- `local_assembler!`: A function that assembles the local stiffness matrix and force vector.

# Returns
- `K::SparseMatrixCSC`: The global stiffness matrix.
- `f::Vector`: The global force vector.
"""
function assemble_global(mesh::Mesh, local_assembler!)
    # Get number of dofs, triangles and basis functions
    n_basefuncs = 3
    n_dofs = get_ndofs(mesh)
    n_tri = get_ntri(mesh)
    # Allocate the element stiffness matrix and element force vector
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)
    # Allocate global force vector f
    f = zeros(n_dofs)
    # Allocate entries for assembling the global matrix as a sparse matrix
    I = zeros(9 * n_tri) # Row indices
    J = zeros(9 * n_tri) # Col indices
    K = zeros(9 * n_tri) # Entries
    loc = 1:9 # Current entries of I, J and A to be modified
    # Loop over all triangles
    for cell_index in 1:n_tri
        # Assemble the local matrices
        local_assembler!(Ke, fe, mesh, cell_index)
        # Get the local-to-global indices
        triangle = mesh.T[:, cell_index]
        # Add the local contribution to the global force vector
        f[triangle] += fe
        # Add the local contribution to the vectors of the assembly of the global stiffness matrix
        irows = repeat(triangle, 1, 3)
        icols = irows'
        I[loc] = reshape(irows, 9)
        J[loc] = reshape(icols, 9)
        K[loc] = reshape(Ke, 9)
        loc = loc .+ 9
    end
    # Assemble K as a sparse matrix
    K = sparse(I, J, K, n_dofs, n_dofs)
    return K, f
end

########################################################################
########################### LOCAL ASSEMBLERS ###########################
########################################################################

########################### POISSON PROBLEM ###########################
"""
    shapef_2DLFE(quadrule::TriQuad)

Compute the shape functions for the Poisson problem.

# Arguments
- `quadrule::TriQuad`: The quadrature rule.

# Returns
- `shapef`: The shape functions evaluated at the quadrature points.
"""
@memoize function shapef_2DLFE(quadrule::TriQuad)
    # Get quadrature points
    points = quadrule.points
    x, y = reshape(points[1, :], 1, :), reshape(points[2, :], 1, :)

    # points = [x; y], with shape (2, n)
    # The basis functions are:
    #   f1(x,y) = 1 - x - y     f2(x,y) = x     f3(x,y) = y
    # Hence we construct a (3, n) matrix shapef = [f1; f2; f3]
    shapef = [1 .- x .- y; x; y]
    return shapef
end

"""
    ∇shapef_2DLFE(quadrule::TriQuad)

Compute the gradients of the shape functions for the Poisson problem.

# Arguments
- `quadrule::TriQuad`: The quadrature rule.

# Returns
- `∇shapef`: The gradients of the shape functions evaluated at the quadrature points.
"""
@memoize function ∇shapef_2DLFE(quadrule::TriQuad)
    # Get quadrature points
    n = size(quadrule.points, 2)
    # points = [x; y], with shape (2, n)
    # The basis functions are:
    #   f1(x,y) = 1 - x - y     f2(x,y) = x          f3(x,y) = y
    # Hence the gradients are
    #   ∇f1(x,y) = [-1;-1]       ∇f2(x,y) = [1;0]     ∇f3(x,y) = [0;1]
    # Hence we construct a (2, 3, n) matrix obtained repeating n times
    # the matrix [-1 1 0;
    #             -1 0 1]
    ∇shapef = repeat([-1.0 1.0 0.0; -1.0 0.0 1.0], 1, 1, n)
    return ∇shapef
end

"""
    poisson_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f)

Assemble the local stiffness matrix and force vector for the Poisson problem.

# Arguments
- `Ke::Matrix`: The local stiffness matrix to be assembled.
- `fe::Vector`: The local force vector to be assembled.
- `mesh::Mesh`: The mesh object.
- `cell_index::Integer`: The index of the current cell.
- `f`: The source term function.

# Returns
- `Ke`: The assembled local stiffness matrix.
- `fe`: The assembled local force vector.
"""
function poisson_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f)
    n_basefuncs = 3
    # Reset to 0
    fill!(Ke, 0)
    fill!(fe, 0)
    # FIXME: It is sufficient to use Q0 quadrule to assemble the stiffness matrix exactly,
    # but here we show how to use a more general quadrature rule like Q2
    quadrule = Q2_ref
    points_e = mesh.Bk[:, :, cell_index] * quadrule.points .+ mesh.ak[:, cell_index]
    # Evaluate basis functions and their gradient
    shapef = shapef_2DLFE(quadrule)
    invBk = mesh.invBk[:, :, cell_index]
    ∇shapef = mapslices(x -> invBk' * x, ∇shapef_2DLFE(quadrule), dims=(1, 2))
    # Loop over quadrature points
    for (q_index, q_point) in enumerate(eachcol(points_e))
        # Get the quadrature weight
        dΩ = quadrule.weights[q_index] * mesh.detBk[cell_index]
        # Loop over test shape functions
        for i in 1:n_basefuncs
            v = shapef[i, q_index]
            ∇v = ∇shapef[:, i, q_index]
            # Add contribution to fe
            fe[i] += f(q_point) * v * dΩ
            # Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇u = ∇shapef[:, j, q_index]
                # Add contribution to Ke
                Ke[i, j] += (∇v ⋅ ∇u) * dΩ
            end
        end
    end
    return Ke, fe
end