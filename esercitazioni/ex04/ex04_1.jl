# Author: Ivan Bioli (https://github.com/IvanBioli)
begin
    import Pkg
    Pkg.activate("elementifinitiunipv_pkg")
    using Revise

    # Load the necessary files
    includet("../../modules/Meshing.jl")
    includet("../../modules/Quadrature_adv.jl")
    includet("../../modules/Assembly.jl")
end

################### CODE FOR SANITY CHECKS ###################
# Common data
begin
    h = 0.1
    out_file = mesh_circle(h)
    T, p = get_nodes_connectivity(out_file)
    msh = Mesh(T, p)
    f(x) = x[1]
end

# Check that your code runs and assembles the matrix
begin
    initialize_assembly!(msh)
    local_assembler(Ke, fe, msh, cell_index) = poisson_assemble_local!(Ke, fe, msh, cell_index, f)
    A, b = assemble_global(msh, local_assembler)
end

# Asseble with Gridap
begin
    # To check the assembly of the matrix with Gridap, use the following code
    using GridapGmsh
    using Gridap
    model = GridapGmsh.GmshDiscreteModel(out_file)
    order = 1
    reffe = ReferenceFE(lagrangian, Float64, order)
    V0 = TestFESpace(model, reffe; conformity=:H1)
    Ug = TrialFESpace(V0)
    degree = 2
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)
    a_(u, v) = ∫(∇(v) ⋅ ∇(u)) * dΩ
    b_(v) = ∫(v * f) * dΩ
    op = AffineFEOperator(a_, b_, Ug, V0)
    A_gridap = get_matrix(op)  # Assembles the stiffness matrix
    b_gridap = get_vector(op)  # Assembles the load vector
end

# Check that the errors are small
begin
    println("-------------------------------------------------------------")
    println("Relative error on A:\t $(norm(A - A_gridap) / norm(A_gridap))")
    println("Relative error on b:\t $(norm(b - b_gridap) / norm(b_gridap))")
end