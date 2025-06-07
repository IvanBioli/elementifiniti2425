# Author: Ivan Bioli (https://github.com/IvanBioli)
begin
    import Pkg
    Pkg.activate("elementifinitiunipv_pkg")
    using Revise

    # Load the necessary files
    includet("../../modules/Meshing.jl")
    includet("../../modules/Quadrature_adv.jl")
    includet("../../modules/Assembly.jl")
    includet("../../modules/Assembly_mixed.jl")
end

h = 0.1
out_file = mesh_circle(h)
edges2nodes, elems2edges, elems2orientation = get_edges_info(out_file)
msh = Mesh(get_nodes_connectivity(out_file)...)
set_edges_info!(msh, edges2nodes, elems2edges, elems2orientation)