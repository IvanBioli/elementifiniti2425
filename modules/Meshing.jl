# Author: Ivan Bioli (https://github.com/IvanBioli)

using Gridap
import Meshes, CairoMakie
import Gmsh: gmsh
using LinearAlgebra
import TriplotRecipes
import PlotlyJS

import Gmsh: gmsh

"""
    mesh_square(h::Number; display::Bool=false, out_file::String="./tmp_square.msh", v1 = [0, 0], v2 = [1, 1])

Generate a square mesh using Gmsh.

# Arguments
- `h::Number`: The mesh size.
- `display::Bool`: Whether to display the mesh using Gmsh's GUI.
- `out_file::String`: The output file path for the mesh.
- `v1`: The coordinates of the first vertex of the square.
- `v2`: The coordinates of the opposite vertex of the square.

# Returns
- `out_file::String`: The output file path for the mesh.
"""
function mesh_square(h::Number; display::Bool=false, out_file::String="./tmp_square.msh", v1 = [0, 0], v2 = [1, 1])
    # Initialize Gmsh
    gmsh.initialize()

    # Define points of the square
    p1 = gmsh.model.geo.addPoint(v1[1], v1[2], 0, h, 1)
    p2 = gmsh.model.geo.addPoint(v2[1], v1[2], 0, h, 2)
    p3 = gmsh.model.geo.addPoint(v2[1], v2[2], 0, h, 3)
    p4 = gmsh.model.geo.addPoint(v1[1], v2[2], 0, h, 4)

    # Define lines for each side of the square
    l1 = gmsh.model.geo.addLine(p1, p2, 1)
    l2 = gmsh.model.geo.addLine(p2, p3, 2)
    l3 = gmsh.model.geo.addLine(p3, p4, 3)
    l4 = gmsh.model.geo.addLine(p4, p1, 4)

    # Create a loop from these lines
    loop1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4], 1)

    # Define the surface within this loop
    surface = gmsh.model.geo.addPlaneSurface([loop1], 1)

    # Synchronize CAD entities
    gmsh.model.geo.synchronize()

    # Label the boundaries with physical groups
    boundary = gmsh.model.addPhysicalGroup(1, [l1, l2, l3, l4])
    gmsh.model.setPhysicalName(1, boundary, "boundary")
    boundary_d0 = gmsh.model.addPhysicalGroup(0, [p1, p2, p3, p4])
    gmsh.model.setPhysicalName(0, boundary_d0, "boundary")

    # Lower edge (bottom)
    lower = gmsh.model.addPhysicalGroup(1, [l1])
    gmsh.model.setPhysicalName(1, lower, "lower")
    lower_d0 = gmsh.model.addPhysicalGroup(0, [p1, p2])
    gmsh.model.setPhysicalName(0, lower_d0, "lower")

    # Right edge
    right = gmsh.model.addPhysicalGroup(1, [l2])
    gmsh.model.setPhysicalName(1, right, "right")
    right_d0 = gmsh.model.addPhysicalGroup(0, [p2, p3])
    gmsh.model.setPhysicalName(0, right_d0, "right")
    # Upper edge (top)
    upper = gmsh.model.addPhysicalGroup(1, [l3])
    gmsh.model.setPhysicalName(1, upper, "upper")
    upper_d0 = gmsh.model.addPhysicalGroup(0, [p3, p4])
    gmsh.model.setPhysicalName(0, upper_d0, "upper")
    # Left edge
    left = gmsh.model.addPhysicalGroup(1, [l4])
    gmsh.model.setPhysicalName(1, left, "left")
    left_d0 = gmsh.model.addPhysicalGroup(0, [p4, p1])
    gmsh.model.setPhysicalName(0, left_d0, "left")

    gmsh.model.addPhysicalGroup(1, [l2])  # right edge
    gmsh.model.setPhysicalName(1, 2, "right_edge")
    gmsh.model.addPhysicalGroup(1, [l3])  # top edge
    gmsh.model.setPhysicalName(1, 3, "upper_edge")
    gmsh.model.addPhysicalGroup(1, [l4])  # left edge
    gmsh.model.setPhysicalName(1, 4, "left_edge")

    # Define physical group for the surface
    surface_group = gmsh.model.addPhysicalGroup(2, [surface])
    gmsh.model.setPhysicalName(2, surface_group, "domain")

    # Generate the 2D mesh and save it
    gmsh.model.mesh.generate(2)
    gmsh.write(out_file)

    # Display mesh if requested
    if display
        gmsh.fltk.run()
    end

    # Close the Gmsh API
    gmsh.finalize()
    return out_file
end

"""
    mesh_circle(h::Number; display::Bool=false, out_file::String="./tmp_circle.msh")

Generate a circular mesh using Gmsh.

# Arguments
- `h::Number`: The mesh size.
- `display::Bool`: Whether to display the mesh using Gmsh's GUI.
- `out_file::String`: The output file path for the mesh.

# Returns
- `out_file::String`: The output file path for the mesh.
"""
function mesh_circle(h::Number; display::Bool=false, out_file::String="./tmp_circle.msh")
    gmsh.initialize()
    # Add vertices of the square (x, y, z, mesh_size_close_to_point, tag)
    p1 = gmsh.model.geo.addPoint(1, 0, 0, h)
    p2 = gmsh.model.geo.addPoint(0, 1, 0, h)
    p3 = gmsh.model.geo.addPoint(-1, 0, 0, h)
    p4 = gmsh.model.geo.addPoint(0, -1, 0, h, 4)
    c = gmsh.model.geo.addPoint(0, 0, 0, h)
    # Add lines (start, center, end, tag)
    arc1 = gmsh.model.geo.addCircleArc(p1, c, p2)
    arc2 = gmsh.model.geo.addCircleArc(p2, c, p3)
    arc3 = gmsh.model.geo.addCircleArc(p3, c, p4)
    arc4 = gmsh.model.geo.addCircleArc(p4, c, p1)
    # Add loop
    loop = gmsh.model.geo.addCurveLoop([arc1, arc2, arc3, arc4])
    surface = gmsh.model.geo.addPlaneSurface([loop])

    # Synchronize CAD entities
    gmsh.model.geo.synchronize()

    # Label the boundaries with physical groups
    boundary = gmsh.model.addPhysicalGroup(1, [arc1, arc2, arc3, arc4])
    gmsh.model.setPhysicalName(1, boundary, "boundary")
    boundary_d0 = gmsh.model.addPhysicalGroup(0, [p1, p2, p3, p4])
    gmsh.model.setPhysicalName(0, boundary_d0, "boundary")
    # Define physical group for the surface
    surface_group = gmsh.model.addPhysicalGroup(2, [surface])
    gmsh.model.setPhysicalName(2, surface_group, "domain")

    # Perform meshing and save 
    gmsh.model.mesh.generate(2)
    gmsh.write(out_file)

    # Display if asked
    if display
        gmsh.fltk.run()
    end

    # Close Gmsh.jl API
    gmsh.finalize()
    return out_file
end

"""
    get_nodes_connectivity(out_file::String)

Retrieve the nodes and the connectivity matrix from the mesh file.

# Arguments
- `out_file::String`: The mesh file path.

# Returns
- `elem_node_tags::Matrix{Int64}`: The element connectivity matrix.
- `node_coords::Matrix{Float64}`: The node coordinates matrix.
"""
function get_nodes_connectivity(out_file)
    # Load the mesh
    gmsh.initialize()
    gmsh.open(out_file)
    # Synchronize to make sure everything is loaded
    gmsh.model.geo.synchronize()

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_tags = convert(Array{Int64}, node_tags)
    nnodes = length(node_tags)
    node_coords = reshape(node_coords, 3, nnodes)
    node_coords = node_coords[1:2, :]

    # Retrieve elements and their connectivity
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2)
    @assert(length(elem_types) == 1)
    elem_tags = convert(Array{Int64}, elem_tags[1])
    nelems = length(elem_tags)
    elem_node_tags = reshape(convert(Array{Int64}, elem_node_tags[1]), 3, nelems)

    gmsh.finalize()

    return elem_node_tags, node_coords
end

"""
    get_physical_groups()

Retrieve the physical groups and their names from the mesh.

# Returns
- `physical_groups::Dict{String,Tuple{Int,Int}}`: A dictionary of physical group names and their dimensions and tags.
"""
function get_physical_groups()
    # Retrieve physical groups and their names
    physical_groups = Dict{String,Tuple{Int,Int}}()
    for dim in 1:2  # 1D for boundaries, 2D for surfaces
        phys_groups = gmsh.model.getPhysicalGroups(dim)
        for (dim, tag) in phys_groups
            phys_name = gmsh.model.getPhysicalName(dim, tag)
            physical_groups[phys_name] = (dim, tag)
        end
    end
    return physical_groups
end

"""
    get_boundary_nodes(out_file::String; labels::Vector{String}=["boundary"])

Retrieve the boundary nodes from the mesh file.

# Arguments
- `out_file::String`: The mesh file path.
- `labels::Vector{String}`: The labels of the physical groups to retrieve.

# Returns
- `all_node_tags::Vector{Int64}`: The tags of the boundary nodes.
- `all_node_coords::Matrix{Float64}`: The coordinates of the boundary nodes.
"""
function get_boundary_nodes(out_file; labels = ["boundary"])
    # Load the mesh
    gmsh.initialize()
    gmsh.open(out_file)
    # Synchronize to make sure everything is loaded
    gmsh.model.geo.synchronize()

    # Initialize arrays to store node tags and coordinates
    all_node_tags = Int64[]
    all_node_coords = Array{Float64}(undef, 2, 0)  # Empty array with 0 rows and 2 columns
    for l in labels
        phys_groups = get_physical_groups()
        node_tags, node_coords = gmsh.model.mesh.getNodesForPhysicalGroup(phys_groups[l]...)
        node_tags = convert(Array{Int64}, node_tags)
        nnodes = length(node_tags)
        node_coords = reshape(node_coords, 3, nnodes)
        node_coords = node_coords[1:2, :]
        all_node_tags = vcat(all_node_tags, node_tags)
        all_node_coords = hcat(all_node_coords, node_coords)
    end
    idx = unique(z -> all_node_tags[z], 1:length(all_node_tags))
    gmsh.finalize()

    return all_node_tags[idx], all_node_coords[:,idx]
end

"""
    to_Meshes(T_mat::Matrix{Integer}, p_mat::Matrix{Real})

Convert the connectivity and point matrices to a Meshes.SimpleMesh object.

# Arguments
- `T_mat::Matrix{Integer}`: The connectivity matrix.
- `p_mat::Matrix{Real}`: The point coordinates matrix.

# Returns
- `mesh::Meshes.SimpleMesh`: The converted mesh object.
"""
function to_Meshes(T_mat, p_mat)
    # Convert into vector of tuples of points 
    points = Meshes.Point.([Tuple{Real,Real}(p_mat[:, i]) for i in 1:size(p_mat, 2)])
    # Connect the points using the provided connectivity matrix T 
    connec = Meshes.connect.([Tuple{Integer,Integer,Integer}(T_mat[:, i]) for i in 1:size(T_mat, 2)])
    # Mesh 
    mesh = Meshes.SimpleMesh(points, connec)
    return mesh
end