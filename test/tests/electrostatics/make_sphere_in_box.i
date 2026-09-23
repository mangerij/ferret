[Mesh]
  [sphere_fine]
    type = SphereMeshGenerator
    radius = 4
    nr = 3
    elem_type = TET4
  []

  [sphere]
    type = XYZDelaunayGenerator
    boundary = sphere_fine
    desired_volume = 2
    output_subdomain_name = sphere
    output_boundary = sphere_surface
  []

  [box]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 2
    xmin = -20
    xmax = 20
    ymin = -20
    ymax = 20
    zmin = -20
    zmax = 20
  []

  [combined]
    type = XYZDelaunayGenerator
    boundary = box
    holes = 'sphere'
    stitch_holes = true
    desired_volume = 250
    output_subdomain_name = medium
    output_boundary = box_surface
    hole_boundaries = 'sphere_iface'
  []

  [blocks]
    type = RenameBlockGenerator
    input = combined
    old_block = 'sphere medium'
    new_block = '1 2'
  []

  [s1]
    type = ParsedGenerateSideset
    input = blocks
    combinatorial_geometry = 'z > 19.999'
    included_boundaries = 'box_surface'
    new_sideset_name = 'zplus'
  []
  [s2]
    type = ParsedGenerateSideset
    input = s1
    combinatorial_geometry = 'z < -19.999'
    included_boundaries = 'box_surface'
    new_sideset_name = 'zminus'
  []
  [s3]
    type = ParsedGenerateSideset
    input = s2
    combinatorial_geometry = 'y < -19.999'
    included_boundaries = 'box_surface'
    new_sideset_name = 'yminus'
  []
  [s4]
    type = ParsedGenerateSideset
    input = s3
    combinatorial_geometry = 'x < -19.999'
    included_boundaries = 'box_surface'
    new_sideset_name = 'xminus'
  []
  [s5]
    type = ParsedGenerateSideset
    input = s4
    combinatorial_geometry = 'y > 19.999'
    included_boundaries = 'box_surface'
    new_sideset_name = 'yplus'
  []
  [s6]
    type = ParsedGenerateSideset
    input = s5
    combinatorial_geometry = 'x > 19.999'
    included_boundaries = 'box_surface'
    new_sideset_name = 'xplus'
  []

  [cleanup]
    type = BoundaryDeletionGenerator
    input = s6
    boundary_names = 'box_surface sphere_surface'
  []
  [ids]
    type = RenameBoundaryGenerator
    input = cleanup
    old_boundary = 'zplus zminus yminus xminus yplus xplus sphere_iface'
    new_boundary = '1 2 3 4 5 6 7'
  []
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]
