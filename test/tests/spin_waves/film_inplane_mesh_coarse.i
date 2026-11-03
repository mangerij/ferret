[Mesh]
  [box]
    type = CartesianMeshGenerator
    dim = 2
    dx = '40.0'
    dy = '80.0 40.0 20.0 10.0 5.0 2.5 1.25 0.5 6.0 0.5 1.25 2.5 5.0 10.0 20.0 40.0 80.0'
    ix = '20'
    iy = '1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1'
    subdomain_id = '0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0'
  []
  [center]
    type = TransformGenerator
    input = box
    transform = TRANSLATE
    vector_value = '0.0 -162.25 0.0'
  []
  [block_names]
    type = RenameBlockGenerator
    input = center
    old_block = '0 1'
    new_block = 'vacuum film'
  []
  [film_left]
    type = ParsedGenerateSideset
    input = block_names
    combinatorial_geometry = 'x < 1e-6'
    included_subdomains = film
    new_sideset_name = film_left
  []
  [film_right]
    type = ParsedGenerateSideset
    input = film_left
    combinatorial_geometry = 'x > 40.0 - 1e-6'
    included_subdomains = film
    new_sideset_name = film_right
  []
  [ends_secondary]
    type = LowerDBlockFromSidesetGenerator
    input = film_right
    sidesets = left
    new_block_id = 11
    new_block_name = ends_secondary
  []
  [ends_primary]
    type = LowerDBlockFromSidesetGenerator
    input = ends_secondary
    sidesets = right
    new_block_id = 12
    new_block_name = ends_primary
  []
  [film_secondary]
    type = LowerDBlockFromSidesetGenerator
    input = ends_primary
    sidesets = film_left
    new_block_id = 13
    new_block_name = film_secondary
  []
  [film_primary]
    type = LowerDBlockFromSidesetGenerator
    input = film_secondary
    sidesets = film_right
    new_block_id = 14
    new_block_name = film_primary
  []
  [phi_pin]
    type = ExtraNodesetGenerator
    input = film_primary
    new_boundary = phi_pin
    coord = '20.0 162.25 0.0'
  []
[]
