[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 20
    ny = 20
    nz = 20
    xmin = 0.0
    xmax = 10.0
    ymin = 0.0
    ymax = 10.0
    zmin = 0.0
    zmax = 10.0
    elem_type = HEX8
  []
[]

[Functions]

  [./fn_2d]
    type = SDFourierNoise
    lambda = 2.5
    range = 0.2
    mid = 1.0
    seed = 1
  [../]
  [./fn_3d]
    type = S3DFourierNoise
    lambda = 2.5
    range = 0.2
    mid = 1.0
    seed = 1
  [../]
[]

[AuxVariables]
  [./noise_2d]
    order = FIRST
    family = LAGRANGE
  [../]
  [./noise_3d]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[ICs]
  [./ic_2d]
    type = FunctionIC
    variable = noise_2d
    function = fn_2d
  [../]
  [./ic_3d]
    type = FunctionIC
    variable = noise_3d
    function = fn_3d
  [../]
[]

[Postprocessors]
  [./avg_2d]
    type = ElementAverageValue
    variable = noise_2d
    execute_on = 'initial'
  [../]
  [./avg_3d]
    type = ElementAverageValue
    variable = noise_3d
    execute_on = 'initial'
  [../]
  [./max_2d]
    type = NodalExtremeValue
    variable = noise_2d
    value_type = max
    execute_on = 'initial'
  [../]
  [./min_2d]
    type = NodalExtremeValue
    variable = noise_2d
    value_type = min
    execute_on = 'initial'
  [../]
  [./max_3d]
    type = NodalExtremeValue
    variable = noise_3d
    value_type = max
    execute_on = 'initial'
  [../]
  [./min_3d]
    type = NodalExtremeValue
    variable = noise_3d
    value_type = min
    execute_on = 'initial'
  [../]
  [./p2d_z0]
    type = PointValue
    variable = noise_2d
    point = '3.0 4.0 0.0'
    execute_on = 'initial'
  [../]
  [./p2d_z5]
    type = PointValue
    variable = noise_2d
    point = '3.0 4.0 5.0'
    execute_on = 'initial'
  [../]
  [./dz_2d]
    type = DifferencePostprocessor
    value1 = p2d_z0
    value2 = p2d_z5
    execute_on = 'initial'
  [../]

  [./p3d_z0]
    type = PointValue
    variable = noise_3d
    point = '3.0 4.0 0.0'
    execute_on = 'initial'
  [../]
  [./p3d_z5]
    type = PointValue
    variable = noise_3d
    point = '3.0 4.0 5.0'
    execute_on = 'initial'
  [../]
  [./dz_3d]
    type = DifferencePostprocessor
    value1 = p3d_z0
    value2 = p3d_z5
    execute_on = 'initial'
  [../]
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  print_linear_residuals = false
  perf_graph = false
  [./out]
    type = Exodus
    file_base = fourier_noise_out
    execute_on = 'initial'
  [../]
[]
