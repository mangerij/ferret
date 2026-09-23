[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 8
    ny = 8
    nz = 8
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
    zmin = 0.0
    zmax = 1.0
    elem_type = HEX8
  []
[]

[AuxVariables]
  [./azimuth_phi1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.3
      max = 0.35
      seed = 2
    [../]
  [../]
  [./polar_theta1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.0001
      max = 0.0002
      seed = 37
    [../]
  [../]

  [./azimuth_phi2]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.3
      max = 0.35
      seed = 2
    [../]
  [../]
  [./polar_theta2]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 3.1415
      max = 3.1416
      seed = 37
    [../]
  [../]
  [./mag1_norm]
    order = FIRST
    family = LAGRANGE
  [../]
  [./mag2_norm]
    order = FIRST
    family = LAGRANGE
  [../]
  [./neel_z]
    order = FIRST
    family = LAGRANGE
  [../]
  [./mtot_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Variables]
  [./mag1_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi1
      theta = polar_theta1
      M0s = 1.0
      component = 0
    [../]
  [../]
  [./mag1_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi1
      theta = polar_theta1
      M0s = 1.0
      component = 1
    [../]
  [../]
  [./mag1_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi1
      theta = polar_theta1
      M0s = 1.0
      component = 2
    [../]
  [../]
  [./mag2_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi2
      theta = polar_theta2
      M0s = 1.0
      component = 0
    [../]
  [../]
  [./mag2_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi2
      theta = polar_theta2
      M0s = 1.0
      component = 1
    [../]
  [../]
  [./mag2_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi2
      theta = polar_theta2
      M0s = 1.0
      component = 2
    [../]
  [../]
[]

[AuxKernels]
  [./k_mag1_norm]
    type = ParsedAux
    variable = mag1_norm
    coupled_variables = 'mag1_x mag1_y mag1_z'
    expression = 'sqrt(mag1_x*mag1_x + mag1_y*mag1_y + mag1_z*mag1_z)'
    execute_on = 'initial'
  [../]
  [./k_mag2_norm]
    type = ParsedAux
    variable = mag2_norm
    coupled_variables = 'mag2_x mag2_y mag2_z'
    expression = 'sqrt(mag2_x*mag2_x + mag2_y*mag2_y + mag2_z*mag2_z)'
    execute_on = 'initial'
  [../]

  [./k_neel_z]
    type = ParsedAux
    variable = neel_z
    coupled_variables = 'mag1_z mag2_z'
    expression = '0.5*(mag1_z - mag2_z)'
    execute_on = 'initial'
  [../]
  [./k_mtot_z]
    type = ParsedAux
    variable = mtot_z
    coupled_variables = 'mag1_z mag2_z'
    expression = '0.5*(mag1_z + mag2_z)'
    execute_on = 'initial'
  [../]
[]

[Postprocessors]
  [./m1_norm_min]
    type = NodalExtremeValue
    variable = mag1_norm
    value_type = min
    execute_on = 'initial'
  [../]
  [./m1_norm_max]
    type = NodalExtremeValue
    variable = mag1_norm
    value_type = max
    execute_on = 'initial'
  [../]
  [./m2_norm_min]
    type = NodalExtremeValue
    variable = mag2_norm
    value_type = min
    execute_on = 'initial'
  [../]
  [./m2_norm_max]
    type = NodalExtremeValue
    variable = mag2_norm
    value_type = max
    execute_on = 'initial'
  [../]
  [./neel_z_avg]
    type = ElementAverageValue
    variable = neel_z
    execute_on = 'initial'
  [../]
  [./mtot_z_avg]
    type = ElementAverageValue
    variable = mtot_z
    execute_on = 'initial'
  [../]
  [./m1x_min]
    type = NodalExtremeValue
    variable = mag1_x
    value_type = min
    execute_on = 'initial'
  [../]
  [./m1x_max]
    type = NodalExtremeValue
    variable = mag1_x
    value_type = max
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
    file_base = random_constrained_vector_out
    execute_on = 'initial'
  [../]
[]
