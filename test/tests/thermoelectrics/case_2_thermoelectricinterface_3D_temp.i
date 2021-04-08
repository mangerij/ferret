[Mesh]
 [gen]
    ############################################
    ##
    ##  Type and dimension of the mesh
    ##
    ############################################

    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 10
    nz = 50
    xmin = -0.7
    xmax = 0.7
    ymin = -0.7
    ymax = 0.7
    zmin = -1.524
    zmax = 1.524
    elem_type = HEX8
  []

  [subdomains]
    type = SubdomainBoundingBoxGenerator
    input = gen
    bottom_left = '-0.7 -0.7 -1.524'
    block_id = 1
    top_right = '0.7 0.7 0.0'
    location = INSIDE
  [../]
  [film_interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = subdomains
    primary_block = 0
    paired_block = 1
    # new_boundary = 52
    new_boundary = 'primary0_interface'
  [../]
  [./break_boundary]
    input = film_interface
    type = BreakBoundaryOnSubdomainGenerator
  [../]
  []


# //3d
# [./subdomain1]
#     input = gen
#     type = SubdomainBoundingBoxGenerator
#     bottom_left = '-0.0007 -0.0007 -0.001524'
#     top_right = '0.0007 0.0007 0'
#     block_id = 1
#     location = INSIDE
#   [../]
#   [./interface]
#     type = SideSetsBetweenSubdomainsGenerator
#     input = subdomain1
#     primary_block = '0'
#     paired_block = '1'
#     new_boundary = 'primary0_interface'
#     # new_boundary = 52
#   [../]
#   [./break_boundary]
#     input = interface
#     type = BreakBoundaryOnSubdomainGenerator
#   [../]
# []


[GlobalParams]
  T = T
  potential_E_int = potential_E_int
[]

[Variables]
  [./potential_E_int]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-5
      max = 0.1e-5
    [../]
  [../]
  [./T]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-5
      max = 0.1e-5
    [../]
  [../]
[]

[AuxVariables]
  [./j_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./j_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./j_z]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./q_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./q_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./q_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./Electric_flux_x]
    type = Electric_flux
    variable = j_x
    T = T
    ecC = 'ecC'
    sbC = 'sbC'
    potential_E_int = 'potential_E_int'
    component = 0
  [../]
  [./Electric_flux_y]
    type = Electric_flux
    variable = j_y
    T = T
    ecC = 'ecC'
    sbC = 'sbC'
    potential_E_int = 'potential_E_int'
    component = 1
  [../]
  [./Electric_flux_z]
    type = Electric_flux
    variable = j_y
    T = T
    ecC = 'ecC'
    sbC = 'sbC'
    potential_E_int = 'potential_E_int'
    component = 2
  [../]

  [./Heat_flux_x]
    type = Heat_flux
    variable = q_x
    T = 'T'
    thC = 'thC'
    ecC = 'ecC'
    sbC = 'sbC'
    potential_E_int = 'potential_E_int'
    component = 0
  [../]
  [./heat_flux_y]
    type = Heat_flux
    variable = q_y
    T = 'T'
    thC = 'thC'
    ecC = 'ecC'
    sbC = 'sbC'
    potential_E_int = 'potential_E_int'
    component = 1
  [../]
  [./heat_flux_z]
    type = Heat_flux
    variable = q_y
    T = 'T'
    thC = 'thC'
    ecC = 'ecC'
    sbC = 'sbC'
    potential_E_int = 'potential_E_int'
    component = 2
  [../]
[]

[Kernels]
  ########BLOCK 1
  [./residualV_00]
    type = ResidualV
    variable = potential_E_int
    T = T
    component = 0
    block = 0
    ecC = ecC
    thC = thC
    sbC = sbC
  [../]
  [./residualT_00]
    type = ResidualT
    variable = T
    T = T
    component = 0
    block = 0
    ecC = ecC
    thC = thC
    sbC = sbC
  [../]
  [./residualV_01]
    type = ResidualV
    variable = potential_E_int
    T = T
    component = 1
    block = 0
    ecC = ecC
    thC = thC
    sbC = sbC
  [../]
  [./residualT_01]
    type = ResidualT
    variable = T
    T = T
    component = 1
    block = 0
    ecC = ecC
    thC = thC
    sbC = sbC
  [../]
  [./residualV_02]
    type = ResidualV
    variable = potential_E_int
    T = T
    component = 2
    block = 0
    ecC = ecC
    thC = thC
    sbC = sbC
  [../]
  [./residualT_02]
    type = ResidualT
    variable = T
    T = T
    component = 2
    block = 0
    ecC = ecC
    thC = thC
    sbC = sbC
  [../]

  ########BLOCK 1
  [./residualV_10]
    type = ResidualV
    variable = potential_E_int
    component = 0
    T = T
    block = 1
    ecC = ecC
    thC = thC
    sbC = sbC
  [../]
  [./residualT_10]
    type = ResidualT
    variable = T
    component = 0
    T = T
    ecC = ecC
    thC = thC
    sbC = sbC
    block = 1
  [../]
  [./residualV_11]
    type = ResidualV
    variable = potential_E_int
    component = 1
    T = T
    block = 1
    ecC = ecC
    thC = thC
    sbC = sbC
  [../]
  [./residualT_11]
    type = ResidualT
    variable = T
    component = 1
    T = T
    ecC = ecC
    thC = thC
    sbC = sbC
    block = 1
  [../]
  [./residualV_12]
    type = ResidualV
    variable = potential_E_int
    component = 2
    T = T
    block = 1
    ecC = ecC
    thC = thC
    sbC = sbC
  [../]
  [./residualT_12]
    type = ResidualT
    variable = T
    component = 2
    T = T
    ecC = ecC
    thC = thC
    sbC = sbC
    block = 1
  [../]
[]


[Materials]
  [./ThermoelectricProperties_block1]
    type = GenericConstantMaterial
    prop_names = 'ecC thC sbC'
    prop_values = '8.422e4 1.612 1.941e-4'
    block = 0
  [../]
  [./ThermoelectricProperties_block2]
    type = GenericConstantMaterial
    prop_names = 'ecC thC sbC'
    prop_values = '8.422e4 1.612 -1.941e-4'
    block = 1
  [../]
[]


[BCs]

  [./side_potential_front]
    type = DirichletBC
    variable = potential_E_int
    boundary = 'front'
    value = 0
  [../]
  [./side_potential_back]
    type = DirichletBC
    variable = potential_E_int
    boundary = 'back'
    value = 0.116
  [../]
  [./side_potential_top]
    type = NeumannBC
    variable = potential_E_int
    boundary = 'top'
    value = 0
  [../]
  [./side_potential_bottom]
    type = NeumannBC
    variable = potential_E_int
    boundary = 'bottom'
    value = 0
  [../]
  [./side_potential_left]
    type = NeumannBC
    variable = potential_E_int
    boundary = 'left'
    value = 0
  [../]
  [./side_potential_right]
    type = NeumannBC
    variable = potential_E_int
    boundary = 'right'
    value = 0
  [../]

  [./side_T_front]
    type = DirichletBC
    variable = T
    boundary = 'front'
    value = 280.0
  [../]
  [./side_T_back]
    type = DirichletBC
    variable = T
    boundary = 'back'
    value = 280.0
  [../]
  [./side_T_left]
    type = DirichletBC
    variable = T
    boundary = 'left'
    value = 280.0
  [../]
  [./top_T_right]
    type = DirichletBC
    variable = T
    boundary = 'right'
    value = 280.0
  [../]
  [./side_T_top]
    type = DirichletBC
    variable = T
    boundary = 'top'
    value = 280.0
  [../]
  [./top_T_bottom]
    type = DirichletBC
    variable = T
    boundary = 'bottom'
    value = 280.0
  [../]
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]

[Debug]
  show_var_residual_norms = true
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
    petsc_options_value = ' lu       NONZERO               1e-10'
  [../]
[]

[Postprocessors]

[]

[Outputs]
  exodus = true
  csv = true
  file_base = thermocouple_3D
[]
