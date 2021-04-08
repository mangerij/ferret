[Mesh]
 [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 30
    ny = 60
    ymin = -0.0007
    ymax = 0.0007
    xmin = -0.001524
    xmax = 0.001524
    elem_type = QUAD4
  []

  [./subdomain1]
      input = gen
      type = SubdomainBoundingBoxGenerator
      bottom_left = '-0.001524 -0.0007 0'
      top_right = '0 0.0007 0'
      block_id = 1
      location = INSIDE
    [../]
    [./interface]
      type = SideSetsBetweenSubdomainsGenerator
      input = subdomain1
      primary_block = '0'
      paired_block = '1'
      new_boundary = 'primary0_interface'
    [../]
    [./break_boundary]
      input = interface
      type = BreakBoundaryOnSubdomainGenerator
    [../]

[]


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
      min = -0.1e-10
      max = 0.1e-10
    [../]
  [../]
  [./T]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0
      max = 0.1e-10
    [../]
  [../]
[]


[Kernels]
  ########BLOCK 0
  [./residualV_0]
    type = ResidualV
    variable = potential_E_int
    T = T
    component = 1
    block = 0
    ecC = ecC
    thC = thC
    sbC = sbC
  [../]
  [./residualT_0]
    type = ResidualT
    variable = T
    T = T
    component = 1
    block = 0
    ecC = ecC
    thC = thC
    sbC = sbC
  [../]


  ########BLOCK 1
  [./residualV_1]
    type = ResidualV
    variable = potential_E_int
    component = 0
    T = T
    block = 1
    ecC = ecC
    thC = thC
    sbC = sbC
  [../]
  [./residualT_1]
    type = ResidualT
    variable = T
    component = 0
    T = T
    ecC = ecC
    thC = thC
    sbC = sbC
    block = 1
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

  [./q_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./q_y]
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
[]

[Functions]
  [./k_func]
    type = ParsedFunction
    value = '1.758 - 5.290e-3 * t + 4.134e-5 * t^2'
  [../]
  [./lam_func]
    type = ParsedFunction
    value = '1.028e5 - 5.369e2 * t + 1.824 * t^2'
  [../]
  [./eps_func]
    type = ParsedFunction
    value = '1.802e-4 + 3.861e-7 * t - 9.582e-10 * t^2'
  [../]
  [./negative_eps_func]
    type = ParsedFunction
    value = '-1.802e-4 - 3.861e-7 * t + 9.582e-10 * t^2'
  [../]
[]


[Materials]

#Temperature dependence
  [./ThermoelectricProperties_0]
     type = ThermoelectricMaterial
     temp = 'T'
     thC_temperature_function = k_func
     ecC_temperature_function = lam_func
     sbC_temperature_function = eps_func
     block = 0
    [../]
  [./ThermoelectricProperties_1]
     type = ThermoelectricMaterial
     temp = 'T'
     thC_temperature_function = k_func
     ecC_temperature_function = lam_func
     sbC_temperature_function = negative_eps_func
     block = 1
    [../]
[]

[BCs]

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
    type = DirichletBC
    variable = potential_E_int
    boundary = 'left'
    value = 0.116
  [../]
  [./side_potential_right]
    type = DirichletBC
    variable = potential_E_int
    boundary = 'right'
    value = 0
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
  file_base = out_thermo_2D_temp
[]
