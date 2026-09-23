[Mesh]
  [gen]

    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 10
    nz = 50
    xmin = -0.0007
    xmax = 0.0007
    ymin = -0.0007
    ymax = 0.0007
    zmin = -0.015
    zmax = 0.015
    elem_type = HEX8
  []
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

[Functions]

  [./sbC_function]
    type = PiecewiseLinear
    data_file = Typek_SeebeckCoefficient_input.csv
  [../]
[]

[Kernels]
  [./Seebeck_1]
    type = SeebeckEffect
    variable = potential_E_int
    potential_E_int = potential_E_int
    T = T
    sbC = sbC
  [../]
  [./q1_1_x]
    type = ThermalDiffusion
    variable = T
    T = T
    thC = thC
    component = 0
  [../]
  [./q1_1_y]
    type = ThermalDiffusion
    variable = T
    T = T
    thC = thC
    component = 1
  [../]
  [./q1_1_z]
    type = ThermalDiffusion
    variable = T
    T = T
    thC = thC
    component = 2
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
  [./zT]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./Electric_flux_x]
    type = ElectricFlux
    variable = j_x
    T = T
    ecC = ecC
    sbC = sbC
    potential_E_int = potential_E_int
    component = 0
  [../]
  [./Electric_flux_y]
    type = ElectricFlux
    variable = j_y
    T = T
    ecC = ecC
    sbC = sbC
    potential_E_int = potential_E_int
    component = 1
  [../]
  [./Electric_flux_z]
    type = ElectricFlux
    variable = j_y
    T = T
    ecC = ecC
    sbC = sbC
    potential_E_int = potential_E_int
    component = 2
  [../]
  [./Heat_flux_x]
    type = HeatFlux
    variable = q_x
    T = T
    thC = thC
    ecC = ecC
    sbC = sbC
    potential_E_int = 'potential_E_int'
    component = 0
  [../]
  [./heat_flux_y]
    type = HeatFlux
    variable = q_y
    T = T
    thC = thC
    ecC = ecC
    sbC = sbC
    potential_E_int = 'potential_E_int'
    component = 1
  [../]
  [./heat_flux_z]
    type = HeatFlux
    variable = q_y
    T = T
    thC = thC
    ecC = ecC
    sbC = sbC
    potential_E_int = 'potential_E_int'
    component = 2
  [../]
  [./zT_value]
    type = ZTAux
    variable = zT
    T = T
  [../]
[]

[Materials]

  [./ThermoelectricProperties_block1]
    type = ThermoelectricMaterial
    temp = T
    sbC = -0.000012
    thC = 2.5
    ecC = 10e4
    block = 0
  [../]

[]

[BCs]
  [./reference_junction_temperature]
    type = DirichletBC
    variable = T
    boundary = 'front'
    value = 273
  [../]
  [./measuring_junction_temperature]
    type = DirichletBC
    variable = T
    boundary = 'back'
    value = 298
  [../]

  [./measuring_junction_ground]
    type = DirichletBC
    variable = potential_E_int
    boundary = 'back'
    value = 0
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
    [./potential]
      type = PointValue
      point = '0 0 0.015'
      variable = potential_E_int
    [../]
    [./temperature]
      type = PointValue
      point = '0 0 -0.015'
      variable = T
    [../]
[]

[Outputs]
  exodus = true
  csv = true
  file_base = out_simpleTC_type_k
[]
