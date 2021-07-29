
  ##########################################
  ##
  ##   This model will simulate the Seebeck effect
  ##   of a simple thermocouple circuit.
  ##   
  ############################################

[Mesh]
    file = thermocouple_concept_v3.e
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
      min = 0.0
      max = 0.1e-10
    [../]
  [../]
[]

[Kernels]

  [./Seebeck_1]
    type = SeebeckEffect
    variable = potential_E_int
    potential_E_int = potential_E_int
    T = T
    sbC = sbC
    block = 1
    component = 0
  [../]

  [./q1_1_x]
    type = ThermalDiffusion
    variable = T
    T = T
    thC = thC
    component = 0
    block = 1
  [../]
  [./q1_1_y]
    type = ThermalDiffusion
    variable = T
    T = T
    thC = thC
    component = 1
    block = 1
  [../]
  [./q1_1_z]
    type = ThermalDiffusion
    variable = T
    T = T
    thC = thC
    component = 2
    block = 1
  [../]


  ########BLOCK 2
  [./Seebeck_2]
    type = SeebeckEffect
    variable = potential_E_int
    potential_E_int = potential_E_int
    T = T
    sbC = sbC
    block = 2
  [../]

  [./q1_2_x]
    type = ThermalDiffusion
    variable = T
    T = T
    thC = thC
    component = 0
    block = 2
  [../]
  [./q1_2_y]
    type = ThermalDiffusion
    variable = T
    T = T
    thC = thC
    component = 1
    block = 2
  [../]
  [./q1_2_z]
    type = ThermalDiffusion
    variable = T
    T = T
    thC = thC
    component = 2
    block = 2
  [../]
[]

[Functions]
  [./chromel_seebeck]
    type = PiecewiseLinear
    x = '100	200	300	425	500'
    y = '0.000011	0.000017	0.00002	0.00002	0.00002'
  [../]
  [./alumel_seebeck]
    type = PiecewiseLinear
    x = '100	200	300	425	500'
    y = '-0.000007	-0.000015	-0.00002	-0.000017	-0.000019'
  [../]
[]

[Materials]
  [./ThermoelectricProperties_chromel]
    type = ThermoelectricMaterial
    thC = 1.4
    electrical_conductivity = 8.422e4
    sbC_temperature_function = chromel_seebeck
    temp = T
    block = 1
  [../]
  [./ThermoelectricProperties_alumel]
    type = ThermoelectricMaterial
    thC = 6.5
    electrical_conductivity = 8.422e2
    sbC_temperature_function = alumel_seebeck
    temp = T
    block = 2
  [../]
[]


[BCs]
  [./measuring_junction]
    type = DirichletBC
    variable = T
    value = 300
    boundary = '2'
  [../]
  [./reference_junction]
    type = DirichletBC
    variable = T
    boundary = 1
    value = 273
  [../]

  [./measuring_junction_potential]
    type = DirichletBC
    variable = potential_E_int
    boundary = '2'
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

[Outputs]
  exodus = true
  csv = true
  file_base = out_thermocouple_concept
[]
