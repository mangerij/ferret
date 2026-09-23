[Mesh]
  [gen]

    type = GeneratedMeshGenerator
    dim = 2

    nx = 150
    ny = 20

    xmin = -30.0
    xmax = 30.0
    ymin = -4.0
    ymax = 4.0

    elem_type = QUAD4
  []

  [subdomains1]

    type = SubdomainBoundingBoxGenerator
    input = gen
    bottom_left = '-30.0 3.0 0.0'
    block_id = 1
    top_right = '30.0 4.0 0.0'
  []
  [subdomains2]
    type = SubdomainBoundingBoxGenerator
    input = subdomains1
    bottom_left = '-30.0 -4.0 0.0'
    block_id = 2
    top_right = '30.0 -3.0 0.0'
  []
[]

[GlobalParams]

  len_scale = 1.0

  polar_x = polar_x
  polar_y = polar_y

  potential_E_int = potential_E_int
[]

[Variables]

  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block = '0 1 2'
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-6
      max = 0.01e-6
      seed = 6
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block = '0 1 2'
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-6
      max = 0.01e-6
      seed = 6
    [../]
  [../]

  [./potential_E_int]
    order = FIRST
    family = LAGRANGE
    block = '0 1 2'
  [../]
[]

[Materials]

  [./Landau_P_FE]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '-0.027722 -0.6381 0.0 7.89 0.0 0.0 0.0 0.0 0.0 0.0'
    block = '0'
  [../]

  [./Landau_P_die_layer]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '8.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    block = '1 2'
  [../]

  [./Landau_G]
    type = GenericConstantMaterial
    prop_names = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '0.5 0.0138 0.0 0.0138 0.0'
  [../]

  [./in_P_susc]
    type = GenericConstantMaterial
    prop_names = 'chi'
    prop_values = '10.0'
  [../]

  [./permitivitty]
    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '1.84167'
    block = '1 2'
  [../]

  [./eps1]

    type = GenericConstantMaterial
    prop_names = 'eps1 eps2 eps3'
    prop_values = '3.692195979 1.5450556315 0.0'
    block = '0'
  [../]
[]

[Kernels]

  [./bed_FE_x]
    type = BulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
    block = '0'
  [../]

  [./bed_FE_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
    block = '0'
  [../]

  [./ip_FE_x]
    type = InPlaneSusceptibilityDerivative
    variable = polar_x
    block = '0'
  [../]

  [./bed_die_x]
    type = BulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
    block = '1 2'
  [../]

  [./bed_die_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
    block = '1 2'
  [../]

  [./walled_x]
    type = WallEnergyDerivative
    variable = polar_x
    component = 0
    block = '0'
  [../]
  [./walled_y]
    type = WallEnergyDerivative
    variable = polar_y
    component = 1
    block = '0'
  [../]

  [./polar_x_electric_E]
    type = PolarElectricEStrong
    variable = potential_E_int
    block = '0 1 2'
  [../]

  [./FE_E_int]
    type = AnisotropicElectrostatics
    variable = potential_E_int
    block = '0'
  [../]

  [./die_E_int]
    type = Electrostatics
    variable = potential_E_int
    block = '1 2'
  [../]

  [./polar_electric_px]
    type = PolarElectricPStrong
    variable = polar_x
    component = 0
    block = '0 1 2'
  [../]
  [./polar_electric_py]
    type = PolarElectricPStrong
    variable = polar_y
    component = 1
    block = '0 1 2'
  [../]

  [./polar_x_time]
    type = TimeDerivativeScaled
    variable = polar_x

    time_scale = 1.0
  [../]
  [./polar_y_time]
    type = TimeDerivativeScaled
    variable=polar_y
    time_scale = 1.0
  [../]
[]

[BCs]

  [./Periodic]
    [./xyz]
      auto_direction = 'x'
      variable = 'polar_x polar_y'
    [../]
  [../]

  [./boundary_top_grounding]
    type = DirichletBC
    boundary = 'top'
    variable = potential_E_int
    value = 0.0
  [../]

  [./boundary_bottom_grounding]
    type = DirichletBC
    boundary = 'bottom'
    variable = potential_E_int
    value = 0.0
  [../]
[]

[Postprocessors]

  [./Fbulk]
    type = BulkEnergyEighth
    execute_on = 'initial timestep_end'
    block = '0'
  [../]
  [./Fwall]
    type = WallEnergy
    execute_on = 'initial timestep_end'
    block = '0'
  [../]
  [./Felec]
    type = ElectrostaticEnergy
    execute_on = 'initial timestep_end'
    block = '0'
  [../]
  [./Ftotal]
    type = LinearCombinationPostprocessor
    pp_names = 'Fbulk Fwall Felec'
    pp_coefs = ' 1 1 1'
    execute_on = 'initial timestep_end'
  [../]
[]

[Preconditioning]

  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type  -build_twosided'
    petsc_options_value = '    160            1e-10      1e-8      1e-6       bjacobi       allreduce'
  [../]
[]

[Executioner]

  type = Transient
  solve_type = 'NEWTON'
  scheme = 'implicit-euler'
  dtmin = 1e-13

  dtmax = 8.0

  l_max_its = 200

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 12
    cutback_factor = 0.75
    linear_iteration_ratio = 1000
    dt = 0.1
  [../]
  verbose = true
  num_steps = 2
[]

[Outputs]

  print_linear_residuals = false
  perf_graph = false

  [./out]
    type = Exodus
    file_base = out_FE_8nmFilm_Die_4nmLayer_U0
    elemental_as_nodal = true
    time_step_interval = 1
  [../]
  [./outCSV]
    type = CSV
    file_base = out_FE_8nmFilm_Die_4nmLayer_U0
  [../]
[]
