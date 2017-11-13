[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 12
  ny = 12
  nz = 6
  xmin = -5
  xmax = 5
  ymin = -5
  ymax = 5
  zmin = -2
  zmax = 2
  elem_type = HEX8
[]

[GlobalParams]
  len_scale = 1.0

  ##################################
  ##--Landau/coupling parameters--##
  ##################################

  alpha1 = -3.362e-1
  alpha11 = 7.155e-2
  alpha12 = 8.852e-2
  alpha111 = -4.448e-3
  alpha112 = 1.966e-3
  alpha123 = -5.324e-2
  alpha1111 = 1.863e-4
  alpha1112 = -4.705e-4
  alpha1122 = 9.552e-4
  alpha1123 = 3.076e-3

  G110 = 0.08
  G11_G110 = 0.6
  G12_G110 = 0
  G44_G110 = 0.3
  G44P_G110 = 0.3

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  potential_int = potential_int

[]



[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-7
      max = 0.1e-7
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-7
      max = 0.1e-7
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-7
      max = 0.1e-7
    [../]
  [../]

  [./potential_int]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]

  ### Operators for the polar field: ###
  [./bed_x]
    type = BulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = BulkEnergyDerivativeEighth
    variable = polar_z
    component = 2
  [../]


  [./walled_x]
    type = WallEnergyDerivative
    variable = polar_x
    component = 0
  [../]
  [./walled_y]
    type = WallEnergyDerivative
    variable = polar_y
    component = 1
  [../]
  [./walled_z]
    type = WallEnergyDerivative
    variable = polar_z
    component = 2
  [../]


  ###Time dependence
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
  [./polar_z_time]
     type = TimeDerivativeScaled
     variable = polar_z
    time_scale = 1.0
  [../]

  ####Electrostatics
  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
  [../]
  [./FE_E_int]
     type = Electrostatics
     permittivity = 0.08854187
     variable = potential_int
  [../]
  [./polar_electric_px]
     type = PolarElectricPStrong
     variable = polar_x
     component = 0
  [../]
  [./polar_electric_py]
     type = PolarElectricPStrong
     variable = polar_y
     component = 1
  [../]
  [./polar_electric_pz]
     type = PolarElectricPStrong
     variable = polar_z
     component = 2
  [../]

[]

[BCs]
  [./potential_cube5]
    type = DirichletBC
    boundary = 'front'
    value = 0.0
    variable = potential_int
  [../]
  [./potential_cube6]
    type = DirichletBC
    boundary = 'back'
    value = 0.00
    variable = potential_int
  [../]
  [./Periodic]
    [./xy]
      auto_direction = 'x y'
      variable = 'polar_x polar_y polar_z potential_int'
    [../]
  [../]
[]

[Postprocessors]
  [./Fbulk]
    type = BulkEnergyEighth
  [../]
  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = Fbulk
  [../]

[]


[UserObjects]
 [./kill]
  type = Terminator
  expression = 'perc_change <= 1.0e-4'
 [../]
[]

[Debug]
  show_var_residual_norms = false
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '     121              1e-10    1e-8         1e-6       bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    optimal_iterations = 4
    growth_factor = 1.4
    linear_iteration_ratio = 100
    cutback_factor =  0.55
  [../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 0.1
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_BFO_polar_Test
    elemental_as_nodal = true
  [../]
[]
