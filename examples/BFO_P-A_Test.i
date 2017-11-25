[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 16
  ny = 16
  nz = 10
  xmin = -4
  xmax = 4
  ymin = -4
  ymax = 4
  zmin = -3
  zmax = 3
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

  beta1 = -3.888e-1
  beta11 = 9.149e-2
  beta12 = 1.07e-2
  beta111 = -7.733e-3
  beta112 = -6.097e-3
  beta123 = -6.926
  beta1111 = 3.961e-4
  beta1112 = 1.294-e4
  beta1122 = 9.67e-4
  beta1123 = 8.115-e4

  t1111 = 1.165e-1
  t1122 = 1.539e-1
  t1212 = -1.925e-1
  t42111111 = 4.432e-3
  t24111111 = -2.662e-2
  t42111122 = -1.695e-2
  t24112222 = -2.157e-2
  t42112233 = -1.577e-2
  t24112233 = -1.133e-2
  t42112211 = 1.272e-2
  t24111122 = 1.214e-2
  t42111212 = -3.652e-2
  t42123312 = 4.972e-2
  t24121112 = -3.891e-3
  t24121233 = -2.554e-2
  t6211111111 = -1.327e-3
  t2611111111 = 3.663e-3
  t6211111122 = 7.066e-4
  t2611222222 = 1.7e-3
  t4411111111 = 6.133e-3
  t4411112222 = 2.438ee-3

  G110 = 0.15
  G11_G110 = 0.5
  G12_G110 = 0
  G44_G110 = 0.5 
  G44P_G110 = 0.5

  H110 = 0.015
  H11_H110 = 0.5
  H12_H110 = 0
  H44_H110 = 0.5
  H44P_H110 = 0.5

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  antiferrodis_A_x = antiferrodis_A_x
  antiferrodis_A_y = antiferrodis_A_y
  antiferrodis_A_z = antiferrodis_A_z

  #potential_int = potential_int
[]



[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-4
      max = 0.5e-4
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-4
      max = 0.5e-4
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-4
      max = 0.5e-4
    [../]
  [../]

  [./antiferrodis_A_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-6
      max = 0.5e-6
    [../]
  [../]
  [./antiferrodis_A_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-6
      max = 0.5e-6
    [../]
  [../]
  [./antiferrodis_A_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-6
      max = 0.5e-6
    [../]
  [../]

  #[./potential_int]
  #  order = FIRST
  #  family = LAGRANGE
  #[../]
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

  [./roto_polar_coupled_x]
    type = RotoPolarCoupledEnergyPolarDerivative
    variable = polar_x
    component = 0
  [../]
  [./roto_polar_coupled_y]
    type = RotoPolarCoupledEnergyPolarDerivative
    variable = polar_y
    component = 1
  [../]
  [./roto_polar_coupled_z]
    type = RotoPolarCoupledEnergyPolarDerivative
    variable = polar_z
    component = 2
  [../]

  #Operators for the AFD field

  [./rbed_x]
    type = RotoBulkEnergyDerivativeEighth
    variable = antiferrodis_A_x
    component = 0
  [../]
  [./rbed_y]
    type = RotoBulkEnergyDerivativeEighth
    variable = antiferrodis_A_y
    component = 1
  [../]
  [./rbed_z]
    type = RotoBulkEnergyDerivativeEighth
    variable = antiferrodis_A_z
    component = 2
  [../]

  [./roto_walled_x]
    type = AFDAntiphaseEnergyDerivative
    variable = antiferrodis_A_x
    component = 0
  [../]
  [./roto_walled_y]
    type = AFDAntiphaseEnergyDerivative
    variable = antiferrodis_A_y
    component = 1
  [../]
  [./roto_walled_z]
    type = AFDAntiphaseEnergyDerivative
    variable = antiferrodis_A_z
    component = 2
  [../]

  [./roto_dis_coupled_x]
    type = RotoPolarCoupledEnergyDistortDerivative
    variable = antiferrodis_A_x
    component = 0
  [../]
  [./roto_dis_coupled_y]
    type = RotoPolarCoupledEnergyDistortDerivative
    variable = antiferrodis_A_y
    component = 1
  [../]
  [./roto_dis_coupled_z]
    type = RotoPolarCoupledEnergyDistortDerivative
    variable = antiferrodis_A_z
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

  [./antiferrodis_A_x_time]
    type = TimeDerivativeScaled
    variable = antiferrodis_A_x
    time_scale = 1.0
  [../]
  [./antiferrodis_A_y_time]
     type = TimeDerivativeScaled
     variable = antiferrodis_A_y
    time_scale = 1.0
  [../]
  [./antiferrodis_A_z_time]
     type = TimeDerivativeScaled
     variable = antiferrodis_A_z
    time_scale = 1.0
  [../]

  ####Electrostatics
  #[./polar_x_electric_E]
  #   type = PolarElectricEStrong
  #   variable = potential_int
  #[../]
  #[./FE_E_int]
  #   type = Electrostatics
  #   permittivity = 0.08854187
  #   variable = potential_int
  #[../]
  #[./polar_electric_px]
  #   type = PolarElectricPStrong
  #   variable = polar_x
  #   component = 0
  #[../]
  #[./polar_electric_py]
  #   type = PolarElectricPStrong
  #   variable = polar_y
  #   component = 1
  #[../]
  #[./polar_electric_pz]
  #   type = PolarElectricPStrong
  #   variable = polar_z
  #   component = 2
  #[../]

[]

[BCs]
  #[./potential_cube5]
  #  type = DirichletBC
  #  boundary = 'front'
  #  value = 0.0
  #  variable = potential_int
  #[../]
  #[./potential_cube6]
  #  type = DirichletBC
  #  boundary = 'back'
  #  value = 0.00
  #  variable = potential_int
  #[../]
  [./Periodic]
    [./xy]
      auto_direction = 'x y z'
      variable = 'polar_x polar_y polar_z antiferrodis_A_x antiferrodis_A_y antiferrodis_A_z'
    [../]
  [../]
[]

[Postprocessors]
  [./FbulkP]
    type = BulkEnergyEighth
    execute_on = 'initial timestep_end'
  [../]
  [./FbulkA]
    type = RotoBulkEnergyEighth
    execute_on = 'initial timestep_end'
  [../]
  [./FwallP]
    type = WallEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./FwallA]
    type = AFDWallEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./FcoupledPA]
    type = RotoPolarCoupledEnergyEighth
    execute_on = 'initial timestep_end'
  [../]
  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = FbulkP
    execute_on = 'initial timestep_end'
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
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type  '
    petsc_options_value = '     121              1e-10    1e-8         1e-5      bjacobi   '
  [../]
[]

[Executioner]
  type = Transient
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.08
    optimal_iterations = 4
    growth_factor = 1.4
    linear_iteration_ratio = 100
    cutback_factor =  0.55
  [../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 0.55
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_BFO_polar_dist_Test
    elemental_as_nodal = true
  [../]
[]
