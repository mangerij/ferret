
[Mesh]
  file = exodus_disk_r8_h1.e
  uniform_refine = 1
[]

[MeshModifiers]
  [./centernodeset_1]
    type = AddExtraNodeset
    new_boundary = 'center_node_1'
    coord = '0.0 0.0 -0.5'
  [../]
  [./centernodeset_2]
    type = AddExtraNodeset
    new_boundary = 'center_node_2'
    coord = '0.0 0.0 0.5'
  [../]
  [./centernodeset_3]
    type = AddExtraNodeset
    new_boundary = 'center_node_3'
    coord = '0.0 0.0 0.166667'
  [../]
  [./centernodeset_4]
    type = AddExtraNodeset
    new_boundary = 'center_node_4'
    coord = '0.0 0.0 -0.166667'
  [../]
[]

[GlobalParams]
  len_scale = 1.0

  alpha1 = -0.09179 #room temp PTO
  alpha11 = 0.0706
  alpha111 = 0.0
  alpha12 = 0.1412
  alpha112 = 0.0
  alpha123 = 0.0

  G110 = 0.141
  G11/G110 = 0.0
  G12/G110 = 0
  G44/G110 = 1.0
  G44P/G110 = 1.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
[]

[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.2e-4
      max = 0.2e-4
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.2e-4
      max = 0.2e-4
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.2e-4
      max = 0.2e-4
    [../]
  [../]
[]


[Kernels]
  #Bulk energy density
  [./bed_x]
    type = BulkEnergyDerivativeSixth
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = BulkEnergyDerivativeSixth
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = BulkEnergyDerivativeSixth
    variable = polar_z
    component = 2
  [../]

  ##Wall energy penalty
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

  [./polar_xEstrong]
     type = PolarElectricEStrong
     variable = polar_x
  [../]
  [./polar_yEstrong]
     type = PolarElectricEStrong
     variable = polar_y
  [../]
  [./polar_zEstrong]
     type = PolarElectricEStrong
     variable = polar_z
  [../]

  [./anis_x]
    type = AnisotropyEnergy
    variable = polar_x
    component = 0
    K = 0.0565 #This sign seems to be right.
  [../]
  [./anis_y]
    type = AnisotropyEnergy
    variable = polar_y
    component = 1
    K = 0.0565
  [../]


  [./depol_z]
    type = DepolEnergy
    permitivitty = 0.008854187
    lambda = 0.0007 #sign I believe affects direction of the field
    variable = polar_z
    avePz = avePz
  [../]


    ##Time dependence
    [./polar_x_time]
       type = TimeDerivativeScaled
       variable=polar_x
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
[]


[BCs]
  [./center_pol_x]
    type = DirichletBC
    variable = 'polar_x'
    value = 0.0
    boundary = 'center_node_1 center_node_2 center_node_3 center_node_4'
  [../]
  [./center_pol_y]
    type = DirichletBC
    variable = 'polar_y'
    value = 0.0
    boundary = 'center_node_1 center_node_2 center_node_3 center_node_4'
  [../]
  [./side_neumann_x]
    variable = 'polar_x'
    type = NeumannBC
    value = 0.0
    boundary = '1'
  [../]
  [./side_neumann_y]
    variable = 'polar_y'
    type = NeumannBC
    value = 0.0
    boundary = '1'
  [../]
  [./side_neumann_z]
    variable = 'polar_z'
    type = NeumannBC
    value = 0.0
    boundary = '1'
  [../]
[]



[Postprocessors]
   [./avePz]
     type = ElementAverageValue
     variable = polar_z
     execute_on = 'initial linear nonlinear timestep_begin timestep_end'
   [../]
   [./Fbulk]
      type = BulkEnergy
      execute_on = 'timestep_end'
    [../]
    [./Fwall]
      type = WallEnergy
      execute_on = 'timestep_end'
    [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    250              1e-10      1e-8      1e-8      bjacobi   '
  [../]
[]

[Executioner]
  type = Transient
    [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    #iteration_window = 3
    optimal_iterations = 6 #should be 5 probably
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.8
[../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 0.7
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    execute_on = 'timestep_end'
    file_base = out_skyrm
    elemental_as_nodal = true
  [../]
[]
