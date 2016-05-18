
[Mesh]
  file = embedded_single_sphere_3.e
  #uniform_refine = 1
[]

[GlobalParams]
  len_scale = 1.0
  alpha1 = -0.1722883 # (3.766(T-765.1)*10^5) C^{-2} nm^2 (T = 293 K)
  alpha11 = -0.07253
  alpha111 = 0.26
  alpha12 = 0.75
  alpha112 = 0.61
  alpha123 = -3.67
  G110 = 0.173
  G11/G110 = 0.6
  G12/G110 = 0
  G44/G110 = 0.3
  G44P/G110 = 0.3
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  displacements = 'disp_x disp_y disp_z'
[]



[Variables]
  [./polar_x]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      block = '1'
      type = FluctuationsIC
      epsilon = 1.0e-5
      q1 = '678 469 744'
      q2 = '328 942 895'
      q3 = '980 1121 343'
      q4 = '1234 970 593'
      h = 2.48
    [../]
  [../]
  [./polar_y]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      block = '1'
      type = FluctuationsIC
      epsilon = 1.0e-5
      q1 = '1011 1184 1329'
      q2 = '512 586 1106'
      q3 = '547 1297 962'
      q4 = '511 833 371'
      h = 1.35
    [../]
  [../]
  [./polar_z]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      block = '1'
      type = FluctuationsIC
      epsilon = 1.0e-5
      q1 = '602 666 338'
      q2 = '678 469 744'
      q3 = '1107 529 1175'
      q4 = '1133 1108 532'
      h = 1.08
    [../]
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
    [./InitialCondition]
      type = FluctuationsIC
      block = '1 2'
      epsilon = 1.0e-5
      q1 = '354 1005 645'
      q2 = '715 1065 1132'
      q3 = '391 305 1106'
      q4 = '1053 1116 627'
      h = 0.25
    [../]
  [../]
  [./disp_x]
    block = '1 2'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      block = '1 2'
      type = FluctuationsIC
      epsilon = 1.0e-5
      q1 = '256 1246 730'
      q2 = '1234 770 1257'
      q3 = '1262 279 467'
      q4 = '467 944 1326'
      h = -0.0017
    [../]
  [../]
  [./disp_y]
    block = '1 2'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FluctuationsIC
      block = '1 2'
      epsilon = 1.0e-5
      q1 = '1522 415 1009'
      q2 = '147 622 351'
      q3 = '678 469 744'
      q4 = '1256 788 402'
      h = 0.00313
    [../]
  [../]
  [./disp_z]
    block = '1 2'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FluctuationsIC
      block = '1 2'
      epsilon = 1.0e-6
      q1 = '1179 521 286'
      q2 = '483 951 322'
      q3 = '723 800 1075'
      q4 = '1114 1269 261'
      h = 0.085
    [../]
  [../]
[]


[AuxVariables]
  [./chernsimonsdensity]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./curlPmag]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./stress_xx_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_xx
  [../]
  [./stress_yy_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_yy
  [../]
  [./stress_xy_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_xy
  [../]
  [./stress_xz_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_xz
  [../]
  [./stress_zz_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_zz
  [../]
  [./stress_yz_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_yz
  [../]
  [./strain_xx_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_xx
  [../]
  [./strain_yy_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_yy
  [../]
  [./strain_xy_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_xy
  [../]
  [./strain_xz_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_xz
  [../]
  [./strain_zz_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_zz
  [../]
  [./strain_yz_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_yz
  [../]
[]


[AuxKernels]
  [./cherm_sim]
    type = ChernSimonsDensity
    variable = chernsimonsdensity
    execute_on = 'timestep_end'
  [../]

  [./curl_P]
    type = CurlP
    variable = curlPmag
    execute_on = 'timestep_end'
  [../]

  [./matl_e11]
    type = RankTwoAux
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 0
    variable = strain_xx_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 1
    variable = strain_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e13]
    type = RankTwoAux
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 2
    variable = strain_xz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 1
    variable = strain_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e23]
    type = RankTwoAux
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 2
    variable = strain_yz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e33]
    type = RankTwoAux
    rank_two_tensor = total_strain
    index_i = 2
    index_j = 2
    variable = strain_zz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s11]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = stress_xx_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s12]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = stress_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s13]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
    variable = stress_xz_elastic
    execute_on = 'timestep_end'
  [../]
 [./matl_s22]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = stress_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s23]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    variable = stress_yz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s33]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    variable = stress_zz_elastic
    execute_on = 'timestep_end'
  [../]
[]

[Materials]
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
    block = '1'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    block = '1'
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
    block = '1'
  [../]

  [./slab_ferroelectric]
    block = '1'
    type = ComputeElectrostrictiveTensor
    Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.03375 0.03375 0.03375'
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
  [../]

  [./elasticity_tensor_2]
    type = ComputeElasticityTensor
    C_ijkl = '319 99.6 99.6 319 99.6 319 109.53 109.53 109.53'
    fill_method = symmetric9
    block = '2'
  [../]
  [./strain_2]
    type = ComputeSmallStrain
    block = '2'
  [../]
  [./stress_2]
    type = ComputeLinearElasticStress
    block = '2'
  [../]
[]


[Kernels]
  #Elastic problem
  [./TensorMechanics]
  #This is an action block
  [../]
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
##Polarization-strain coupling

  [./ferroelectriccouplingp_xx]
    block = '1'
    type = FerroelectricCouplingP
    variable = polar_x
    component = 0
  [../]
  [./ferroelectriccouplingp_yy]
    block = '1'
    type = FerroelectricCouplingP
    variable = polar_y
    component = 1
  [../]
  [./ferroelectriccouplingp_zz]
    block = '1'
    type = FerroelectricCouplingP
    variable = polar_z
    component = 2
  [../]

  [./ferroelectriccouplingX_xx]
    type = FerroelectricCouplingX
    block = '1'
    variable = disp_x
    component = 0
  [../]
  [./ferroelectriccouplingX_yy]
    type = FerroelectricCouplingX
    block = '1'
    variable = disp_y
    component = 1
  [../]
  [./ferroelectriccouplingX_zz]
    type = FerroelectricCouplingX
    block = '1'
    variable = disp_z
    component = 2
  [../]
  ##Electrostatics
  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
     permittivity = 0.08854187
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     permittivity = 0.08854187
  [../]

  [./DIE_E_int]
     type = Electrostatics
     variable = potential_int
     permittivity = 2.65626 # 300 * \epsilon_r
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

  #Short circuit BC
  [./potential_int_1]
    type = DirichletBC
    variable = potential_int
    boundary = '1'
    value = -0.00001
  [../]

  [./potential_int_2]
    type = DirichletBC
    variable = potential_int
    boundary = '2'
    value = -0.00001
  [../]

  #Let elastic deformation die off far from sphere
  [./sides_disp_x]
    variable = disp_x
    type = DirichletBC
    value = 0
    boundary = '1 2 3 4 5 6'
  [../]
  [./sides_disp_y]
    variable = disp_y
    type = DirichletBC
    value = 0
    boundary = '1 2 3 4 5 6'
  [../]
  [./sides_disp_z]
    variable = disp_z
    type = DirichletBC
    value = 0
    boundary = '1 2 3 4 5 6'
  [../]
[]



[Postprocessors]
   [./Fbulk]
      type = BulkEnergy
      block = '1'
      execute_on = 'timestep_end'
    [../]
    [./Fwall]
      type = WallEnergy
      block = '1'
      execute_on = 'timestep_end'
    [../]
    [./Felastic]
      type = ElasticEnergy
      block = '1'
      execute_on = 'timestep_end'
    [../]
    [./Fcoupled]
      block = '1'
      type = CoupledEnergy
      execute_on = 'timestep_end'
    [../]
    [./Felec]
      block = '1'
      type = ElectrostaticEnergy
      permittivity = 0.08854187
      execute_on = 'timestep_end'
    [../]
    [./Ftotal]
      type = TotalEnergyFlow
      Fbulk = Fbulk
      Fwall = Fwall
      Fcoupled = Fcoupled
      Felec = Felec
      execute_on = 'timestep_end'
    [../]
    [./perc_change]
     type = PercentChangePostprocessor
     postprocessor = Ftotal
   [../]
   [./num_NLin]
    type = NumNonlinearIterations
   [../]
  # [./num_Lin]
  #  type = NumLinearIterations
  # [../]
  []
[]


[UserObjects]
 [./kill]
  type = Terminator
  expression = 'perc_change <= 5.0e-5'
 [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121                1e-6      1e-8    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
    [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.35
    #iteration_window = 3
    optimal_iterations = 5 #should be 5 probably
    growth_factor = 1.4
    linear_iteration_ratio = 200
    cutback_factor =  0.8
[../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  #dt = 0.5
  dtmin = 1e-13
  dtmax = 0.35
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = outDef
    elemental_as_nodal = true
    interval = 1
  [../]
[]

