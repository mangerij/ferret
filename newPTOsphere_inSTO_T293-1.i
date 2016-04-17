
[Mesh]
  file = embedded_single_sphere_1.e

  #sphere = block 1
  #dielectric = block 2
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
  Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.03375 0.03375 0.03375'
  #permittivity = 5.8
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
  #potential_ext = potential_ext
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  displacements = 'disp_x disp_y disp_z'
  #use_displaced_mesh = false
[]



[Variables]
  [./polar_x]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      seed = 1
      min = -0.5e-5
      max = 0.5e-5
    [../]
    #[./InitialCondition]
    #  type = FluctuationsIC
    #  epsilon = 1.0e-5
    #  q1 = '1011 1184 1329'
    #  q2 = '512 586 1106'
    #  q3 = '547 1297 962'
    #  q4 = '511 833 371'
    #  h = 1.35
    #[../]
  [../]
  [./polar_y]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      seed = 1
      min = -0.5e-5
      max = 0.5e-5
    [../]
    #[./InitialCondition]
    #  type = FluctuationsIC
    #  epsilon = 1.0e-5
    #  q1 = '678 469 744'
    #  q2 = '328 942 895'
    #  q3 = '980 1121 343'
    #  q4 = '1234 970 593'
    #  h = 2.48
    #[../]
  [../]
  [./polar_z]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      seed = 1
      min = -0.5e-5
      max = 0.5e-5
    [../]
    #[./InitialCondition]
    #  type = FluctuationsIC
    #  epsilon = 1.0e-5
    #  q1 = '602 666 538'
    #  q2 = '678 469 744'
    #  q3 = '1107 529 1175'
    #  q4 = '1133 1108 532'
    #  h = 1.08
    #[../]
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
    [./InitialCondition]
      type = RandomIC
      seed = 1
      min = -0.5e-5
      max = 0.5e-5
    [../]
    #[./InitialCondition]
    #  type = FluctuationsIC
    #  block = '1 2'
    #  epsilon = 1.0e-5
    #  q1 = '354 1005 645'
    #  q2 = '715 1065 1132'
    #  q3 = '391 305 1106'
    #  q4 = '1053 1116 627'
    #  h = 0.25
    #[../]
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      seed = 1
      min = -0.5e-5
      max = 0.5e-5
    [../]
    #[./InitialCondition]
    #  block = '1'
    #  type = FluctuationsIC
    #  epsilon = 1.0e-5
    #  q1 = '256 1246 730'
    #  q2 = '1234 770 1257'
    #  q3 = '1262 279 467'
    #  q4 = '467 944 1326'
    #  h = -0.017
    #[../]
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      seed = 1
      min = -0.5e-5
      max = 0.5e-5
    [../]
    #[./InitialCondition]
    #  type = FluctuationsIC
    #  block = '1'
    #  epsilon = 1.0e-5
    #  q1 = '1522 415 1009'
    #  q2 = '147 622 351'
    #  q3 = '678 469 744'
    #  q4 = '1256 788 402'
    #  h = 0.0313
    #[../]
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      seed = 1
      min = -0.5e-5
      max = 0.5e-5
    [../]
    #[./InitialCondition]
    #  type = FluctuationsIC
    #  block = '1'
    #  epsilon = 1.0e-6
    #  q1 = '1179 521 286'
    #  q2 = '483 951 322'
    #  q3 = '723 800 1075'
    #  q4 = '1114 1269 261'
    #  h = 0.085
    #[../]
  [../]
[]

[AuxVariables]
  [./stress_xx_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xx_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./curlmag]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./csnumber]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]

  [./curlP]
    type = CurlP
    variable = curlmag
    execute_on = 'timestep_end'
  [../]

  [./CSnumber]
    type = ChernSimonsDensity
    variable = csnumber
    execute_on = 'timestep_end'
  [../]


  [./matl_s11]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    block = '1 2'
    variable = stress_xx_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s12]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    block = '1 2'
    variable = stress_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s13]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
    block = '1 2'
    variable = stress_xz_elastic
    execute_on = 'timestep_end'
  [../]
 [./matl_s22]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    block = '1 2'
    variable = stress_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s23]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    block = '1 2'
    variable = stress_yz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s33]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    block = '1 2'
    variable = stress_zz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e11]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    block = '1 2'
    variable = strain_xx_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    block = '1 2'
    variable = strain_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e13]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 2
    block = '1 2'
    variable = strain_xz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    block = '1 2'
    variable = strain_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e23]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 2
    block = '1 2'
    variable = strain_yz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e33]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 2
    index_j = 2
    block = '1 2'
    variable = strain_zz_elastic
    execute_on = 'timestep_end'
  [../]
[]

[Kernels]

    #Elastic problem
    [./TensorMechanics]
       #This is an action block but adds the appropriate elastic problem kernels (stress divergence)
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
       type=WallEnergyDerivative
       variable = polar_x
       component = 0
    [../]
   [./walled_y]
       type=WallEnergyDerivative
       variable = polar_y
       component = 1
    [../]
    [./walled_z]
       type=WallEnergyDerivative
       variable = polar_z
       component = 2
    [../]

    ##Polarization-strain coupling
    [./ferroelectriccouplingp_xx]
       type = FerroelectricCouplingP
       variable = polar_x
       block = '1'
       component = 0
    [../]
    [./ferroelectriccouplingp_yy]
       type = FerroelectricCouplingP
       variable = polar_y
       block = '1'
       component = 1
    [../]
    [./ferroelectriccouplingp_zz]
       type = FerroelectricCouplingP
       variable = polar_z
       block = '1'
       component = 2
    [../]


    ##Electrostatics
    [./polar_x_electric_E]
       type=PolarElectricEStrong
       variable = potential_int
       permittivity = 0.0885
       block = '1'
    [../]
    [./FE_E_int]
       type=Electrostatics
       variable = potential_int
       block = '1'
       permittivity = 0.0885
    [../]
    [./DIE_E_int]
       type=Electrostatics
       variable = potential_int
       block = '2'
       permittivity = 2.64
    [../]
    [./polar_electric_px]
       type=PolarElectricPStrong
       variable = polar_x
       component = 0
    [../]
    [./polar_electric_py]
       type=PolarElectricPStrong
       variable = polar_y
       component = 1
    [../]
    [./polar_electric_pz]
       type=PolarElectricPStrong
       variable = polar_z
       component = 2
    [../]

    ##Time dependence
    [./polar_x_time]
      type=TimeDerivativeScaled
      variable=polar_x
      # time_scale = 1.8e-7
      time_scale = 1.0
    [../]
    [./polar_y_time]
      type=TimeDerivativeScaled
      variable=polar_y
      # time_scale = 1.8e-7
      time_scale = 1.0
    [../]
    [./polar_z_time]
      type=TimeDerivativeScaled
      variable = polar_z
      # time_scale = 1.8e-7
      time_scale = 1.0
  [../]
[]



[BCs]

  #Short circuit boundary condition:

  [./potential_int_1]
    type = DirichletBC
    variable = potential_int
    boundary = '1'
    value = 0.0001
  [../]
  [./potential_int_2]
    type = DirichletBC
    variable = potential_int
    boundary = '2'
    value = 0.0001
  [../]




  [./disp_x_1]
    type = DirichletBC
    variable = disp_x
    boundary = '1'
    value = 0
  [../]
  [./disp_y_1]
    type = DirichletBC
    variable = disp_y
    boundary = '1'
    value = 0
  [../]
  [./disp_z_1]
    type = DirichletBC
    variable = disp_z
    boundary = '1'
    value = 0
  [../]

  [./disp_x_2]
    type = DirichletBC
    variable = disp_x
    boundary = '2'
    value = 0
  [../]
  [./disp_y_2]
    type = DirichletBC
    variable = disp_y
    boundary = '2'
    value = 0
  [../]
  [./disp_z_2]
    type = DirichletBC
    variable = disp_z
    boundary = '2'
    value = 0
  [../]

  [./disp_x_3]
    type = DirichletBC
    variable = disp_x
    boundary = '3'
    value = 0
  [../]
  [./disp_y_3]
    type = DirichletBC
    variable = disp_y
    boundary = '3'
    value = 0
  [../]
  [./disp_z_3]
    type = DirichletBC
    variable = disp_z
    boundary = '3'
    value = 0
  [../]

  [./disp_x_4]
    type = DirichletBC
    variable = disp_x
    boundary = '4'
    value = 0
  [../]
  [./disp_y_4]
    type = DirichletBC
    variable = disp_y
    boundary = '4'
    value = 0
  [../]
  [./disp_z_4]
    type = DirichletBC
    variable = disp_z
    boundary = '4'
    value = 0
  [../]

  [./disp_x_5]
    type = DirichletBC
    variable = disp_x
    boundary = '5'
    value = 0
  [../]
  [./disp_y_5]
    type = DirichletBC
    variable = disp_y
    boundary = '5'
    value = 0
  [../]
  [./disp_z_5]
    type = DirichletBC
    variable = disp_z
    boundary = '5'
    value = 0
  [../]

  [./disp_x_6]
    type = DirichletBC
    variable = disp_x
    boundary = '6'
    value = 0
  [../]
  [./disp_y_6]
    type = DirichletBC
    variable = disp_y
    boundary = '6'
    value = 0
  [../]
  [./disp_z_6]
    type = DirichletBC
    variable = disp_z
    boundary = '6'
    value = 0
  [../]

[]

[Materials]
    #amorphous Si #C_ijkl = '62.543 5.735 5.735 62.543 5.735 62.543 28.404 28.404 28.404'
    #STO #C_ijkl = '319 99.6 99.6 319 99.6 319 109.53 109.53 109.53' #this is averaged..
    #PTO #C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
    #BTO #C_ijkl = '178 96.4 96.4 178 96.4 178 122 122 122'
  [./spheres_ferroelectric]
    type=LinearFerroelectricMaterial
    block = '1'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
    #euler_angle_1 = 0.0 #currently will only rotate C_ijkl but Q_ijkl and C_ijkl are to be collinear
    #euler_angle_2 = 0.0
    #euler_angle_3 = 0.0
  [../]

  [./elasticity_tensor1]
    type = ComputeElasticityTensor
    block = '1'
    fill_method = symmetric9
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
  [../]
  [./strain1]
    type = ComputeSmallStrain
    block = '1'
  [../]
  [./stress1]
    type = ComputeLinearElasticStress
    block = '1'
  [../]

  [./elasticity_tensor2]
    type = ComputeElasticityTensor
    block = '2'
    fill_method = symmetric9
    C_ijkl = '319 99.6 99.6 319 99.6 319 109.53 109.53 109.53'
  [../]
  [./strain2]
    type = ComputeSmallStrain
    block = '2'
  [../]
  [./stress2]
    type = ComputeLinearElasticStress
    block = '2'
  [../]
[]

[Postprocessors]
  [./bulk_energy]
     type = BulkEnergy
     block = '1'
     execute_on = 'timestep_end'
    [../]
   [./wall_energy]
    type = WallEnergy
    block = '1'
    execute_on = 'timestep_end'
   [../]
    [./elastic_energy]
    type = ElasticEnergy
    block = '1'
    execute_on = 'timestep_end'
    [../]
    [./coupled_energy]
     block = '1'
    type = CoupledEnergy
    execute_on = 'timestep_end'
    [../]
    [./electrostatic_energy]
     block = '1'
     type = ElectrostaticEnergy
     execute_on = 'timestep_end'
     permittivity = 0.08854187
    [../]
    [./total_energy_noelastic]
      type = TotalEnergyFlow
      bulk_energy = bulk_energy
      wall_energy = wall_energy
      bulk_energy_fourth = bulk_energy_fourth
      coupled_energy = coupled_energy
      electrostatic_energy = electrostatic_energy
      execute_on = 'timestep_end'
    [../]
    [./perc_change]
     type = PercentChangePostprocessor
     postprocessor = total_energy_noelastic
   [../]
  []



[UserObjects]
 [./kill]
  type = Terminator
  expression = 'perc_change <= 1.0e-4'
 [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -options_left'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type '
    petsc_options_value = '    121               1e-8      1e-8     bjacobi '
  [../]
[]

[Executioner]
  type = Transient
    [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.15
    #iteration_window = 3
    optimal_iterations = 5
    growth_factor = 1.2
    linear_iteration_ratio = 1000
    cutback_factor =  0.95
[../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 0.2
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_PTOsphere_inSTO_1
    elemental_as_nodal = true
    interval = 8
  [../]
[]
