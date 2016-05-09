
[Mesh]
  file = exodus_thinfilm_12_8_10.e
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
  #permittivity = 5.8
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
  #potential_ext = potential_ext
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  kappa = 1.0
  displacements = 'disp_x disp_y disp_z'
  #artificial = 1.0 #this is an artificial scaling parameter for coupling.
  #use_displaced_mesh = false
  prefactor = 0.005 #negative = tension, positive = compression for stress-free strain along Z
[]



[Variables]
  [./polar_x]
    block = '1'
    order = FIRST  #comment: Third order Hermites ->? dense problem -> pc fails
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./polar_y]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./polar_z]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-6
      max = 0.5e-6
    [../]
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-6
      max = 0.5e-6
    [../]
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-6
      max = 0.5e-6
    [../]
  [../]
[]


[AuxVariables]
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
  [./wallenergydensity]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]


[AuxKernels]

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
  [./wallenergydensity]
    type = WallEnergyDensity
    block = '1'
    variable = wallenergydensity
    execute_on = 'timestep_end'
  [../]
[]


[Kernels]
  #Elastic problem
  [./TensorMechanics]
    block = '1 2'
  #This is an action block
  [../]

  #[./sd_x]
  #  type = ModifiedStressDivergenceTensors
  #  variable = disp_x
  #  component = 0
  #  block = '1'
  #[../]
  #[./sd_y]
  #  type = ModifiedStressDivergenceTensors
  #  variable = disp_y
  #  component = 1
  #  block = '1'
  #[../]
  #[./sd_z]
  #  type = ModifiedStressDivergenceTensors
  #  variable = disp_z
  #  component = 2
  #  block = '1'
  #[../]


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
    type = FerroelectricCouplingP
    variable = polar_x
    component = 0
  [../]
  [./ferroelectriccouplingp_yy]
    type = FerroelectricCouplingP
    variable = polar_y
    component = 1
  [../]
  [./ferroelectriccouplingp_zz]
    type = FerroelectricCouplingP
    variable = polar_z
    component = 2
  [../]

  [./ferroelectriccouplingq_xx]
    type = FerroelectricCouplingQ
    variable = polar_x
    component = 0
  [../]
  [./ferroelectriccouplingq_yy]
    type = FerroelectricCouplingQ
    variable = polar_y
    component = 1
  [../]
  [./ferroelectriccouplingq_zz]
    type = FerroelectricCouplingQ
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

  #[./disp_x_time]
  #   type = TimeDerivativeScaled
  #   variable=disp_x
  #  time_scale = 1.0
  #[../]
  #[./disp_y_time]
  #   type = TimeDerivativeScaled
  #   variable=disp_y
  #  time_scale = 1.0
  #[../]
  #[./disp_z_time]
  #   type = TimeDerivativeScaled
  #   variable = disp_z
  #  time_scale = 1.0
  #[../]
[]


[BCs]

  #[./disp_x_SF]
  #  type = StressFreeBC
  #  variable = disp_x
  #  component = 0
  #  boundary = '1'
  #
  #[../]
  #[./disp_y_1]
  #  type = StressFreeBC
  #  variable = disp_y
  #  component = 1
  #  boundary = '1'
  #
  #[../]
  #[./disp_z_1]
  #  type = StressFreeBC
  #  variable = disp_z
  #  component = 2
  #  boundary = '1'
  #[../]

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

  [./bot_disp_x]
    variable = disp_x
    type = DirichletBC
    value = 0
    boundary = '7'
  [../]
  [./bot_disp_y]
    variable = disp_y
    type = DirichletBC
    value = 0
    boundary = '7'
  [../]
  [./bot_disp_z]
    variable = disp_z
    type = DirichletBC
    value = 0
    boundary = '7'
  [../]

  [./Periodic]
    [./TB_disp_x_pbc]
      variable = disp_x
      primary = '3'
      secondary = '5'
      translation = '0 12 0'
    [../]
    [./TB_disp_y_pbc]
      variable = disp_y
      primary = '3'
      secondary = '5'
      translation = '0 12 0'
    [../]
    [./TB_disp_z_pbc]
      variable = disp_z
      primary = '3'
      secondary = '5'
      translation = '0 12 0'
    [../]

    [./TB_polar_x_pbc]
      variable = polar_x
      primary = '3'
      secondary = '5'
      translation = '0 12 0'
    [../]
    [./TB_polar_y_pbc]
      variable = polar_y
      primary = '3'
      secondary = '5'
      translation = '0 12 0'
    [../]
    [./TB_polar_z_pbc]
      variable = polar_z
      primary = '3'
      secondary = '5'
      translation = '0 12 0'
    [../]
    [./TB_potential_int_pbc]
      variable = potential_int
      primary = '3'
      secondary = '5'
      translation = '0 12 0'
    [../]
  ##
    [./TBsub_disp_x_pbc]
      variable = disp_x
      primary = '8'
      secondary = '10'
      translation = '0 12 0'
    [../]
    [./TBsub_disp_y_pbc]
      variable = disp_y
      primary = '8'
      secondary = '10'
      translation = '0 12 0'
    [../]
    [./TBsub_disp_z_pbc]
      variable = disp_z
      primary = '8'
      secondary = '10'
      translation = '0 12 0'
    [../]

    [./RL_disp_x_pbc]
      variable = disp_x
      primary = '4'
      secondary = '6'
      translation = '12 0 0'
    [../]
    [./RL_disp_y_pbc]
      variable = disp_y
      primary = '4'
      secondary = '6'
      translation = '12 0 0'
    [../]
    [./RL_disp_z_pbc]
      variable = disp_z
      primary = '4'
      secondary = '6'
      translation = '12 0 0'
    [../]

    [./RL_polar_x_pbc]
      variable = polar_x
      primary = '4'
      secondary = '6'
      translation = '12 0 0'
    [../]
    [./RL_polar_y_pbc]
      variable = polar_y
      primary = '4'
      secondary = '6'
      translation = '12 0 0'
    [../]
    [./RL_polar_z_pbc]
      variable = polar_z
      primary = '4'
      secondary = '6'
      translation = '12 0 0'
    [../]
    [./RL_potential_int_pbc]
      variable = potential_int
      primary = '4'
      secondary = '6'
      translation = '12 0 0'
    [../]

    [./RLsub_disp_x_pbc]
      variable = disp_x
      primary = '9'
      secondary = '11'
      translation = '12 0 0'
    [../]
    [./RLsub_disp_y_pbc]
      variable = disp_y
      primary = '9'
      secondary = '11'
      translation = '12 0 0'
    [../]
    [./RLsub_disp_z_pbc]
      variable = disp_z
      primary = '9'
      secondary = '11'
      translation = '12 0 0'
    [../]
  [../]
[]

[Materials]
   [./eigen_strain_xx_yy] #Use for stress-free strain (which approach is correct?)
    type = ComputeEigenstrain
    block = '1'
   # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
    eigen_base = '0 0 0 0 0 0 0 0 1'
  [../]
 # [./eigen_strain_xx_yy] #Use for stress-free strain (ie epitaxial)
 #  type = ComputeEigenstrain
 #  block = '1'
 # # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
 #  eigen_base = '1 0 0 0 1 0 0 0 0'
 #[../]

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

[Postprocessors]
#  [./volume]
#    type = VolumePostprocessor
#    block = '1'
#    use_displaced_mesh = true
#  [../]
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
   [./nodes]
     type = NumNodes
    [../]
   [./num_NLin]
    type = NumNonlinearIterations
   [../]
   [./num_Lin]
    type = NumLinearIterations
   [../]
[]

#[Adaptivity]
#  marker = combo
#  initial_steps = 0
#  max_h_level = 3
#  [./Markers]
#    [./EFM_3]
#      type = ErrorFractionMarker
#      coarsen = 0.025
#      refine = 0.45
#      indicator = GJI_3
#    [../]
#
#    [./combo]
#      type = ComboMarker
#      markers = 'EFM_3'
#    [../]
#  [../]
#
#  [./Indicators]
#    [./GJI_3]
#     type = GradientJumpIndicator
#     variable = potential_int
#    [../]
#  [../]
#[]

[UserObjects]
 [./kill]
  type = Terminator
  expression = 'perc_change <= 1.0e-5'
 [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -ksp_snes_ew'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121                1e-8      1e-8    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.3
    optimal_iterations = 6 #should be 5 probably or 1?
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.55
  [../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  #dt = 0.5
  dtmin = 0.05
  dtmax = 0.3
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out1]
    type = Exodus
    file_base = outPTO_TF_c05_modSF_E_89_oneeightminus_CTD
    elemental_as_nodal = true
    interval = 1
  [../]

  #outPTO_TF_c05_modSF_E_2 been diverging at step 87 if we don't let dt decrease so that E flattens

  #[./out2] #need for AMR -- [postprocessors don't print correctly in .e files]
  #  type = CSV
  #  file_base = outPTO_thinfilm_c05_modSF_AMRcsv_3200
  #  interval = 1
  #[../]
[]
