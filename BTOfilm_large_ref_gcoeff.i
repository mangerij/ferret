
[Mesh]
  file = outBTO_thinfilm_085_50_50_18_c001.e
[]

[GlobalParams]
  len_scale = 1.0
  alpha1 = -0.03
  alpha11 = -0.671
  alpha111 = 8.28
  alpha12 = 0.323
  alpha112 = 4.47
  alpha123 = 4.919
  G110 = 1.0
  G11_G110 = 0.51
  G12/G110 = 0.0
  G44/G110 = 0.02
  G44P/G110 = 0.02
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  displacements = 'disp_x disp_y disp_z'
  prefactor = -0.002

  n_a = 2.42 #cubic refractive index of BTO
  n_b = 2.42
  n_g = 2.42
[]




[Variables]
  [./polar_x]
    block = '1'
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_x
    initial_from_file_timestep = 'LATEST'
  [../]
  [./polar_y]
    block = '1'
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_y
    initial_from_file_timestep = 'LATEST'
  [../]
  [./polar_z]
    block = '1'
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_z
    initial_from_file_timestep = 'LATEST'
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
    initial_from_file_var = potential_int
    initial_from_file_timestep = 'LATEST'
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = disp_x
    initial_from_file_timestep = 'LATEST'
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = disp_y
    initial_from_file_timestep = 'LATEST'
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = disp_z
    initial_from_file_timestep = 'LATEST'
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

  [./dn_1]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./dn_2]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./dn_3]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./dn_bire_12]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./dn_bire_23]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./dn_bire_13]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./n_1]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]

  [./n_2]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]

  [./n_3]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
[]


[AuxKernels]
  [./matl_e11]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    variable = strain_xx_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    variable = strain_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e13]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 2
    variable = strain_xz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    variable = strain_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e23]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 2
    variable = strain_yz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e33]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
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


 [./dn_s1]
    type = ChangeInRefractiveIndexWithGCoeffPolar
    index_i = 0
    index_j = 0
    index_k = 0
    index_l = 0
    variable = dn_1
    execute_on = 'timestep_end'
  [../]

  [./dn_s2]
    type = ChangeInRefractiveIndexWithGCoeffPolar
    index_i = 1
    index_j = 1
    index_k = 1
    index_l = 1
    variable = dn_2
    execute_on = 'timestep_end'
  [../]

  [./dn_s3]
    type = ChangeInRefractiveIndexWithGCoeffPolar
    index_i = 2
    index_j = 2
    index_k = 2
    index_l = 2
    variable = dn_3
    execute_on = 'timestep_end'
  [../]

  [./dn_bire_s12]
    type = Birefringence
    variable = dn_bire_12
    per1 = n_1
    per2 = n_2
    execute_on = 'timestep_end'
  [../]

  [./dn_bire_s23]
    type = Birefringence
    variable = dn_bire_23
    per1 = n_2
    per2 = n_3
    execute_on = 'timestep_end'
  [../]

  [./dn_bire_s13]
    type = Birefringence
    variable = dn_bire_13
    per1 = n_1
    per2 = n_3
    execute_on = 'timestep_end'
  [../]


  [./n_1_c]
    type = RefractiveIndex
    variable = n_1
    component = 0
    var1 = dn_1
    execute_on = 'timestep_end'
  [../]

  [./n_2_c]
    type = RefractiveIndex
    variable = n_2
    component = 1
    var1 = dn_2
    execute_on = 'timestep_end'
  [../]

  [./n_3_c]
    type = RefractiveIndex
    variable = n_3
    component = 2
    var1 = dn_3
    execute_on = 'timestep_end'
  [../]

[]

[Materials]
  [./eigen_strain_zz] #Use for stress-free strain (ie epitaxial)
    type = ComputeEigenstrain
    block = '1'
    # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
    eigen_base = '1 0 0 0 1 0 0 0 0'
    eigenstrain_name = eigenstrain
  [../]

  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #from MaterialsProject
    C_ijkl = '282.23 103.7 103.7 282.23 103.7 282.23 120.24 120.24 120.24'
    block = '1'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    block = '1'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
    block = '1'
  [../]

  [./slab_ferroelectric]
    block = '1'
    type = ComputeElectrostrictiveTensor
    Q_mnkl = '0.11 -0.045 -0.045 0.11 -0.045 0.11 0.059 0.059 0.059'
    C_ijkl = '282.23 103.7 103.7 282.23 103.7 282.23 120.24 120.24 120.24'
    compute_electrostrictive_coeff = true
  [../]
  [./photoelastic_tensor_1]
    type = ComputeElastoopticTensor
    fill_method = symmetric9
    # Use BaTiO3, use cubic crystal symmetry. see J. Appl. Phys., 78, 4, 1995
    P_mnkl = '0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    #P_mnkl = '0.37 0.11 0.11 0.37 0.11 0.37 -0.3 -0.3 -0.3'
    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
    block = '1'
  [../]
  [./gtensor_1]
    type = ComputeGCoeffTensor
    fill_method = symmetric9
    # Use BaTiO3, use cubic crystal symmetry. see J. Appl. Phys., 78, 4, 1995 (converted to nm^4/aC^2
    g_ijkl = '0.17.0 0.013 0.013 0.17 0.013 0.17 0.062 0.062 0.062'
    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
    block = '1'
  [../]
  [./beta_tensor_1]
    type = ComputeIndicatrix
    block = '1'
    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]
  [./delta_beta_tensor_1]
    type = ComputeDeltaIndicatrix
    block = '1'
  [../]
  [./PO_tensor_1]
    type = ComputePolarOpticGCoeffTensor
    block = '1'
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
     permittivity = 0.0664064025 #er = 7.5
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     block = '1'
     permittivity = 0.0664064025
  [../]

  [./DIE_E_int]
     type = Electrostatics
     variable = potential_int
     block  = '2'
     permittivity = 2.6562561
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
  [./bot_disp_x]
    variable = disp_x
    type = DirichletBC
    value = 0
    boundary = '2'
  [../]
  [./bot_disp_y]
    variable = disp_y
    type = DirichletBC
    value = 0
    boundary = '2'
  [../]
  [./bot_disp_z]
    variable = disp_z
    type = DirichletBC
    value = 0
    boundary = '2'
  [../]

  [./bot_potential_int]
    variable = potential_int
    type = DirichletBC
    value = 0.00001
    boundary = '2'
  [../]
  [./top_potential_int]
    variable = potential_int
    type = DirichletBC
    value = 0.00001
    boundary = '7'
  [../]

  [./Periodic]
    [./TB_disp_x_pbc]
      variable = disp_x
      primary = '9'
      secondary = '11'
      translation = '0 50 0'
    [../]
    [./TB_disp_y_pbc]
      variable = disp_y
      primary = '9'
      secondary = '11'
      translation = '0 50 0'
    [../]
    [./TB_disp_z_pbc]
      variable = disp_z
      primary = '9'
      secondary = '11'
      translation = '0 50 0'
    [../]

    [./TB_polar_x_pbc]
      variable = polar_x
      primary = '9'
      secondary = '11'
      translation = '0 50 0'
    [../]
    [./TB_polar_y_pbc]
      variable = polar_y
      primary = '9'
      secondary = '11'
      translation = '0 50 0'
    [../]
    [./TB_polar_z_pbc]
      variable = polar_z
      primary = '9'
      secondary = '11'
      translation = '0 50 0'
    [../]
    [./TB_potential_int_pbc]
      variable = potential_int
      primary = '9'
      secondary = '11'
      translation = '0 50 0'
    [../]
  #
    [./TBsub_disp_x_pbc]
      variable = disp_x
      primary = '3'
      secondary = '5'
      translation = '0 50 0'
    [../]
    [./TBsub_disp_y_pbc]
      variable = disp_y
      primary = '3'
      secondary = '5'
      translation = '0 50 0'
    [../]
    [./TBsub_disp_z_pbc]
      variable = disp_z
      primary = '3'
      secondary = '5'
      translation = '0 50 0'
    [../]

    [./RL_disp_x_pbc]
      variable = disp_x
      primary = '10'
      secondary = '12'
      translation = '50 0 0'
    [../]
    [./RL_disp_y_pbc]
      variable = disp_y
      primary = '10'
      secondary = '12'
      translation = '50 0 0'
    [../]
    [./RL_disp_z_pbc]
      variable = disp_z
      primary = '10'
      secondary = '12'
      translation = '50 0 0'
    [../]

    [./RL_polar_x_pbc]
      variable = polar_x
      primary = '10'
      secondary = '12'
      translation = '50 0 0'
    [../]
    [./RL_polar_y_pbc]
      variable = polar_y
      primary = '10'
      secondary = '12'
      translation = '50 0 0'
    [../]
    [./RL_polar_z_pbc]
      variable = polar_z
      primary = '10'
      secondary = '12'
      translation = '50 0 0'
    [../]
    [./RL_potential_int_pbc]
      variable = potential_int
      primary = '10'
      secondary = '12'
      translation = '50 0 0'
    [../]

    [./RLsub_disp_x_pbc]
      variable = disp_x
      primary = '4'
      secondary = '6'
      translation = '50 0 0'
    [../]
    [./RLsub_disp_y_pbc]
      variable = disp_y
      primary = '4'
      secondary = '6'
      translation = '50 0 0'
    [../]
    [./RLsub_disp_z_pbc]
      variable = disp_z
      primary = '4'
      secondary = '6'
      translation = '50 0 0'
    [../]
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
      block = '1 2'
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
     execute_on = 'timestep_end'
   [../]
[]


[UserObjects]
 [./kill]
  type = Terminator
  expression = 'perc_change <= 1.0e-3'
  execute_on = 'timestep_end'
 [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    240               1e-10       1e-8      1e-6    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
    [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    #iteration_window = 3
    optimal_iterations = 5 #should be 5 probably
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.92
[../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 0.8
  num_steps = 2
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    execute_on = 'timestep_end'
    file_base = outBTO_thinfilm_gcoeff_2
    elemental_as_nodal = true
    interval = 1
  [../]
  [./outcsv]
    type = CSV
    file_base = outBTO_thinfilm_gcoeff_2
    execute_on = 'timestep_end'
  [../]
[]

