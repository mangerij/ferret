[Mesh]
  file = supercell.e
[]

[GlobalParams]
  len_scale = 1.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_E_int = potential_E_int

  displacements = 'u_x u_y u_z'

  ##############################################
  ##=
  ##  IMPORTANT(!): Units in Ferret are nm, kg,
  ##                seconds, and attocoulombs
  ##
  ##############################################

  u_x = u_x
  u_y = u_y
  u_z = u_z
[]

[Variables]

  #################################
  ##
  ##  Variable definitions
  ##    P, u, phi
  ##  and their initial conditions
  ##
  #################################

  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-4
      max = 0.1e-4
      seed = 1
    [../]
    block = '1 2'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-4
      max = 0.1e-4
      seed = 2
    [../]
    block = '1 2'
  [../]

  [./potential_E_int]
    order = FIRST
    family = LAGRANGE
    block = '1 2 3'
  [../]

  [./u_x]
    order = FIRST
    family = LAGRANGE
    block = '1 2 3'
  [../]
  [./u_y]
    order = FIRST
    family = LAGRANGE
    block = '1 2 3'
  [../]
  [./u_z]
    order = FIRST
    family = LAGRANGE
    block = '1 2 3'
  [../]
[]

[ICs]
    [./py1]
      type = RandomIC
      variable = polar_y
      min = 0.58
      max = 0.6
      seed = 3
      block = '1'
    [../]
    [./py2]
      type = RandomIC
      variable = polar_y
      min = -0.6
      max = -0.58
      seed = 3
      block = '2'
    [../]
[]

[AuxVariables]

  ######################################
  ##
  ##  Auxiarilly variable definitions
  ##   (can be intermediate variables
  ##   or for postprocessed quantities)
  ##
  ######################################


  ######################################
  ##
  ##  Stress/strain tensor components
  ##
  ######################################

  [./e00]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3'
  [../]
  [./e01]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3'
  [../]
  [./e10]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3'
  [../]
  [./e11]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3'
  [../]
  [./e12]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3'
  [../]
  [./e22]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3'
  [../]

  [./s00]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3'
  [../]
  [./s01]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3'
  [../]
  [./s10]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3'
  [../]
  [./s11]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3'
  [../]
  [./s12]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3'
  [../]
  [./s22]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2 3'
  [../]

[]

[AuxKernels]

  ######################################
  ##
  ##  Auxiarilly Kernel definitions
  ##   (can be intermediate "operations"
  ##   or for postprocessed quantities)
  ##
  ######################################

  [./e00]
    type = RankTwoAux
    variable = e00
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 0
  [../]
  [./e01]
    type = RankTwoAux
    variable = e01
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 1
  [../]
  [./e10]
    type = RankTwoAux
    variable = e10
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 0
  [../]
  [./e12]
    type = RankTwoAux
    variable = e12
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 2
  [../]
  [./e11]
    type = RankTwoAux
    variable = e11
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 1
  [../]
  [./e22]
    type = RankTwoAux
    variable = e22
    rank_two_tensor = total_strain
    index_i = 2
    index_j = 2
  [../]

  [./s00]
    type = RankTwoAux
    variable = s00
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
  [../]
  [./s01]
    type = RankTwoAux
    variable = s01
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
  [../]
  [./s10]
    type = RankTwoAux
    variable = s10
    rank_two_tensor = stress
    index_i = 1
    index_j = 0
  [../]
  [./s12]
    type = RankTwoAux
    variable = s12
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
  [../]
  [./s11]
    type = RankTwoAux
    variable = s11
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
  [../]
  [./s22]
    type = RankTwoAux
    variable = s22
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
  [../]
[]

[Materials]

  #################################################
  ##
  ##
  ## NOTE: there might be some Legendre transforms
  ##        depending on what approach you use
  ##        -i.e. inhomogeneous strain vs
  ##            homogeneous strain [renormalized]
  ##
  ##################################################

  [./Landau_P_FE]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '-0.1722883 -0.073 0.75 0.26 0.61 -3.67 0.0 0.0 0.0 0.0'  #corresponds to T = 298K
    block = '1 2'
  [../]

  [./Landau_G_FE]
    type = GenericConstantMaterial
    prop_names = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '0.173 0.6 0.0 0.3 0.3'
    block = '1 2'
  [../]

  [./mat_C_FE]
    type = GenericConstantMaterial
    prop_names = 'C11 C12 C44'
    prop_values = '175.0 79.4 111.1'
    block = '1 2'
  [../]
  [./mat_C_sub]
    type = GenericConstantMaterial
    prop_names = 'C11 C12 C44'
    prop_values = '121.0 42.0 37.0'
    block = '3'
  [../]

  ##################################################
  ##=
  ## NOTE: Sign convention in Ferret for the
  ##        electrostrictive coeff. is multiplied by
  ##        an overall factor of (-1)
  ##
  ##################################################

  [./mat_Q]
    type = GenericConstantMaterial
    prop_names = 'Q11 Q12 Q44'
    prop_values = '-0.089 0.026 -0.0084375'
    block = '1 2'
  [../]

  [./mat_q]
    type = GenericConstantMaterial
    prop_names = 'q11 q12 q44'
    prop_values = '-11.4 -0.01438 -7.5'
    block = '1 2'
  [../]

  [./eigen_strain]
    type = ComputeEigenstrain
    # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
    eigen_base = '1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0'
    eigenstrain_name = eigenstrain
    prefactor = 0.0
  [../]

  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9

   ###############################################
   ##
   ## symmetric9 fill_method is (default)
   ##     C11 C12 C13 C22 C23 C33 C44 C55 C66
   ##
   ###############################################

    C_ijkl = '175.0 79.4 79.4 175.0 79.4 175.0 111.1 111.1 111.1'
    block = '1 2'
  [../]
  [./elasticity_tensor_2]
    type = ComputeElasticityTensor
    fill_method = symmetric9

   ###############################################
   ##
   ## symmetric9 fill_method is (default)
   ##     C11 C12 C13 C22 C23 C33 C44 C55 C66
   ##
   ###############################################

    C_ijkl = '121.0 0.0 0.0 121.0 0.0 121.0 0.0 0.0 0.0'
    block = '3'
  [../]


  [./strain_12]
    type = ComputeSmallStrain
    eigenstrain_names = eigenstrain
  [../]
  [./stress_12]
    type = ComputeLinearElasticStress
  [../]


  [./permitivitty_1]

    ###############################################
    ##
    ##  so-called background dielectric constant
    ##  (it encapsulates the motion of core electrons
    ##  at high frequency) = e_b*e_0 (here we use
    ##  e_b = 10), see PRB. 74, 104014, (2006)
    ##
    ###############################################

    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '0.08854187'
    block = '1 2'
  [../]

  [./permitivitty_2]

    ###############################################
    ##
    ##  so-called background dielectric constant
    ##  (it encapsulates the motion of core electrons
    ##  at high frequency) = e_b*e_0 (here we use
    ##  e_b = 10), see PRB. 74, 104014, (2006)
    ##
    ###############################################

    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '1.4'
    block = '3'
  [../]
[]


[Kernels]

  ###############################################
  ##
  ## Physical Kernel operators
  ## to enforce TDLGD evolution
  ##
  ###############################################


  #Elastic problem
  [./TensorMechanics]
    use_displaced_mesh = false
    eigenstrain_names = eigenstrain
  [../]

  [./bed_x]
    type = BulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
    block = '1 2'
  [../]
  [./bed_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
    block = '1 2'
  [../]
  [./bed_z]
    type = BulkEnergyDerivativeEighth
    variable = polar_z
    component = 2
    block = '1 2'
  [../]

  [./walled_x]
    type = WallEnergyDerivative
    variable = polar_x
    component = 0
    block = '1 2'
  [../]
  [./walled_y]
    type = WallEnergyDerivative
    variable = polar_y
    component = 1
    block = '1 2'
  [../]
  [./walled_z]
    type = WallEnergyDerivative
    variable = polar_z
    component = 2
    block = '1 2'
  [../]

  [./electrostr_ux]
    type = ElectrostrictiveCouplingDispDerivative
    variable = u_x
    component = 0
    block = '1 2'
  [../]
  [./electrostr_uy]
    type = ElectrostrictiveCouplingDispDerivative
    variable = u_y
    component = 1
    block = '1 2'
  [../]
  [./electrostr_uz]
    type = ElectrostrictiveCouplingDispDerivative
    variable = u_z
    component = 2
    block = '1 2'
  [../]

  [./electrostr_polar_coupled_x]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_x
    component = 0
    block = '1 2'
  [../]
  [./electrostr_polar_coupled_y]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_y
    component = 1
    block = '1 2'
  [../]
  [./electrostr_polar_coupled_z]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_z
    component = 2
    block = '1 2'
  [../]
  [./polar_x_electric_E]
    type = PolarElectricEStrong
    variable = potential_E_int
    block = '1 2'
  [../]
  [./FE_E_int]
    type = Electrostatics
    variable = potential_E_int
    block = '1 2 3'
  [../]

  [./polar_electric_px]
    type = PolarElectricPStrong
    variable = polar_x
    component = 0
    block = '1 2'
  [../]
  [./polar_electric_py]
    type = PolarElectricPStrong
    variable = polar_y
    component = 1
    block = '1 2'
  [../]
  [./polar_electric_pz]
    type = PolarElectricPStrong
    variable = polar_z
    component = 2
    block = '1 2'
  [../]

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

  [./u_x_time]
    type = TimeDerivativeScaled
    variable=u_x
    time_scale = 1.0
  [../]
  [./u_y_time]
    type = TimeDerivativeScaled
    variable=u_y
    time_scale = 1.0
  [../]
  [./u_z_time]
    type = TimeDerivativeScaled
    variable = u_z
    time_scale = 1.0
  [../]
[]


[BCs]
  [./Periodic]
    [./xy1]
      auto_direction = 'x y z'
      variable = 'u_x u_y u_z polar_x polar_y polar_z potential_E_int'
    [../]
    [./xy2]
      auto_direction = 'x y z'
      variable = 'u_x u_y u_z potential_E_int'
    [../]
  [../]


[]

[Postprocessors]

  ###############################################
  ##=
  ##  Postprocessors (integrations over the
  ##  computational domain) to calculate the total energy
  ##  decomposed into linear combinations of the
  ##  different physics.
  ##
  ###############################################

  [./Fbulk]
    type = BulkEnergyEighth
    execute_on = 'timestep_end'
    block = '1 2'
  [../]
  [./Fwall]
    type = WallEnergy
    execute_on = 'timestep_end'
    block = '1 2'
  [../]
  [./Felastic]
    type = ElasticEnergy
    execute_on = 'timestep_end'
    use_displaced_mesh = false
    block = '1 2'
  [../]
  [./Fcoupled]
    type = ElectrostrictiveCouplingEnergy
    execute_on = 'timestep_end'
    block = '1 2'
  [../]
  [./Felec]
    type = ElectrostaticEnergy
    execute_on = 'timestep_end'
    block = '1 2'
  [../]
  [./Ftotal]
    type = LinearCombinationPostprocessor
    pp_names = 'Fbulk Fwall Fcoupled Felec'
    pp_coefs = ' 1 1 1 1'
    execute_on = 'timestep_end'
  [../]

  [./perc_change]
    type = ChangeOverTimePostprocessor
    postprocessor = Ftotal
    execute_on = 'timestep_end'
    compute_relative_change = true
  [../]

[]

[UserObjects]
  ###############################################
  ##
  ##  terminator to end energy evolution when the energy difference
  ##  between subsequent time steps is lower than 1e-5
  ##
  ##  NOTE: can fail if the time step is small
  ##
  ###############################################

  [./kill]
    type = Terminator
    expression = 'perc_change <= 1.0e-4'
   [../]
[]

[Preconditioning]

  ###############################################
  ##
  ##  Numerical preconditioning/solver options
  ##=
  ###############################################

  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type  -build_twosided'
    petsc_options_value = '    160             1e-8        1e-6      1e-5        bjacobi       allreduce'
  [../]
[]

[Executioner]

  ##########################################
  ##
  ##  Time integration/solver options
  ##
  ##########################################

  type = Transient
  solve_type = 'NEWTON'
  scheme = 'implicit-euler'
  dtmin = 1e-13
  dtmax = 1.5

  l_max_its = 200

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    growth_factor = 1.2
    cutback_factor = 0.75
    linear_iteration_ratio = 1000
    dt = 0.6
  [../]
  verbose = true
  nl_max_its = 20
  num_steps = 1
[]

[Outputs]
  print_linear_residuals = false
  perf_graph = true
  [./out]
    type = Exodus
    file_base = out_bcc_pbc
    elemental_as_nodal = true
    interval = 1
  [../]
[]
