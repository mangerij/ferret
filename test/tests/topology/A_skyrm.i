[Mesh]
  file = exodus_disk_r8_h1.e
[]

[GlobalParams]
  len_scale = 1.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_E_int = potential_E_int
  vol = vol

  displacements = 'u_x u_y u_z'


  ##############################################
  ##
  ##  IMPORTANT(!): Units in Ferret are nm, kg,
  ##                seconds, and attocoulombs
  ##
  ##############################################

[]

[Functions]
  [./parsed_function_x_skyrm]
    type = ParsedFunction
    value = '-(0.738217-0.00686984*(x^2+y^2)^(0.5)+0.00644497*(x^2+y^2)-0.0188174*(x^2+y^2)^(1.5)+0.00441745*(x^2+y^2)^2-0.000274842*(x^2+y^2)^(5/2))*sin(-0.028395+0.267482*(x^2+y^2)^(0.5)-0.146762*(x^2+y^2)+0.0632932*(x^2+y^2)^(1.5)-0.00790942*(x^2+y^2)^(2)+0.000294936*(x^2+y^2)^(5/2))*sin(atan(y/x))'
  [../]
  [./parsed_function_y_skyrm]
    type = ParsedFunction
    value = '(0.738217-0.00686984*(x^2+y^2)^(0.5)+0.00644497*(x^2+y^2)-0.0188174*(x^2+y^2)^(1.5)+0.00441745*(x^2+y^2)^2-0.000274842*(x^2+y^2)^(5/2))*sin(-0.028395+0.267482*(x^2+y^2)^(0.5)-0.146762*(x^2+y^2)+0.0632932*(x^2+y^2)^(1.5)-0.00790942*(x^2+y^2)^(2)+0.000294936*(x^2+y^2)^(5/2))*cos(atan(y/x))'
  [../]
  [./parsed_function_z_skyrm]
    type = ParsedFunction
    value = '(0.738217-0.00686984*(x^2+y^2)^(0.5)+0.00644497*(x^2+y^2)-0.0188174*(x^2+y^2)^(1.5)+0.00441745*(x^2+y^2)^2-0.000274842*(x^2+y^2)^(5/2))*cos(-0.028395+0.267482*(x^2+y^2)^(0.5)-0.146762*(x^2+y^2)+0.0632932*(x^2+y^2)^(1.5)-0.00790942*(x^2+y^2)^(2)+0.000294936*(x^2+y^2)^(5/2))'
  [../]

[]

[Variables]

  #################################
  ##
  ##  Variable definitions
  ##    P, u, phi, e^global_ij  
  ##  and their initial conditions
  ##
  #################################

  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = parsed_function_x_skyrm
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = parsed_function_y_skyrm
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = parsed_function_z_skyrm
    [../]
  [../]

  [./potential_E_int]
    order = FIRST
    family = LAGRANGE
  [../]

  [./u_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./u_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./u_z]
    order = FIRST
    family = LAGRANGE
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
  [../]
  [./e01]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e22]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./s00]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s01]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s22]
    order = CONSTANT
    family = MONOMIAL
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
  ## Bulk free energy and electrostrictive
  ## coefficients gleaned from
  ## Marton and Hlinka
  ##    Phys. Rev. B. 74, 104014, (2006)
  ##
  ## NOTE: there might be some Legendre transforms
  ##        depending on what approach you use
  ##        -i.e. inhomogeneous strain vs 
  ##            homogeneous strain [renormalized]
  ##
  ##################################################

  [./Landau_P]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '-0.027721 -0.64755 0.323 8.004 4.47 4.91 0.0 0.0 0.0 0.0'
  [../]

  ############################################
  ##
  ## Gradient coefficients from
  ## Marton and Hlinka
  ##    Phys. Rev. B. 74, 104014, (2006)
  ##
  ############################################

  [./Landau_G]
    type = GenericConstantMaterial
    prop_names = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '0.5 0.51 -0.02 0.02 0.0'
  [../]

  [./mat_C]
    type = GenericConstantMaterial
    prop_names = 'C11 C12 C44'
    prop_values = '275.0 179.0 54.3'
  [../]

  ##################################################
  ##
  ## NOTE: Sign convention in Ferret for the 
  ##        electrostrictive coeff. is multiplied by
  ##        an overall factor of (-1)
  ##
  ##################################################

  [./mat_Q]
    type = GenericConstantMaterial
    prop_names = 'Q11 Q12 Q44'
    prop_values = '-0.110448 0.0451767 -0.115654'
  [../]
  [./mat_q]
    type = GenericConstantMaterial
    prop_names = 'q11 q12 q44'
    prop_values = '-14.2 0.74 -6.28'
  [../]

  [./eigen_strain]
    type = ComputeEigenstrain
    eigen_base = '0. 0 0 0 0 0 0 0 0'
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

    C_ijkl = '275.0 179.0 179.0 275.0 179.0 275.0 54.3 54.3 54.3'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    eigenstrain_names = eigenstrain
  [../]

  [./stress_1]
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
  [./SolidMechanics]
    use_displaced_mesh = false
  [../]

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

  [./electrostr_ux]
    type = ElectrostrictiveCouplingDispDerivative
    variable = u_x
    component = 0

  [../]
  [./electrostr_uy]
    type = ElectrostrictiveCouplingDispDerivative
    variable = u_y
    component = 1
  [../]
  [./electrostr_uz]
    type = ElectrostrictiveCouplingDispDerivative
    variable = u_z
    component = 2
  [../]

  [./electrostr_polar_coupled_x]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_x
    component = 0
  u_x = u_x
  u_y = u_y
  u_z = u_z
  [../]
  [./electrostr_polar_coupled_y]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_y
    component = 1
  u_x = u_x
  u_y = u_y
  u_z = u_z
  [../]
  [./electrostr_polar_coupled_z]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_z
    component = 2
  u_x = u_x
  u_y = u_y
  u_z = u_z
  [../]


  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_E_int
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_E_int
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
     variable = u_x
     time_scale = 1.0
  [../]
  [./u_y_time]
     type = TimeDerivativeScaled
     variable = u_y
     time_scale = 1.0
  [../]
  [./u_z_time]
     type = TimeDerivativeScaled
     variable = u_z
     time_scale = 1.0
  [../]
[]


[BCs]
   [./top_electrode_top]
      type = DirichletBC
      variable = 'potential_E_int'
      value = 0.0001
      boundary = '2'
   [../]

   [./top_electrode_bottom]
      type = DirichletBC
      variable = 'potential_E_int'
      value = 0.0001
      boundary = '3'
   [../]
[]

[Postprocessors]

  ###############################################
  ##
  ##  Postprocessors (integrations over the 
  ##  computational domain) to calculate the total energy
  ##  decomposed into linear combinations of the
  ##  different physics.
  ##
  ###############################################

  [./Fb]
    type = BulkEnergyEighth
    execute_on = 'initial timestep_end'
  [../]
  [./Fw]
    type = WallEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./Fela]
    type = ElasticEnergy
    execute_on = 'initial timestep_end'
    use_displaced_mesh = false
  [../]
  [./Fc]
    type = ElectrostrictiveCouplingEnergy
    execute_on = 'initial timestep_end'
  u_x = u_x
  u_y = u_y
  u_z = u_z
  [../]
  [./Fele]
    type = ElectrostaticEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'Fb Fw Fc Fele'
    pp_coefs = ' 1 1 1 1'
    execute_on = 'initial timestep_end'
  [../]
  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = Ftot
    execute_on = 'initial timestep_end'
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
  ##
  ###############################################

  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type  -build_twosided'
    petsc_options_value = '    160               1e-10      1e-8      1e-6          bjacobi       allreduce'
  [../]
[]

[Executioner]

  ##########################################
  ##
  ##  Time integration/solver options
  ##
  ##########################################

  type = Transient
  solve_type = 'PJFNK'
  scheme = 'implicit-euler'
  dtmin = 1e-13

  ###########################################
  ##
  ##  dtmax is material dependent!
  ##  for PTO is about 0.8 but BTO more like 3-10
  ##
  ###########################################

  dtmax = 8.0

  l_max_its = 200

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 8
    cutback_factor = 0.75
    growth_factor = 1.25
    linear_iteration_ratio = 1000
    dt = 0.000008
  [../]
  verbose = true
  num_steps = 1
[]



[Outputs]

  ###############################################
  ##
  ##  Output options
  ##
  ###############################################

  print_linear_residuals = false
  perf_graph = false

  [./out]
    type = Exodus
    file_base = out_bto_skyrm_disk
    elemental_as_nodal = true
  [../]
[]
