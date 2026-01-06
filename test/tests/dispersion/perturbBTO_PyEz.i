freq = 177827941.00389227
amplitude = 1e-3

[Mesh]
   file = out_PyEz.e
[]

[GlobalParams]
  len_scale = 1.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_E_int = potential_E_int

  displacements = 'u_x u_y u_z'


  ##############################################
  ##
  ##  IMPORTANT(!): Units in Ferret are nm, kg,
  ##                seconds, and attocoulombs
  ##
  ##############################################

  u_x = u_x
  u_y = u_y
  u_z = u_z

[]


[Functions]
  ##############################
  ##
  ## Define the electric field
  ## expression to be used below
  ##
  ##############################

  [./bc_func_1]
    type = ParsedFunction
    expression = 'amplitude*sin(freq*t)'
    symbol_names = 'freq amplitude'
    symbol_values = '${freq}  ${amplitude}'
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

  [./global_strain]
    order = SIXTH
    family = SCALAR
  [../]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_x
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_y
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_z
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
  ##  Global displacements
  ##
  ######################################

  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]

  ######################################
  ##
  ##  Stress/strain tensor components
  ##
  ######################################

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

  [./Ez]
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

  [./disp_x]
    type = GlobalDisplacementAux
    variable = disp_x
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
    component = 0
    use_displaced_mesh = false
  [../]
  [./disp_y]
    type = GlobalDisplacementAux
    variable = disp_y
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
    component = 1
    use_displaced_mesh = false
  [../]
  [./disp_z]
    type = GlobalDisplacementAux
    variable = disp_z
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
    component = 2
    use_displaced_mesh = false
  [../]
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

  [./cEy]
    type = QuasistaticFieldAux
    component = 2
    potential_int = potential_E_int
    variable = Ez
  [../]

[]

[ScalarKernels]

  ######################################
  ##
  ##  Necessary for PBC system
  ##
  ######################################

  [./global_strain]
    type = GlobalStrain
    variable = global_strain
    global_strain_uo = global_strain_uo
    use_displaced_mesh = false
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

  ##############################################################
  ##
  ## NOTE: Sign convention in **this implementation**
  ##       for the electrostrictive coeff. is multiplied by
  ##       an overall factor of (-1). Note that other elastic
  ##       coupling Kernels/Materials in Ferret DO NOT have the 
  ##       (-1) prefactor. Please be careful here.
  ##
  ###############################################################

  [./mat_Q]
    type = GenericConstantMaterial
    prop_names = 'Q11 Q12 Q44'
    prop_values = '-0.11 0.045 -0.029'
  [../]
  [./mat_q]
    type = GenericConstantMaterial
    prop_names = 'q11 q12 q44'
    prop_values = '-14.2 0.74 -1.57'
  [../]

  [./eigen_strain]
    type = ComputeEigenstrain
    eigen_base = '0.0 0.0 0 0 0 0 0 0 0'
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
    global_strain = global_strain
    eigenstrain_names = eigenstrain
  [../]

  [./stress_1]
    type = ComputeLinearElasticStress
  [../]

  [./global_strain]
    type = ComputeGlobalStrain
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
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
  [../]
  [./electrostr_polar_coupled_y]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_y
    component = 1
  [../]
  [./electrostr_polar_coupled_z]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_z
    component = 2
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
     variable = polar_x
     # Time scale estimate for BTO, from Hlinka (2007)
     # We use seconds here
     time_scale = 1e-12
  [../]
  [./polar_y_time]
     type = TimeDerivativeScaled
     variable = polar_y
     time_scale = 1e-12
  [../]
  [./polar_z_time]
     type = TimeDerivativeScaled
     variable = polar_z
     time_scale = 1e-12
  [../]
[]


[BCs]
  [./Periodic]
    [./xyz]
      auto_direction = 'x y z'
      variable = 'u_x u_y u_z polar_x polar_y polar_z'
    [../]
  [../]


  [./front_pot]
    type = FunctionDirichletBC
    variable = potential_E_int
    boundary = 'front'
    function = bc_func_1
  [../]

  [./boundary_grounding]
    type = DirichletBC
    boundary = 'left right top bottom back'
    variable = potential_E_int
    value = 0.0
  [../]


  # fix center point location
  [./centerfix_x]
    type = DirichletBC
    boundary = 100
    variable = u_x
    value = 0
  [../]
  [./centerfix_y]
    type = DirichletBC
    boundary = 100
    variable = u_y
    value = 0
  [../]
  [./centerfix_z]
    type = DirichletBC
    boundary = 100
    variable = u_z
    value = 0
  [../]
[]

[Postprocessors]

  [./avePz]
    type = ElementAverageValue
    variable = polar_z
    execute_on = 'initial timestep_end'
  [../]
  [./Ea]
    type = ElementAverageValue
    variable = Ez
    execute_on = 'initial timestep_end'
  [../]


  ###############################################
  ##
  ##  Postprocessors (integrations over the
  ##  computational domain) to calculate the total energy
  ##  decomposed into linear combinations of the
  ##  different physics.
  ##
  ###############################################

  [./Fbulk]
    type = BulkEnergyEighth
    execute_on = 'initial timestep_end'
  [../]
  [./Fwall]
    type = WallEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./Felastic]
    type = ElasticEnergy
    execute_on = 'initial timestep_end'
    use_displaced_mesh = false
  [../]
  [./Fcoupled]
    type = ElectrostrictiveCouplingEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./Felec]
    type = ElectrostaticEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./Ftotal]
    type = LinearCombinationPostprocessor
    pp_names = 'Fbulk Fwall Fcoupled Felec'
    pp_coefs = ' 1 1 1 1'
    execute_on = 'initial timestep_end'
  [../]
  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = Ftotal
    execute_on = 'initial timestep_end'
  [../]
[]

[UserObjects]

  ###############################################
  ##
  ##  GlobalStrain system to enforce periodicity
  ##  in the anisotropic strain field
  ##
  ###############################################

  [./global_strain_uo]
    type = GlobalATiO3MaterialRVEUserObject
    use_displaced_mesh = false
    execute_on = 'Initial Linear Nonlinear'
    applied_stress_tensor = '0.0 0.0 0.0 0.0 0.0 0.0'
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

  dtmin = 1e-16

  dt = 3.5332947520558995e-10

  dtmax = 1e-10

  verbose = true

  num_steps = 5
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
    file_base = out_perturbBTO_PyEz0
    elemental_as_nodal = true
  [../]
  [./outCSV]
    type = CSV
    new_row_tolerance = 1e-16
    file_base = out_perturbBTO_PyEz0
  [../]
[]
