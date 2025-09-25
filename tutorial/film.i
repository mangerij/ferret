[Mesh]
  [gen]
    ############################################
    ##
    ##  Type and dimension of the mesh
    ##
    ############################################

    type = GeneratedMeshGenerator
    dim = 3

    nx = 32
    ny = 32
    nz = 30

    xmin = -16.0
    xmax = 16.0
    ymin = -16.0
    ymax = 16.0
    zmin = -10.0
    zmax = 20.0

    #############################################
    ##
    ##  FE type/order (hexahedral, tetrahedral
    ##
    #############################################

    elem_type = HEX8
  []
  [./cnode]
    input = gen

    ############################################
    ##
    ##   additional boundary sideset (one node)
    ##   to zero one of the elastic displacement vectors
    ##   vectors and eliminates rigid body translations
    ##   from the degrees of freedom
    ##
    ##   NOTE: This must conform with the about
    ##         [Mesh] block settings
    ##
    ############################################

    type = ExtraNodesetGenerator
    coord = '-16.0 -16.0 -10.0'
    new_boundary = 100
  [../]

  [subdomains]
    type = SubdomainBoundingBoxGenerator
    input = cnode
    bottom_left = '-16.0 -16.0 -10.0'
    block_id = 1
    top_right = '16.0 16.0 0'
    location = INSIDE
  []
  [film_interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = subdomains
    primary_block = 0
    paired_block = 1
    new_boundary = 52
  []
  [film_surface]
    type = SideSetsFromNormalsGenerator
    input = film_interface
    normals = '0  0  1'
    fixed_normal = true
    new_boundary = '107'
  []
  [substrate_bottom]
    type = SideSetsFromNormalsGenerator
    input = film_surface
    normals = '0  0  -1'
    fixed_normal = true
    new_boundary = '108'
  []
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

  vol = vol

  u_x = u_x
  u_y = u_y
  u_z = u_z
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
    [./InitialCondition]
      type = RandomIC
      min = -1e-2
      max = 1e-2
    [../]
    block = '0'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -1e-2
      max = 1e-2
    [../]
    block = '0'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -1e-2
      max = 1e-2
    [../]
    block = '0'
  [../]

  [./potential_E_int]
    order = FIRST
    family = LAGRANGE
    block = '0 1'
  [../]

  [./u_x]
    order = FIRST
    family = LAGRANGE
    block = '0 1'
  [../]
  [./u_y]
    order = FIRST
    family = LAGRANGE
    block = '0 1'
  [../]
  [./u_z]
    order = FIRST
    family = LAGRANGE
    block = '0 1'
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
    block = '0'
  [../]
  [./bed_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
    block = '0'
  [../]
  [./bed_z]
    type = BulkEnergyDerivativeEighth
    variable = polar_z
    component = 2
    block = '0'
  [../]

  [./walled_x]
    type = WallEnergyDerivative
    variable = polar_x
    component = 0
    block = '0'
  [../]
  [./walled_y]
    type = WallEnergyDerivative
    variable = polar_y
    component = 1
    block = '0'
  [../]
  [./walled_z]
    type = WallEnergyDerivative
    variable = polar_z
    component = 2
    block = '0'
  [../]

  [./walled2_x]
    type = Wall2EnergyDerivative
    variable = polar_x
    component = 0
    block = '0'
  [../]
  [./walled2_y]
    type = Wall2EnergyDerivative
    variable = polar_y
    component = 1
    block = '0'
  [../]
  [./walled2_z]
    type = Wall2EnergyDerivative
    variable = polar_z
    component = 2
    block = '0'
  [../]

  [./electrostr_ux]
    type = ElectrostrictiveCouplingDispDerivative
    variable = u_x
    component = 0
    block = '0'
  [../]
  [./electrostr_uy]
    type = ElectrostrictiveCouplingDispDerivative
    variable = u_y
    component = 1
    block = '0'
  [../]
  [./electrostr_uz]
    type = ElectrostrictiveCouplingDispDerivative
    variable = u_z
    component = 2
    block = '0'
  [../]

  [./electrostr_polar_coupled_x]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_x
    component = 0
    block = '0'
  [../]
  [./electrostr_polar_coupled_y]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_y
    component = 1
    block = '0'
  [../]
  [./electrostr_polar_coupled_z]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_z
    component = 2
    block = '0'
  [../]


  [./polar_x_electric_E]
    type = PolarElectricEStrong
    variable = potential_E_int
    block = '0'
  [../]
  [./FE_E_int]
    type = Electrostatics
    variable = potential_E_int
    block = '0 1'
  [../]

  [./polar_electric_px]
    type = PolarElectricPStrong
    variable = polar_x
    component = 0
    block = '0'
  [../]
  [./polar_electric_py]
    type = PolarElectricPStrong
    variable = polar_y
    component = 1
    block = '0'
  [../]
  [./polar_electric_pz]
    type = PolarElectricPStrong
    variable = polar_z
    component = 2
    block = '0'
  [../]

  [./polar_x_time]
    type = TimeDerivativeScaled
    variable=polar_x
    time_scale = 1.0
    block = '0'
  [../]
  [./polar_y_time]
    type = TimeDerivativeScaled
    variable = polar_y
    time_scale = 1.0
    block = '0'
  [../]
  [./polar_z_time]
    type = TimeDerivativeScaled
    variable = polar_z
    time_scale = 1.0
    block = '0'
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
    block = '0 1'
  [../]
  [./disp_y]
    block = '0 1'
  [../]
  [./disp_z]
    block = '0 1'
  [../]

######################################
  ##
  ##  Stress/strain tensor components
  ##
  ######################################

  [./e00]
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
  [../]
  [./e01]
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
  [../]
  [./e10]
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
  [../]
  [./e11]
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
  [../]
  [./e12]
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
  [../]
  [./e22]
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
  [../]

  [./s00]
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
  [../]
  [./s01]
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
  [../]
  [./s10]
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
  [../]
  [./s11]
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
  [../]
  [./s12]
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
  [../]
  [./s22]
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
  [../]

  [./E_x]
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
  [../]
  [./E_y]
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
  [../]
  [./E_z]
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
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
  ## Landau coefficients from Li et al (2001)
  ##
  ##################################################

  [./Landau_P_FE]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '-0.1722883 -0.073 0.75 0.26 0.61 -3.67 0.0 0.0 0.0 0.0'
    block = '0'
  [../]

  [./Landau_P_substr]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '10.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    block = '1'
  [../]

  [./Landau_G_FE]
    type = GenericConstantMaterial
    prop_names = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '0.173 0.6 0.0 0.3 0.3'
    block = '0'
  [../]

  [./mat_C_FE]
    type = GenericConstantMaterial
    prop_names = 'C11 C12 C44'
    prop_values = '175.0 79.4 111.1'
    block = '0'
  [../]
  [./mat_C_sub]
    type = GenericConstantMaterial
    prop_names = 'C11 C12 C44'
    prop_values = '220.0 34.4 161.1'
    block = '1'
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
    prop_values = '-0.089 0.026 -0.03375'
    block = '0 1'
  [../]

  [./mat_q]
    type = GenericConstantMaterial
    prop_names = 'q11 q12 q44'
    prop_values = '-11.4 0.01438 -7.5'
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


[BCs]
  [./Periodic]
    [./xy]
      auto_direction = 'x y'
      variable = 'u_x u_y u_z polar_x polar_y polar_z potential_E_int'
    [../]
  [../]



  [./boundary_interface_grounding]
    type = DirichletBC
    boundary = '52'
    variable = potential_E_int
    value = 0.0
  [../]


  # fix center point location
  [./centerfix_x]
    type = DirichletBC
    boundary = '108'
    variable = u_x
    value = 0
  [../]
  [./centerfix_y]
    type = DirichletBC
    boundary = '108'
    variable = u_y
    value = 0
  [../]
  [./centerfix_z]
    type = DirichletBC
    boundary = '108'
    variable = u_z
    value = 0
  [../]
[]

[Postprocessors]

  ###############################################
  ##=
  ##  Postprocessors (integrations over the
  ##  computational domain) to calculate the total
  ##  energy decomposed into linear combinations of
  ##  the different physics.
  ##
  ###############################################

  [./Fbulk]
    type = BulkEnergyEighth
    execute_on = 'timestep_end'
    block = '0'
  [../]
  [./Fwall]
    type = WallEnergy
    execute_on = 'timestep_end'
    block = '0'
  [../]
  [./Felastic]
    type = ElasticEnergy
    execute_on = 'timestep_end'
    use_displaced_mesh = false
    block = '0'
  [../]
  [./Fcoupled]
    type = ElectrostrictiveCouplingEnergy
    execute_on = 'timestep_end'
    block = '0'
  [../]
  [./Felec]
    type = ElectrostaticEnergy
    execute_on = 'timestep_end'
    block = '0'
  [../]
  [./Ftotal]
    type = LinearCombinationPostprocessor
    pp_names = 'Fbulk Fwall Fcoupled Felec'
    pp_coefs = '0.160218 0.160218 0.160218 0.160218' #converted to eV
    execute_on = 'timestep_end'
  [../]
  [./vol]
    type = VolumePostprocessor
    execute_on = 'timestep_end'
  [../]

  [./px]
    type = DomainVariantPopulation
    execute_on = 'timestep_end'
    component = 0
    block = '0'
  [../]
  [./py]
    type = DomainVariantPopulation
    execute_on = 'timestep_end'
    component = 1
    block = '0'
  [../]
  [./pz]
    type = DomainVariantPopulation
    execute_on = 'timestep_end'
    component = 2
    block = '0'
  [../]


  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = Ftotal
    execute_on = 'timestep_end'
  [../]
  [./elapsed]
    type = PerfGraphData
    section_name = "Root"  # for profiling the problem [on]
    data_type = total
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
    applied_stress_tensor = '2.1 2.1 1.9056 0.0 0.0 0.0'
    block = '0'
  [../]

  ###############################################
  ##
  ##  terminator to end energy evolution when the energy difference
  ##  between subsequent time steps is lower than 5e-6
  ##
  ##  NOTE: can fail if the time step is small
  ##
  ###############################################

  [./kill]
    type = Terminator
    expression = 'perc_change <= 5.0e-4'
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
    petsc_options_value = '    80             1e-8        1e-6      1e-5       bjacobi      allreduce'
  [../]
[]

[Executioner]

  ##########################################
  ##
  ##  Time integ=ration/solver options
  ##
  ##########################################

  type = Transient
  solve_type = 'PJFNK'
  scheme = 'bdf2'
  dtmin = 1e-13
  dtmax = 0.6

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
[]

[Outputs]

  ###############################################
  ##==
  ##  Output options
  ##
  ###############################################

  print_linear_residuals = false
  perf_graph = false

  [./out]
    type = Exodus
    file_base = out_PTOfilm_e12_T298K_E0_E0
    elemental_as_nodal = true
    time_step_interval = 1
  [../]
[]
