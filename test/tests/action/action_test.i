a1temp = -0.172197

[Mesh]
  [gen]
    ############################################
    ##
    ##  Type and dimension of the mesh 
    ##
    ############################################

    type = GeneratedMeshGenerator
    dim = 3

    #############################################
    ##
    ##  Grid definition. Note that it should be
    ##  nJ = 2*(Jmax-Jmin) for J = x, y, z
    ##
    #############################################

    nx = 5
    ny = 5
    nz = 5

    #############################################
    ##
    ##   Actual spatial coordinates of mesh. 
    ##   Jmax - Jmin = nJ/2 for J = x, y, z
    ##   Units are in nanometers
    ##
    #############################################

    xmin = -3.0
    xmax = 3.0
    ymin = -3.0
    ymax = 3.0
    zmin = -3.0
    zmax = 3.0

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
    coord = '-3.0 -3.0 -3.0'
    new_boundary = 100
  [../]
[]

[GlobalParams]
  len_scale = 1.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_E_int = potential_E_int

  ##############################################
  ##
  ##  IMPORTANT(!): Units in Ferret are nm, kg,
  ##                seconds, and attocoulombs
  ##
  ##############################################

  displacements = 'u_x u_y u_z'

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
      min = -0.1e-2
      max = 0.1e-2
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-2
      max = 0.1e-2
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-2
      max = 0.1e-2
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

  [./ez]
    type = QuasistaticFieldAux
    variable = Ez
    potential_int = potential_E_int
    component = 2
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

  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #C11 C12 C13 C22 C23 C33 C44 C55 C66 
    C_ijkl = '175.0 79.4 79.4 175.0 79.4 175.0 111.1 111.1 111.1'
  [../]

  [./strain_1]
    type = ComputeSmallStrain
    global_strain = global_strain
  [../]

  [./stress_1]
    type = ComputeLinearElasticStress
  [../]

  [./global_strain]
    type = ComputeGlobalStrain
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
  [../]
[]

[Ferret]

  ###############################################
  ##
  ## MOOSE and Ferret Kernel, Materials, and
  ## Postprocessors for TDLGD evolution 
  ##
  ###############################################

  [./ABO3CoupledPhaseField]

    #  {Variables to solve for}

    variables = 'polar_x polar_y polar_z potential_E_int u_x u_y u_z'

    ###############################################
    #
    #  {Governing Equation Flags}
    #
    ###############################################

    coupled_problem = true
    polar_time_dependence = true
    u_time_dependence = false
    phi_time_dependence = false


    ###############################################
    #
    #  {Materials: PbTiO3}
    #
    ###############################################

    alpha_ijkl = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    alpha_ijkl_val = '${a1temp} -0.073 0.75 0.26 0.61 -3.67 0.0 0.0 0.0 0.0'

    G_ij = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    G_ij_val = '0.173 0.6 0.0 0.3 0.3'

    is_permittivity_anisotropic = false
    permittivity = 'permittivity'
    permittivity_val = '0.08854187'

    Q_ij = 'Q11 Q12 Q44'
    Q_ij_val = '-0.089 0.026 -0.03375'

    q_ij = 'q11 q12 q44'
    q_ij_val = '-11.4 -0.01438 -7.5'
  
    C_ij = 'C11 C12 C44'
    C_ij_val = '175.0 79.4 111.1'
   [../]
[]



[Kernels]
  [./SolidMechanics]
  [../]
[]

[Functions]
  [./bc_func_1]
    type = ParsedFunction
    value = -5.0*sin(0.005*t)
  [../]
[]


[BCs]
  [./Periodic]
    [./xyz]
      auto_direction = 'x y z'
      variable = 'polar_x polar_y polar_z u_x u_y u_z '
    [../]
    [./xy]
      auto_direction = 'x y'
      variable = 'potential_E_int'
    [../]
  [../]


  [./boundary_grounding1]
    type = FunctionDirichletBC
    boundary = 'front'
    variable = potential_E_int
    function = bc_func_1
  [../]

  [./boundary_grounding2]
    type = DirichletBC
    boundary = 'back'
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

  ###############################################
  ##
  ##  Postprocessors (integrations over the 
  ##  computational domain) to calculate the total energy
  ##  decomposed into linear combinations of the
  ##  different physics.
  ##
  ###############################################

  [./avePx]
    type = ElementAverageValue
    variable = polar_x
    execute_on = 'timestep_end'
  [../]
  [./avePy]
    type = ElementAverageValue
    variable = polar_y
    execute_on = 'timestep_end'
  [../]
  [./avePz]
    type = ElementAverageValue
    variable = polar_z
    execute_on = 'timestep_end'
  [../]

  [./aveEz]
    type = ElementAverageValue
    variable = Ez
    execute_on = 'timestep_end'
  [../]


  [./ave_e00]
    type = ElementAverageValue
    variable = e00
    execute_on = 'timestep_end'
  [../]
  [./ave_e11]
    type = ElementAverageValue
    variable = e11
    execute_on = 'timestep_end'
  [../]
  [./ave_e22]
    type = ElementAverageValue
    variable = e22
    execute_on = 'timestep_end'
  [../]

  [./Fb]
    type = BulkEnergyEighth
    execute_on = 'timestep_end'
  [../]
  [./Fw]
    type = WallEnergy
    execute_on = 'timestep_end'
  [../]
  [./Fela]
    type = ElasticEnergy
    execute_on = 'timestep_end'
    use_displaced_mesh = false
  [../]
  [./Fc]
    type = ElectrostrictiveCouplingEnergy
    execute_on = 'timestep_end'
  [../]
  [./Fele]
    type = ElectrostaticEnergy
    execute_on = 'timestep_end'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'Fb Fw Fc Fele'
    pp_coefs = '0.160218 0.160218 0.160218 0.160218'
    execute_on = 'timestep_end'
  [../]
  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = Ftot
    execute_on = 'timestep_end'
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
  [../]

  ###############################################
  ##
  ##  terminator to end energy evolution when the energy difference 
  ##  between subsequent time steps is lower than 5e-6
  ##
  ##  NOTE: can fail if the time step is small
  ##
  ###############################################

  #[./kill]
  #  type = Terminator
  #  expression = 'perc_change <= 1.0e-4'
  # [../]
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
    petsc_options_value = '    160               1e-8      1e-6      1e-6          bjacobi       allreduce'
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
  ##
  ###########################################

  dtmax = 1.0

  l_max_its = 200

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 8
    cutback_factor = 0.75
    linear_iteration_ratio = 1000
    dt = 0.3
  [../]
  verbose = true
  num_steps = 2
[]



[Outputs]

  ###############################################
  ##
  ##  Output options
  ##
  ###############################################

  print_linear_residuals = false
  perf_graph = false

  [./outExo]
    type = Exodus
    file_base = out_pto-md_action
  [../]

  [./outCSV]
    type = CSV
    file_base = out_pto-md_action
  [../]
[]
