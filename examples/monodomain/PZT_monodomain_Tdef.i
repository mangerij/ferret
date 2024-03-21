a1temp = a1def

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

    nx = 8
    ny = 8
    nz = 8

    #############################################
    ##
    ##   Actual spatial coordinates of mesh. 
    ##   Jmax - Jmin = nJ/2 for J = x, y, z
    ##   Units are in nanometers
    ##
    #############################################

    xmin = -2.0
    xmax = 2.0
    ymin = -2.0
    ymax = 2.0
    zmin = -2.0
    zmax = 2.0

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
    coord = '0.0 0.0 0.0'
    new_boundary = 100
  [../]
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
      min = -0.1e-6
      max = 0.1e-6
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-6
      max = 0.1e-6
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.05
      max = 0.1
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


  [./eigs00]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eigs11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eigs22]
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


  [./eigs00]
    type = LocalABO3EigenstressAux
    variable = eigs00
    index_i = 0
    index_j = 0
  [../]
  [./eigs11]
    type = LocalABO3EigenstressAux
    variable = eigs11
    index_i = 1
    index_j = 1
  [../]
  [./eigs22]
    type = LocalABO3EigenstressAux
    variable = eigs22
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
  ## Landau, electrostrictive, elastic coefficients  
  ##  M. Mtebwa, A. K. Tagantsev, and N. Setter
  ##    AIP Adv. 4, 127150 (2014).
  ##
  ##################################################

  [./Landau_P_FE]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '${a1temp} 0.04764 0.1336 0.1735 0.6128 -2.894 0.0 0.0 0.0 0.0'
    block = '0'
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
    prop_values = '179.073 66.71 82.6446'
    block = '0'
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
    prop_values = '-0.0966 0.046 -0.0819'
    block = '0'
  [../]

  [./mat_q]
    type = GenericConstantMaterial
    prop_names = 'q11 q12 q44'
    prop_values = '-11.1608 4.86165 -6.76859'
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

    C_ijkl = '179.073 66.71 66.71 179.073 66.71 179.073 82.6446 82.6446 82.6446'
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
    eigenstrains_name = eigenstrain
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
  [./Periodic]
    [./xyz]
      auto_direction = 'x y z'
      variable = 'u_x u_y u_z polar_x polar_y polar_z'
    [../]
  [../]

  [./boundary_grounding]
    type = DirichletBC
    boundary = '0 1 2 3 4 5'
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
    pp_coefs = ' 1 1 1 1'
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
  [../]
  [./py]
    type = DomainVariantPopulation
    execute_on = 'timestep_end'
    component = 1
  [../]
  [./pz]
    type = DomainVariantPopulation
    execute_on = 'timestep_end'
    component = 2
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
  ##  NOTE: can fail if the time step is smallhotkey for tilde
  ##
  ###############################################

  [./kill]
    type = Terminator
    expression = 'perc_change <= 1.0e-6'
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
[]



[Outputs]

  ###############################################
  ##
  ##  Output options
  ##
  ###############################################

  print_linear_residuals = false
  perf_graph = false

  [./outCSV]
    type = CSV
    file_base = out_pzt_monodomain_Tdef
  [../]
[]
