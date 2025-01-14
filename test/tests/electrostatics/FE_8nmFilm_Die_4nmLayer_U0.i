
[Mesh]
  [gen]
    #######################################################
    ##
    ##  Type and dimension of the mesh
    ##
    #######################################################

    type = GeneratedMeshGenerator
    dim = 2

    #######################################################
    ##
    ##  Finite element grid definition. For ferroelectric
    ##  calculations, this should be limited by the
    ##  domain wall width. We suggest ~1.2-2 elements per
    ##  1 nm of spatial resolution)
    ##
    #######################################################

    nx = 150
    ny = 20

    #######################################################
    ##
    ##   Actual spatial coordinates of mesh.
    ##   Units are in nanometers. Note the above constraint
    ##
    #######################################################

    xmin = -30.0
    xmax = 30.0
    ymin = -4.0
    ymax = 4.0

    #######################################################
    ##
    ##  Finite element type/order (hexahedral, tetrahedral)
    ##
    #######################################################

    elem_type = QUAD4
  []

  [subdomains1]
    #######################################################
    ##
    ##  One of the many mesh modifiers in MOOSE
    ##  This one splits the computational domain into
    ##  two (and then three) different layers
    ##  These layers are where we change the physics to
    ##  reflect the different materials
    ##
    #######################################################

    type = SubdomainBoundingBoxGenerator
    input = gen
    bottom_left = '-30.0 3.0 0.0'
    block_id = 1
    top_right = '30.0 4.0 0.0'
  []
  [subdomains2]
    type = SubdomainBoundingBoxGenerator
    input = subdomains1
    bottom_left = '-30.0 -4.0 0.0'
    block_id = 2
    top_right = '30.0 -3.0 0.0'
  []
[]

[GlobalParams]

  ##############################################
  ##
  ##  IMPORTANT(!): Units in Ferret are nm, kg,
  ##                seconds, and attocoulombs
  ##
  ##############################################

  len_scale = 1.0

  polar_x = polar_x
  polar_y = polar_y

  potential_E_int = potential_E_int
[]

[Variables]

  #################################
  ##
  ##  Variable definitions
  ##    Px, Py, \Phi
  ##  and their initial conditions
  ##
  #################################

  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block = '0 1 2'
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-6
      max = 0.01e-6
      seed = 6
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block = '0 1 2'
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-6
      max = 0.01e-6
      seed = 6
    [../]
  [../]

  [./potential_E_int]
    order = FIRST
    family = LAGRANGE
    block = '0 1 2'
  [../]
[]

[Materials]

  #################################################
  ##
  ## Bulk free energy and gradient
  ## coefficients retrieved from
  ## Park and Hwang
  ##    Adv. Mater. 2019, 31, 1805266
  ##
  ##################################################

  [./Landau_P_FE]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '-0.027722 -0.6381 0.0 7.89 0.0 0.0 0.0 0.0 0.0 0.0'
    block = '0'
  [../]

  [./Landau_P_die_layer]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '8.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0'   # arbitrary positive value set
    block = '1 2'
  [../]

  [./Landau_G]
    type = GenericConstantMaterial
    prop_names = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '0.5 0.0138 0.0 0.0138 0.0'   #note the 0.5 out front is set with G110 instead of in the expression.
  [../]

  [./in_P_susc]
    type = GenericConstantMaterial
    prop_names = 'chi'
    prop_values = '10.0'   # arbitrary positive value set
  [../]

  ###############################################
  ##
  ##  Background dielectric constant of the
  ##  dielectric layer. It is isotropic e_d = 208
  ##  where e_0 = 0.008854187 in Ferret's units
  ##
  ##############################################

  [./permitivitty]
    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '1.84167'
    block = '1 2'
  [../]

  [./eps1]
    ###############################################
    ##
    ##  background dielectric constant of the
    ##  ferroelectric layer is anisotropic with two
    ##  components. Note that e_0 is 0.008854187 in
    ##  Ferret's native units
    ##
    ##############################################

    type = GenericConstantMaterial
    prop_names = 'eps1 eps2 eps3'
    prop_values = '3.692195979 1.5450556315 0.0'
    block = '0'
  [../]
[]


[Kernels]

  ###############################################
  ##
  ## Physical Kernel operators
  ## to enforce TDLGD evolution
  ##
  ###############################################

  [./bed_FE_x]
    type = BulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
    block = '0'
  [../]

  [./bed_FE_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
    block = '0'
  [../]

  [./ip_FE_x]
    type = InPlaneSusceptibilityDerivative
    variable = polar_x
    block = '0'
  [../]


  [./bed_die_x]
    type = BulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
    block = '1 2'
  [../]

  [./bed_die_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
    block = '1 2'
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

  [./polar_x_electric_E]
    type = PolarElectricEStrong
    variable = potential_E_int
    block = '0 1 2'
  [../]

  [./FE_E_int]
    type = AnisotropicElectrostatics
    variable = potential_E_int
    block = '0'
  [../]

  [./die_E_int]
    type = Electrostatics
    variable = potential_E_int
    block = '1 2'
  [../]

  [./polar_electric_px]
    type = PolarElectricPStrong
    variable = polar_x
    component = 0
    block = '0 1 2'
  [../]
  [./polar_electric_py]
    type = PolarElectricPStrong
    variable = polar_y
    component = 1
    block = '0 1 2'
  [../]


  [./polar_x_time]
    type = TimeDerivativeScaled
    variable = polar_x

    ###############################################
    ##
    ##  Comments on time_scale:
    ##      time_scale is equal to 1/Gamma
    ##      currently it is set to 1.0
    ##      which means the time scale of the
    ##      problem is "quasi-static". It can be
    ##      turned into real units easily.
    ##
    ##############################################

    time_scale = 1.0
  [../]
  [./polar_y_time]
    type = TimeDerivativeScaled
    variable=polar_y
    time_scale = 1.0
  [../]
[]


[BCs]
  ####################################################
  ##
  ##  Simple BC: periodicity along x direction in
  ##  the polar field and Dirichlet conditions on the
  ##  electrostatic potential at the electrodes
  ##
  ####################################################

  [./Periodic]
    [./xyz]
      auto_direction = 'x'
      variable = 'polar_x polar_y'
    [../]
  [../]

  [./boundary_top_grounding]
    type = DirichletBC
    boundary = 'top'
    variable = potential_E_int
    value = 0.0
  [../]

  [./boundary_bottom_grounding]
    type = DirichletBC
    boundary = 'bottom'
    variable = potential_E_int
    value = 0.0
  [../]
[]

[Postprocessors]

  ####################################################
  ##
  ##  Postprocessors (integrations over the
  ##  computational domain) to calculate the total
  ##  energy decomposed into linear combinations of
  ##  the different physics. Note that we only track
  ##  the energy of the ferroelectric layer
  ##  The dielectric layer also can be tracked if you
  ##  want.
  ##
  ####################################################

  [./Fbulk]
    type = BulkEnergyEighth
    execute_on = 'initial timestep_end'
    block = '0'
  [../]
  [./Fwall]
    type = WallEnergy
    execute_on = 'initial timestep_end'
    block = '0'
  [../]
  [./Felec]
    type = ElectrostaticEnergy
    execute_on = 'initial timestep_end'
    block = '0'
  [../]
  [./Ftotal]
    type = LinearCombinationPostprocessor
    pp_names = 'Fbulk Fwall Felec'
    pp_coefs = ' 1 1 1'
    execute_on = 'initial timestep_end'
  [../]
  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = Ftotal
    execute_on = 'initial timestep_end'
  [../]
[]

[UserObjects]

  ####################################################
  ##
  ##  terminator to end energy evolution when the
  ##  energy difference between subsequent time steps
  ##  is lower than 5e-6
  ##
  ##  NOTE: can fail if the time step is small
  ##
  ####################################################

  [./kill]
    type = Terminator
    expression = 'perc_change <= 1.0e-5'
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
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type  -build_twosided'
    petsc_options_value = '    160            1e-10      1e-8      1e-6       bjacobi       allreduce'
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

  ###########################################
  ##
  ##  dtmax is material dependent!
  ##  for PTO is about 0.8 but BTO/HFO more
  ##  like dt = 3-10
  ##
  ###########################################

  dtmax = 8.0

  l_max_its = 200

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 12
    cutback_factor = 0.75
    linear_iteration_ratio = 1000
    dt = 0.1
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

  [./out]
    type = Exodus
    file_base = out_FE_8nmFilm_Die_4nmLayer_U0
    elemental_as_nodal = true
    interval = 1
  [../]
  [./outCSV]
    type = CSV
    file_base = out_FE_8nmFilm_Die_4nmLayer_U0
  [../]
[]
