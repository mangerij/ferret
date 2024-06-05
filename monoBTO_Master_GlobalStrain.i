[Mesh]
  [gen]
    ############################################
    ##
    ##  Type and dimension of the mesh
    ##
    ############################################

    type = GeneratedMeshGenerator
    dim = 3

    nx = 3
    ny = 3
    nz = 3

    xmin = -0.5
    xmax = 0.5
    ymin = -0.5
    ymax = 0.5
    zmin = -0.5
    zmax = 0.5

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
    ##   NOTE: This must conform with the above
    ##         [Mesh] block settings
    ##
    ############################################

    type = ExtraNodesetGenerator
    coord = '-0.5 -0.5 -0.5'
    new_boundary = 100
  [../]
[]

[GlobalParams]
  len_scale = 1.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

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

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [./all]
        strain = SMALL
        #incremental = true
        add_variables = true
        global_strain = global_strain
        eigenstrain_names = 'ferro'
        generate_output = ' strain_xx elastic_strain_xx elastic_strain_yy elastic_strain_zz elastic_strain_xy elastic_strain_yz'
      [../]
    []
    # GlobalStrain action for generating the objects associated with the global
    # strain calculation and associated displacement visualization
    [./GlobalStrain]
      [./global_strain]
        scalar_global_strain = global_strain
        displacements = 'u_x u_y u_z'
        auxiliary_displacements = 'disp_x disp_y disp_z'
        global_displacements = 'ug_x ug_y ug_z'
      [../]
    [../]
  []
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
      min = -1e-4
      max = 1e-4
    [../]
    block = '0'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -1e-4
      max = 1e-4
    [../]
    block = '0'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 1e-3
      max = 1e-2
    [../]
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

[Kernels]

  ###############################################
  ##
  ## Physical Kernel operators
  ## to enforce TDLGD evolution
  ##
  ###############################################


  #Elastic Problem

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
  #
  # [./electrostr_ux]
  #   type = ElectrostrictiveCouplingDispDerivative
  #   variable = u_x
  #   component = 0
  #   block = '0'
  # [../]
  # [./electrostr_uy]
  #   type = ElectrostrictiveCouplingDispDerivative
  #   variable = u_y
  #   component = 1
  #   block = '0'
  # [../]
  # [./electrostr_uz]
  #   type = ElectrostrictiveCouplingDispDerivative
  #   variable = u_z
  #   component = 2
  #   block = '0'
  # [../]

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

[Materials]

  #################################################
  ##
  ## Landau coefficients from Hlinka and Marton
  ##     Phys. Rev. B. 74, 104014, (2006)
  ##################################################

  [./Landau_P_FE]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '-0.027721 -0.64755 0.323 8.004 4.47 4.91 0.0 0.0 0.0 0.0'
    block = '0'
  [../]

  [./Landau_G_FE]
    type = GenericConstantMaterial
    prop_names = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '0.5 0.51 -0.02 0.02 0.0'
    block = '0'
  [../]

  [./mat_C_FE]
    type = GenericConstantMaterial
    prop_names = 'C11 C12 C44'
    prop_values = '275.0 179.0 54.3'
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
    prop_values = '-0.11 0.045 -0.029'
    block = '0'
  [../]

  [./mat_q]
    type = GenericConstantMaterial
    prop_names = 'q11 q12 q44'
    prop_values = '-14.2 0.74 -1.57'
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

  [./stress_1]
    type = ComputeLinearElasticStress
  [../]

  [./ferroelectric_eigenstrain]
  type = ComputeFerroelectricStrain
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  eigenstrain_name = 'ferro'
[../]
[]


[BCs]
  [./Periodic]
    [./xy]
      auto_direction = 'x y z'
      variable = 'polar_x polar_y polar_z u_x u_y u_z'
    [../]
  [../]

  [./cnode_ux]
    type = DirichletBC
    boundary = '100'
    variable = u_x
    value = 0.0
  [../]
  [./cnode_uy]
    type = DirichletBC
    boundary = '100'
    variable = u_y
    value = 0.0
  [../]
  [./cnode_uz]
    type = DirichletBC
    boundary = '100'
    variable = u_z
    value = 0.0
  [../]
[]

[Postprocessors]

  [./avePz]
    type = ElementAverageValue
    variable = polar_z
    execute_on = 'initial timestep_end'
  [../]


  ###############################################
  ##
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
  [./Ftotal]
    type = LinearCombinationPostprocessor
    pp_names = 'Fbulk Fwall Fcoupled'
    pp_coefs = '0.160218 0.160218 0.160218' #converted to eV
    execute_on = 'timestep_end'
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
  ##  terminator to end energy evolution when the energy difference
  ##  between subsequent time steps is lower than 5e-6
  ##
  ##  NOTE: can fail if the time step is small
  ##
  ###############################################

  [./kill]
    type = Terminator
    expression = 'perc_change <= 1.0e-8'
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
  dtmax = 3.6

  l_max_its = 200
  automatic_scaling = true
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    growth_factor = 1.2
    cutback_factor = 0.75
    linear_iteration_ratio = 1000
    dt = 0.5
  [../]
  verbose = true
  nl_max_its = 20
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
    file_base = out_monoBTO
    elemental_as_nodal = true
    interval = 1
  [../]
  [./outCSV]
    type = CSV
    file_base = out_monoBTO
    execute_on = 'timestep_end'
  [../]
[]
