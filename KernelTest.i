[Mesh]
  [gen]
    ############################################
    ##
    ##  Type and dimension of the mesh
    ##
    ############################################

    type = GeneratedMeshGenerator
    dim = 3

    nx = 6
    ny = 6
    nz = 6

    xmin = -1
    xmax = 1
    ymin = -1
    ymax = 1
    zmin = -1
    zmax = 1

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
    coord = '-1 -1 -1'
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
        eigenstrain_names = 'spont_polar'
        global_strain = global_strain
        generate_output = ' strain_xx  strain_yy strain_zz strain_xy strain_yz'
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
      min = -1e-8
      max = 1e-8
    [../]
    block = '0'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -1e-7
      max = 1e-7
    [../]
    block = '0'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 1e-6
      max = 1e-5
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

  [./electrostr_polar_coupled_x]
    type = ElectrostrictiveCouplingPolarDerivativeRefactor
    variable = polar_x
    component = 0
    block = '0'
  [../]
  [./electrostr_polar_coupled_y]
    type = ElectrostrictiveCouplingPolarDerivativeRefactor
    variable = polar_y
    component = 1
    block = '0'
  [../]
  [./electrostr_polar_coupled_z]
    type = ElectrostrictiveCouplingPolarDerivativeRefactor
    variable = polar_z
    component = 2
    block = '0'
  [../]
  # [./electrostr_polar_coupled_x]
  #   type = ElectrostrictiveCouplingPolarDerivative
  #   variable = polar_x
  #   component = 0
  #   block = '0'
  # [../]
  # [./electrostr_polar_coupled_y]
  #   type = ElectrostrictiveCouplingPolarDerivative
  #   variable = polar_y
  #   component = 1
  #   block = '0'
  # [../]
  # [./electrostr_polar_coupled_z]
  #   type = ElectrostrictiveCouplingPolarDerivative
  #   variable = polar_z
  #   component = 2
  #   block = '0'
  # [../]


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
  ## Electrostrictive coupling
  ##
  ##################################################

  [./mat_Q]
    type = GenericConstantMaterial
    prop_names = 'Q11 Q12 Q44'
    prop_values = '0.11 -0.045 0.029'
    block = '0'
  [../]

  [./mat_q]
    type = GenericConstantMaterial
    prop_names = 'q11 q12 q44'
    prop_values = '14.2 -0.74 1.57'
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
    eigenstrain_name = 'spont_polar'
    block = '0'
[../]
[]


[BCs]
  [./Periodic]
    [./xyz]
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
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type -pc_factor_mat_solver_type'
    petsc_options_value = '    120             1e-10        1e-6      1e-6       lu   superlu_dist'
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
  dtmax = 3.6

  l_max_its = 200
  #automatic_scaling = true
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

  num_steps = 100
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
    file_base = out_monoBTO_eig
    elemental_as_nodal = true
    interval = 5
  [../]
[]
