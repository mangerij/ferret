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
    ##=
    #############################################

    nx = 40
    ny = 40
    nz = 15

    #############################################
    ##
    ##   Actual spatial coordinates of mesh.
    ##   Jmax - Jmin = nJ/2 for J = x, y, z
    ##   Units are in nanometers
    ##
    #############################################

    xmin = -20.0
    xmax = 20.0
    ymin = -20.0
    ymax = 20.0
    zmin = -5.0
    zmax = 10.0

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
    ##=
    ############################################

    type = ExtraNodesetGenerator
    coord = '-20.0 -20.0 -5.0'
    new_boundary = 100
  [../]

  [subdomains]
    type = SubdomainBoundingBoxGenerator
    input = cnode
    bottom_left = '-20.0 -20.0 -5.0'
    block_id = 1
    top_right = '20.0 20.0 0.0'
    location = INSIDE
  []
  [film_interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = subdomains
    primary_block = 0
    paired_block = 1
    new_boundary = 52
  []
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
[]

[Variables]

  #################################
  ##
  ##  Variable definitions
  ##    P, u, phi, e^global_ij
  ##  and their initial conditions
  ##
  #################################

#  [./global_strain]
#    order = SIXTH
#    family = SCALAR
#  [../]
  [./global_strain_film]
    order = SIXTH
    family = SCALAR
    block = '0'
  [../]

  [./global_strain_sub]
    order = SIXTH
    family = SCALAR
    block = '1'
  [../]


  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -1.0e-5
      max = 1.0e-5
      seed = 17
    [../]
    block = '0'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -1.0e-5
      max = 1.0e-5
      seed = 17
    [../]
    block = '0'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -1.0e-5
      max = 1.0e-5
      seed = 17
    [../]
    block = '0'
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

[Materials]

  #################################################
  ##
  ## add comments
  ##
  ## NOTE: there might be some Legendre transforms
  ##        depending on what approach you use
  ##        -i.e. inhomogeneous strain vs
  ##            homogeneous strain [renormalized]
  ##=
  ##################################################

  [./Landau_P_FE]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '-0.1722883 -0.073 0.75 0.26 0.61 -3.67 0.0 0.0 0.0 0.0'
    block = '0'
  [../]

  ############################################
  ##
  ## add comments
  ##
  ############################################

  [./Landau_G_FE]
    type = GenericConstantMaterial
    prop_names = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '1 0.6 0.0 0.3 0.3'  #dx = 1
#    prop_values = '1 0.42 0.0 0.21 0.21' #dx = 1.2
#    prop_values = '1 0.15 0.0 0.075 0.075' #dx = 1.2
    block = '0'
  [../]

  [./mat_C_FE]
    type = GenericConstantMaterial
    prop_names = 'C11 C12 C44'
    prop_values = '175.0 79.4 111.1'
    block = '0 1'
  [../]
#  [./mat_C_sub]
#    type = GenericConstantMaterial
#    prop_names = 'C11 C12 C44'
#    prop_values = '220.0 34.4 161.1'
#    block = '1'
#  [../]

  [./mat_Q]
    type = GenericConstantMaterial
    prop_names = 'Q11 Q12 Q44'
    prop_values = '0.089 -0.026 0.03375'
    block = '0'
  [../]

  [./mat_q]
    type = GenericConstantMaterial
    prop_names = 'q11 q12 q44'
#    prop_values = '11.4 -0.01438 7.5'
    prop_values = '15.523 0.4522 7.5'
    block = '0'
  [../]

  [./eigen_strain]
    type = ComputeEigenstrain
#    eigen_base = '0 0 0 0 0 0'
#    eigen_base = 'x y z yz xz xy'
    eigen_base = '-0.002 -0.002 0.0012 0  0  0'
    eigenstrain_name = epitaxy
    block = '0'
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
    block = '0 1'
  [../]
  [./film_eigenstrain]
    type = CompositeEigenstrain
    tensors = 'ferro  epitaxy'
    weights = 'weight1 weight2'
    args = 'polar_x polar_y polar_z'
    eigenstrain_name = total_eigenstrain
    block = 0
  [../]
  [./weight1]
    type = DerivativeParsedMaterial
    block = 0
#    expression = '-0.1*(sqrt(polar_x^2 + polar_y^2 + polar_z^2))^2'
    expression = '-0.1'
    property_name = weight1
    coupled_variables = 'polar_x polar_y polar_z'
  [../]
  [./weight2]
    type = DerivativeParsedMaterial
    block = 0
    expression = '1'
    property_name = weight2
    coupled_variables = u_x
  [../]


#  [./elasticity_tensor_0]
#    type = ComputeElasticityTensor
#    fill_method = symmetric9

   ###############################################
   ##
   ## symmetric9 fill_method is (default)
   ##     C11 C12 C13 C22 C23 C33 C44 C55 C66
   ##
   ###############################################

#    C_ijkl = '220.0 34.4 34.4 220.0 34.4 220.0 161.1 161.1 161.1'
#    block = '1'
#  [../]

  [./stress_1]
    type = ComputeLinearElasticStress
  [../]

  [./ferroelectric_eigenstrain]
    type = ComputeFerroelectricStrain
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
    eigenstrain_name = 'ferro'
    block = '0'
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
    block = '0 1'
  [../]
[]


[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [./all_FE]
        strain = SMALL
        add_variables = true
        incremental = false
        eigenstrain_names = 'total_eigenstrain'
        global_strain = global_strain
        generate_output = 'strain_xx strain_xy strain_xz strain_yy strain_yz strain_zz vonmises_stress'
	block = '0'
      [../]

      [./all_substrate]
        strain = SMALL
        add_variables = true
#        eigenstrain_names = 'epitaxy'
        global_strain = global_strain
        incremental = false
        generate_output = 'strain_xx strain_xy strain_xz strain_yy strain_yz strain_zz vonmises_stress'
        block = '1'
      [../]
    []

    # GlobalStrain action for generating the objects associated with the global
    # strain calculation and associated displacement visualization


    [./GlobalStrain]
      [./global_strain_film]
        scalar_global_strain = global_strain_film
        displacements = 'u_x u_y u_z'
        auxiliary_displacements = 'disp_x disp_y disp_z'
        global_displacements = 'ug_x ug_y ug_z'
        block = '0'
      [../]
      [./global_strain_sub]
        scalar_global_strain = global_strain_sub
        displacements = 'u_x u_y u_z'
        auxiliary_displacements = 'disp_x disp_y disp_z'
        global_displacements = 'ug_x ug_y ug_z'
        block = '1'
      [../]
    [../]


#    [./GlobalStrain]
#      [./global_strain]
#        scalar_global_strain = global_strain
#        displacements = 'u_x u_y u_z'
#        auxiliary_displacements = 'disp_x disp_y disp_z'
#        global_displacements = 'ug_x ug_y ug_z'
#      [../]
#    [../]
  []
[]

[Kernels]

  ###############################################
  ##
  ## Physical Kernel operators
  ## to enforce TDLGD evolution
  ##
  ###############################################

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
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_x
    component = 0
    block = '0'
    u_x = u_x
    u_y = u_y
    u_z = u_z
  [../]
  [./electrostr_polar_coupled_y]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_y
    component = 1
    block = '0'
    u_x = u_x
    u_y = u_y
    u_z = u_z
  [../]
  [./electrostr_polar_coupled_z]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_z
    component = 2
    block = '0'
    u_x = u_x
    u_y = u_y
    u_z = u_z
  [../]

  [./polar_x_time]
    type = TimeDerivativeScaled
    variable = polar_x
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


[BCs]
  [./Periodic]
    [./xy]
      auto_direction = 'x y'
      variable = 'u_x u_y u_z polar_x polar_y polar_z'
    [../]
  [../]


  # fix center point location
  [./centerfix_x]
    type = DirichletBC
    boundary = 'back'
    variable = u_x
    value = 0
  [../]
  [./centerfix_y]
    type = DirichletBC
    boundary = 'back'
    variable = u_y
    value = 0
  [../]
  [./centerfix_z]
    type = DirichletBC
    boundary = 'back'
    variable = u_z
    value = 0
  [../]
# [./test_z]
#   type = DirichletBC
#   boundary = 'front'
#   variable = u_z
#   value = -0.12
# [../]
#  [./Pressure]
#    [./top]
#      boundary = 'front'
#      factor = -0.1
#      displacements = 'u_x u_y u_z'
#    [../]
#  [../]
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
    u_x = u_x
    u_y = u_y
    u_z = u_z
    block = '0'
  [../]
  [./Ftotal]
    type = LinearCombinationPostprocessor
    pp_names = 'Fbulk Fwall Fcoupled'
    pp_coefs = ' 1 1 1'
    execute_on = 'timestep_end'
  [../]

  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = Ftotal
    execute_on = 'timestep_end'
  [../]

  [./max]
    type = ElementExtremeValue
    variable = 'polar_z'
    execute_on = 'timestep_end'
    block = '0'
  []
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

  # [./kill]
  #   type = Terminator
  #   expression = 'perc_change <= 1.0e-8'
  #  [../]
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
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    160             1e-8        1e-6      1e-5        bjacobi'
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
  dtmax = 2.0
  end_time = 500
  l_max_its = 200
  automatic_scaling = true
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 8
    growth_factor = 1.1
    cutback_factor = 0.8
    linear_iteration_ratio = 1000
    dt = 0.2
  [../]
  verbose = true
  nl_max_its = 16
[]
#
# [Debug]
#   show_var_residual_norms = true
# []

[Outputs]

  ###############################################
  ##===============
  ##  Output options
  ##
  ###############################################

  print_linear_residuals = false
  perf_graph = false

  [./out]
    type = Exodus
    file_base = out_pto_film_LQChen_neg0
    #elemental_as_nodal = true
    time_step_interval = 3
  [../]
[]
