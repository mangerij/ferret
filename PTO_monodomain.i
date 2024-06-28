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
    prop_values = '-0.1722883 -0.073 0.75 0.26 0.61 -3.67 0.0 0.0 0.0 0.0'
#    prop_values = '-0.1722883 0.42 0.735 0.26 0.61 -3.67 0.0 0.0 0.0 0.0'

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
    prop_values = '0.173 0.6 0.0 0.3 0.3'
  [../]

  [./mat_C]
    type = GenericConstantMaterial
    prop_names = 'C11 C12 C44'
    prop_values = '175.0 79.4 111.1'
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
    prop_values = '-0.089 0.026 -0.03375'
  [../]
  [./mat_q]
    type = GenericConstantMaterial
    prop_names = 'q11 q12 q44'
    prop_values = '-11.4 -0.01438 -7.5'
  [../]

    [./eigen_strain]
      type = ComputeEigenstrain
      eigen_base = '0 0 0 0 0 0'
  #    eigen_base = 'x y z yz xz xy'
  #    eigen_base = '0.004 0.004 -0.0024 0  0  0'
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
      block = '0'
    [../]
    [./film_eigenstrain]
      type = CompositeEigenstrain
      tensors = 'ferro  epitaxy'
      weights = 'weight1 weight2'
      eigenstrain_name = total_eigenstrain
      block = 0
      coupled_variables = u_x
    [../]
    [./weight1]
      type = DerivativeParsedMaterial
      block = 0
      expression = '-1'
      property_name = weight1
      coupled_variables = u_x
    [../]
    [./weight2]
      type = DerivativeParsedMaterial
      block = 0
      expression = '1'
      property_name = weight2
      coupled_variables = u_x
    [../]

      [./ferroelectric_eigenstrain]
        type = ComputeFerroelectricStrain
        polar_x = polar_x
        polar_y = polar_y
        polar_z = polar_z
        eigenstrain_name = 'ferro'
        block = '0'
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
    []

    # GlobalStrain action for generating the objects associated with the global
    # strain calculation and associated displacement visualization

    [./GlobalStrain]
      [./global_strain_film]
        scalar_global_strain = global_strain
        displacements = 'u_x u_y u_z'
        auxiliary_displacements = 'disp_x disp_y disp_z'
        global_displacements = 'ug_x ug_y ug_z'
        block = '0'
      [../]
    [../]
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
    [./E2_x]
      type = CorrectedElectrostrictiveCouplingPolarDerivative
      variable = polar_x
      component = 0
      block = '0'
    [../]
    [./E2_y]
      type = CorrectedElectrostrictiveCouplingPolarDerivative
      variable = polar_y
      component = 1
      block = '0'
    [../]
    [./E2_z]
      type = CorrectedElectrostrictiveCouplingPolarDerivative
      variable = polar_z
      component = 2
      block = '0'
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
    file_base = out_pto_monodomain
  [../]
[]
