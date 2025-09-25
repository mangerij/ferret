
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 4
    ny = 4
    nz = 4
    xmin = -6.0
    xmax = 6.0
    ymin = -6.0
    ymax = 6.0
    zmin = -2.0
    zmax = 2.0
    elem_type = HEX8
  []
  [./cnode]
    input = gen
    type = ExtraNodesetGenerator
    coord = '-6.0 -6.0 -2.0'
    new_boundary = 100
  [../]

  # additional boundary sideset (one node) to zero one of the elastic displacement vectors - eliminates rigid body translations from the degrees of freedom
[]

[GlobalParams]
  len_scale = 1.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_E_int = potential_E_int

  displacements = 'u_x u_y u_z'

  u_x = u_x
  u_y = u_y
  u_z = u_z

[]



[Variables]
  [./global_strain]
    order = SIXTH
    family = SCALAR
  [../]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-4
      max = 0.01e-4
      seed = 6
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-4
      max = 0.01e-4
      seed = 6
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-4
      max = 0.01e-4
      seed = 6
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
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
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

  [./Fb_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Fb_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Fb_z]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./Felstr_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Felstr_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Felstr_z]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./Felec_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Felec_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Felec_z]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./Fw_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Fw_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Fw_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
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

  [./fbx]
    type = MicroforceBulkEnergy
    variable = Fb_x
    component = 0
  [../]
  [./fby]
    type = MicroforceBulkEnergy
    variable = Fb_y
    component = 1
  [../]
  [./fbz]
    type = MicroforceBulkEnergy
    variable = Fb_z
    component = 2
  [../]

  [./felstrx]
    type = MicroforceElectrostrictiveCouplingEnergy
    variable = Felstr_x
    component = 0
  [../]
  [./felstry]
    type = MicroforceElectrostrictiveCouplingEnergy
    variable = Felstr_y
    component = 1
  [../]
  [./felstrz]
    type = MicroforceElectrostrictiveCouplingEnergy
    variable = Felstr_z
    component = 2
  [../]

  [./felecx]
    type = MicroforceElectrostaticEnergy
    variable = Felec_x
    component = 0
  [../]
  [./felecy]
    type = MicroforceElectrostaticEnergy
    variable = Felec_y
    component = 1
  [../]
  [./felecz]
    type = MicroforceElectrostaticEnergy
    variable = Felec_z
    component = 2
  [../]

  [./fwx]
    type = MicroforceWallEnergy
    variable = Fw_x
    component = 0
  [../]
  [./fwy]
    type = MicroforceWallEnergy
    variable = Fw_y
    component = 1
  [../]
  [./fwz]
    type = MicroforceWallEnergy
    variable = Fw_z
    component = 2
  [../]
[]

[ScalarKernels]
  [./global_strain]
    type = GlobalStrain
    variable = global_strain
    global_strain_uo = global_strain_uo
    use_displaced_mesh = false
  [../]
[]

[Materials]

  [./Landau_P]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '-0.1722883 0.42 0.735 0.26 0.61 -3.67 0 0 0 0'
  [../]

  [./Landau_G]
    type = GenericConstantMaterial
    prop_names = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '0.173 0.6 0.0 0.3 0.3'
  [../]

  [./mat_C]
    type = GenericConstantMaterial
    prop_names = 'C11 C12 C44'
    prop_values = '175 79.4 111.1'
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
  [../]
  [./mat_q]
    type = GenericConstantMaterial
    prop_names = 'q11 q12 q44'
    prop_values = '-11.4 -0.01438 -7.5'
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

  [./slab_ferroelectric]
    type = ComputeElectrostrictiveTensor
    Q_mnkl = '-0.089 0.026 0.026 -0.089 0.026 -0.089 -0.03375 -0.03375 -0.03375'
    C_ijkl = '175.0 79.4 79.4 175.0 79.4 175.0 111.1 111.1 111.1'
  [../]

  [./permitivitty_1]
    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '0.08854187'
  [../]
[]


[Kernels]

  #note below we use a strain-renormalized functional for lead titanate (this is different than the stress-free functionals typically used)
  # they are related by a Legendre transformation

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
      variable = 'u_x u_y u_z polar_x polar_y polar_z potential_E_int'
    [../]
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
  #[./elapsed]
  #  type = PerfGraphData
  #  section_name = "Root"
  #  data_type = total
  #[../]
  #[./nodes]
  #  type = NumNodes
  #[../]
[]

[UserObjects]
  [./global_strain_uo]
    type = GlobalATiO3MaterialRVEUserObject
    use_displaced_mesh = false
    execute_on = 'Initial Linear Nonlinear'
  [../]
  [./kill]
    type = Terminator
    expression = 'perc_change <= 1.0e-6'
   [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type  -build_twosided'
    petsc_options_value = '    160               1e-8      1e-8      1e-8          bjacobi       allreduce'
  [../]
[]

#[Debug]
#  show_var_residual_norms = true
#[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  scheme = 'implicit-euler'
  dtmin = 1e-13
  dtmax = 0.8
  l_max_its = 200
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 8
    cutback_factor = 0.75
    linear_iteration_ratio = 1000
    dt = 0.3
  [../]
  num_steps = 10
[]

[Outputs]
  print_linear_residuals = false
  perf_graph = true
  [./out]
    type = Exodus
    file_base = out_microforce_test
    elemental_as_nodal = true
    time_step_interval = 1
  [../]
  [./outCSV]
    type = CSV
    file_base = out_microforce_test
  [../]
[]
