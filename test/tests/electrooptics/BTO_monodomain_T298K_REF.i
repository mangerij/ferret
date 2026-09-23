[Mesh]
  [gen]

    type = GeneratedMeshGenerator
    dim = 3

    nx = 2
    ny = 2
    nz = 2

    xmin = -0.5
    xmax = 0.5
    ymin = -0.5
    ymax = 0.5
    zmin = -0.5
    zmax = 0.5

    elem_type = HEX8
  []
  [./cnode]
    input = gen

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
  potential_E_int = potential_E_int

  displacements = 'u_x u_y u_z'

  u_x = u_x
  u_y = u_y
  u_z = u_z

[]

[Functions]
  [./stripe1]
    type = ParsedFunction
    expression = 0.01*cos(0.08975979010256552*(x))
  [../]
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
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-4
      max = 0.01e-4
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.01
      max = 0.02
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

  [./dpn_1]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./dpn_2]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./dpn_3]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./den_1]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./den_2]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./den_3]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./n_1]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./n_2]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./n_3]
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

 [./dn_p1]
    type = PolarOpticChangeInRefractiveIndex
    component = 0
    variable = dpn_1
    execute_on = 'timestep_end'
  [../]

  [./dn_p2]
    type = PolarOpticChangeInRefractiveIndex
    component = 1
    variable = dpn_2
    execute_on = 'timestep_end'
  [../]

  [./dn_p3]
    type = PolarOpticChangeInRefractiveIndex
    component = 2
    variable = dpn_3
    execute_on = 'timestep_end'
  [../]

 [./dn_e1]
    type = ElastoChangeInRefractiveIndex
    component = 0
    variable = den_1
    u_x = u_x
    u_y = u_y
    u_z = u_z
    execute_on = 'timestep_end'
  [../]

  [./dn_e2]
    type = ElastoChangeInRefractiveIndex
    component = 1
    variable = den_2
    u_x = u_x
    u_y = u_y
    u_z = u_z
    execute_on = 'timestep_end'
  [../]

  [./dn_e3]
    type = ElastoChangeInRefractiveIndex
    component = 2
    variable = den_3
    u_x = u_x
    u_y = u_y
    u_z = u_z
    execute_on = 'timestep_end'
  [../]

  [./n_1_c]
    type = ReworkedRefractiveIndex
    variable = n_1
    elasto = true
    electro = false
    polar = true
    component = 0
    var1 = dpn_1
    var2 = den_1
    execute_on = 'timestep_end'
  [../]

  [./n_2_c]
    type = ReworkedRefractiveIndex
    variable = n_2
    elasto = true
    electro = false
    polar = true
    component = 1
    var1 = dpn_2
    var2 = den_2
    execute_on = 'timestep_end'
  [../]

  [./n_3_c]
    type = ReworkedRefractiveIndex
    variable = n_3
    elasto = true
    electro = false
    polar = true
    component = 2
    var1 = dpn_3
    var2 = den_3
    execute_on = 'timestep_end'
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
    prop_values = '-0.027721 -0.64755 0.323 8.004 4.47 4.91 0.0 0.0 0.0 0.0'
  [../]

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
    eigen_base = '0. 0 0 0 0 0 0 0 0'
    eigenstrain_name = eigenstrain
    prefactor = 0.0
  [../]

  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9

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

  [./slab_ferroelectric]
    type = ComputeElectrostrictiveTensor
    Q_mnkl = '-0.11 0.045 0.045 -0.11 0.045 -0.11 -0.029 -0.029 -0.029'
    C_ijkl = '275.0 179.0 179.0 275.0 179.0 275.0 54.3 54.3 54.3'
  [../]

  [./ref_index]

    type = GenericConstantMaterial
    prop_names = 'n1 n2 n3 n4 n5 n6'
    prop_values = '2.4272 2.4272 2.4272 0.0 0.0 0.0'
  [../]

  [./po_tensor]

    type = GenericConstantMaterial
    prop_names = 'g1111 g1122 g1212'
    prop_values = '0.15 0.038 0.07'
  [../]

  [./eo_tensor]

    type = GenericConstantMaterial
    prop_names = 'p1111 p1122 p1212'
    prop_values = '0.37 0.11 -0.30'
  [../]
[]

[Kernels]

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

     permittivity = 0.08854187
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
[]

[UserObjects]

  [./global_strain_uo]
    type = GlobalATiO3MaterialRVEUserObject
    use_displaced_mesh = false
    execute_on = 'Initial Linear Nonlinear'
  [../]

[]

[Preconditioning]

  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type  -build_twosided'
    petsc_options_value = '    160               1e-10      1e-8      1e-6          bjacobi       allreduce'
  [../]
[]

[Executioner]

  type = Transient
  solve_type = 'PJFNK'
  scheme = 'implicit-euler'
  dtmin = 1e-13

  dtmax = 3.0

  l_max_its = 200

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 8
    cutback_factor = 0.75
    linear_iteration_ratio = 1000
    dt = 0.3
  [../]
  verbose = true
  num_steps = 5
[]

[Outputs]

  print_linear_residuals = false
  perf_graph = false

  [./out]
    type = Exodus
    file_base = test_bire
    elemental_as_nodal = true
  [../]
[]
