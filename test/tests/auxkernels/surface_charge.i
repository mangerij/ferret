[Mesh]
  [gen]

    type = GeneratedMeshGenerator
    dim = 3

    nx = 4
    ny = 4
    nz = 8

    xmin = -1.0
    xmax = 1.0
    ymin = -1.0
    ymax = 1.0
    zmin = -2.0
    zmax = 2.0

    elem_type = HEX8
  []
  [./cnode]
    input = gen

    type = ExtraNodesetGenerator
    coord = '-1.0 -1.0 -2.0'
    new_boundary = 100
  [../]

  [subdomains]
    type = SubdomainBoundingBoxGenerator
    input = cnode
    bottom_left = '-1.0 -1.0 -2.0'
    block_id = 1
    top_right = '1.0 1.0 0'
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

[AuxVariables]

  [./disp_x]
    block = '0 1'
  [../]
  [./disp_y]
    block = '0 1'
  [../]
  [./disp_z]
    block = '0 1'
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

  [./divP]
    order = CONSTANT
    family = MONOMIAL
    block = '0'
  [../]

  [./surfP]
    order = CONSTANT
    family = MONOMIAL
  [../]

  ##  P.n on the film/substrate interface (52), the remaining piece of
  ##  the closed surface of block 0.  It needs its OWN variable: a face
  ##  value is written into the element dof, so two boundaries sharing
  ##  an element cannot share a CONSTANT MONOMIAL variable.
  [./surfP_52]
    order = CONSTANT
    family = MONOMIAL
    block = '0'

    ##  This variable exists only to feed the divergence-theorem assertion
    ##  in [UserObjects], so keep it out of the Exodus file.  That makes the
    ##  whole check strictly additive and leaves the gold untouched.
    outputs = none
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

  [./divP]
    type = DivP
    variable = divP
  [../]

  [./surfP]
    type = SurfaceChargeP
    variable = surfP
    boundary = '107'
  [../]
  [./surfP_52]
    type = SurfaceChargeP
    variable = surfP_52
    boundary = '52'
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

  [./mat_Q]
    type = GenericConstantMaterial
    prop_names = 'Q11 Q12 Q44'
    prop_values = '0.089 -0.026 0.03375'
    block = '0 1'
  [../]

  [./ferro]
    type = ComputeCubicParentElectrostrictiveStrain
    eigenstrain_name = ferro
    block = '0'
  [../]

  [./eigen_strain]
    type = ComputeEigenstrain
    eigen_base = '1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0'
    eigenstrain_name = eigenstrain
    prefactor = 0.0
    block = '1'
  [../]

  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9

    C_ijkl = '175.0 79.4 79.4 175.0 79.4 175.0 111.1 111.1 111.1'
  [../]
  [./strain_film]
    type = ComputeSmallStrain
    global_strain = global_strain
    eigenstrain_names = 'ferro'
    block = '0'
  [../]
  [./strain_substrate]
    type = ComputeSmallStrain
    global_strain = global_strain
    eigenstrain_names = 'eigenstrain'
    block = '1'
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

    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '0.08854187'
  [../]
[]

[Kernels]

  [./div_x]
    type = StressDivergenceTensors
    variable = u_x
    component = 0
    use_displaced_mesh = false
  [../]
  [./div_y]
    type = StressDivergenceTensors
    variable = u_y
    component = 1
    use_displaced_mesh = false
  [../]
  [./div_z]
    type = StressDivergenceTensors
    variable = u_z
    component = 2
    use_displaced_mesh = false
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

  [./electrostr_polar_coupled_x]
    type = CubicParentElasticPDerivative
    variable = polar_x
    component = 0
    displacements = 'u_x u_y u_z'
    block = '0'
  [../]
  [./electrostr_polar_coupled_y]
    type = CubicParentElasticPDerivative
    variable = polar_y
    component = 1
    displacements = 'u_x u_y u_z'
    block = '0'
  [../]
  [./electrostr_polar_coupled_z]
    type = CubicParentElasticPDerivative
    variable = polar_z
    component = 2
    displacements = 'u_x u_y u_z'
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
    variable=polar_y
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
    type = CubicParentElasticEnergy
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
    pp_names = 'Fbulk Fwall Felastic Felec'
    pp_coefs = '0.160218 0.160218 0.160218 0.160218'
    execute_on = 'timestep_end'
  [../]


  ###############################################
  ##
  ##  Gold-independent correctness check on the bound-charge
  ##  capability: the divergence theorem
  ##
  ##      int_V div(P) dV  ==  oint_dV P.n dS
  ##
  ##  over the FILM (block 0, z in [0,2]).  The closed boundary of
  ##  block 0 is the free surface 107 (z=+2), the film/substrate
  ##  interface 52 (z=0), and the four lateral faces, which cancel
  ##  pairwise under the x/y periodic BCs.
  ##
  ##  NOTE: boundary 108 is the SUBSTRATE bottom (z=-2) and lies on
  ##        block 1, where P is not defined.  It is NOT part of this
  ##        surface -- adding it is an error, not a refinement.
  ##
  ##  Sign convention: DivP returns div(P) and SurfaceChargeP returns
  ##  P.n, so the identity carries no minus sign.  (The bound charges
  ##  themselves are rho_b = -div(P) and sigma_b = +P.n.)
  ##
  ##  Verified 2026-09-22 on the full closed surface: residual ~1e-18,
  ##  i.e. ~1e-17 relative to the sum of the absolute face contributions.
  ##  The 1e-6 tolerance below is loose because the lateral faces are
  ##  omitted here; they cancel only to solver tolerance (~1e-8), not to
  ##  round-off.  Restricting sidesets to block 0 so they could be
  ##  included would break Periodic/auto_direction, which is not worth it.
  ##
  ##  outputs = none on purpose: this check must not perturb the Exodus
  ##  gold file, so it adds no global variables to the output.
  ##
  ###############################################

  [./volDivP]
    type = ElementIntegralVariablePostprocessor
    variable = divP
    block = '0'
    execute_on = 'timestep_end'
    outputs = none
  [../]
  [./surfInt_107]
    type = SideIntegralVariablePostprocessor
    variable = surfP
    boundary = '107'
    execute_on = 'timestep_end'
    outputs = none
  [../]
  [./surfInt_52]
    type = SideIntegralVariablePostprocessor
    variable = surfP_52
    boundary = '52'
    execute_on = 'timestep_end'
    outputs = none
  [../]
  [./closedSurfP]
    type = LinearCombinationPostprocessor
    pp_names = 'surfInt_107 surfInt_52'
    pp_coefs = '1 1'
    execute_on = 'timestep_end'
    outputs = none
  [../]
  [./divThmResidual]
    type = LinearCombinationPostprocessor
    pp_names = 'closedSurfP volDivP'
    pp_coefs = '1 -1'
    execute_on = 'timestep_end'
    outputs = none
  [../]

[]

[UserObjects]

  [./global_strain_uo]
    type = GlobalATiO3MaterialRVEUserObject
    use_displaced_mesh = false
    execute_on = 'Initial Linear Nonlinear'
    applied_stress_tensor = '0.0 0.0 0.0 0.0 0.0 0.0'
    block = '0'
  [../]


  ###############################################
  ##
  ##  Assert the divergence theorem every step.  Unlike the Exodiff
  ##  this does not care what the solution IS, only that DivP and
  ##  SurfaceChargeP are consistent, so it survives a regold.
  ##  error_level = ERROR is required: fail_mode = HARD on its own
  ##  stops the run but still exits 0, which a test would not catch.
  ##
  ###############################################

  [./divthm_assert]
    type = Terminator
    expression = 'abs(divThmResidual) > 1.0e-6'
    fail_mode = HARD
    error_level = ERROR
    execute_on = 'TIMESTEP_END'
    message = 'Divergence theorem violated: oint P.n dS != int div(P) dV over block 0. DivP/SurfaceChargeP are inconsistent.'
  [../]
[]

[Preconditioning]

  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type  -build_twosided'
    petsc_options_value = '    80             1e-8        1e-5      1e-5       bjacobi      allreduce'
  [../]
[]

[Executioner]

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
  num_steps = 6
[]

[Outputs]

  print_linear_residuals = false
  perf_graph = false

  [./out]
    type = Exodus
    file_base = out_surface_charge
    elemental_as_nodal = true
    time_step_interval = 1
  [../]
[]
