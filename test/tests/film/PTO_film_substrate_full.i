TC   = 25.0
l0   = 1.0
L    = 4.0
dx   = 1.0
nx   = ${fparse int(L/dx + 0.5)}

h_film = 4.0
h_sub  = 2.0
nz     = ${fparse int((h_film + h_sub)/dx + 0.5)}

alpha1 = ${fparse 3.8e-4*(TC - 479.0)}
G110   = ${fparse l0*l0*abs(alpha1)}

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${nx}
    nz = ${nz}
    xmin = 0.0
    xmax = ${L}
    ymin = 0.0
    ymax = ${L}
    zmin = ${fparse -h_sub}
    zmax = ${h_film}
    elem_type = HEX8
  []

  [subdomains]
    type = SubdomainBoundingBoxGenerator
    input = gen
    block_id = 1
    bottom_left = '0.0 0.0 ${fparse -h_sub}'
    top_right   = '${L} ${L} 0.0'
    location = INSIDE
  []

  [film_interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = subdomains
    primary_block = 0
    paired_block = 1
    new_boundary = 'film_interface'
  []
[]

[GlobalParams]
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  displacements = 'u_x u_y u_z'
[]

[Variables]

  [polar_x]
    order = FIRST
    family = LAGRANGE
    block = 0
    [InitialCondition]
      type = RandomIC
      min = -0.5
      max = 0.5
      seed = 1
    []
  []
  [polar_y]
    order = FIRST
    family = LAGRANGE
    block = 0
    [InitialCondition]
      type = RandomIC
      min = -0.5
      max = 0.5
      seed = 2
    []
  []
  [polar_z]
    order = FIRST
    family = LAGRANGE
    block = 0
    [InitialCondition]
      type = RandomIC
      min = -0.5
      max = 0.5
      seed = 3
    []
  []

  [u_x]
    order = FIRST
    family = LAGRANGE
  []
  [u_y]
    order = FIRST
    family = LAGRANGE
  []
  [u_z]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [Pmag]
    order = FIRST
    family = LAGRANGE
    block = 0
  []
  [strain_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_zz]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [Pmag]
    type = ParsedAux
    variable = Pmag
    coupled_variables = 'polar_x polar_y polar_z'
    expression = 'sqrt(polar_x^2 + polar_y^2 + polar_z^2)'
    block = 0
    execute_on = 'initial timestep_end'
  []

  [strain_xx]
    type = RankTwoAux
    variable = strain_xx
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 0
    execute_on = 'initial timestep_end'
  []
  [strain_zz]
    type = RankTwoAux
    variable = strain_zz
    rank_two_tensor = total_strain
    index_i = 2
    index_j = 2
    execute_on = 'initial timestep_end'
  []
[]

[Materials]

  [Landau_P]
    type = GenericConstantMaterial
    prop_names  = 'alpha1     alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '${alpha1}  -0.073  0.75    0.26     0.61     -3.7     0.0       0.0       0.0       0.0'
    block = 0
  []
  [Landau_G]
    type = GenericConstantMaterial
    prop_names  = 'G110     G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '${G110}  0.6      0.0      0.3      0.3'
    block = 0
  []
  [mat_C]
    type = GenericConstantMaterial
    prop_names  = 'C11     C12    C44'
    prop_values = '174.269 79.029 47.62'
    block = 0
  []
  [mat_Q]
    type = GenericConstantMaterial
    prop_names  = 'Q11   Q12    Q44'
    prop_values = '0.089 -0.026 0.03375'
    block = 0
  []
  [elasticity_tensor_film]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '174.269 79.029 79.029 174.269 79.029 174.269 47.62 47.62 47.62'
    block = 0
  []
  [ferro]
    type = ComputeCubicParentElectrostrictiveStrain
    eigenstrain_name = ferro
    block = 0
  []

  [strain_film]
    type = ComputeSmallStrain
    eigenstrain_names = 'ferro'
    block = 0
  []
  [stress_film]
    type = ComputeLinearElasticStress
    block = 0
  []

  [elasticity_tensor_sub]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '318.0 103.0 103.0 318.0 103.0 318.0 124.0 124.0 124.0'
    block = 1
  []
  [strain_sub]
    type = ComputeSmallStrain
    block = 1
  []
  [stress_sub]
    type = ComputeLinearElasticStress
    block = 1
  []
[]

[Kernels]

  [div_x]
    type = StressDivergenceTensors
    variable = u_x
    component = 0
  []
  [div_y]
    type = StressDivergenceTensors
    variable = u_y
    component = 1
  []
  [div_z]
    type = StressDivergenceTensors
    variable = u_z
    component = 2
  []

  [time_x]
    type = TimeDerivative
    variable = polar_x
    block = 0
  []
  [time_y]
    type = TimeDerivative
    variable = polar_y
    block = 0
  []
  [time_z]
    type = TimeDerivative
    variable = polar_z
    block = 0
  []

  [bed_x]
    type = BulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
    block = 0
  []
  [bed_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
    block = 0
  []
  [bed_z]
    type = BulkEnergyDerivativeEighth
    variable = polar_z
    component = 2
    block = 0
  []

  [walled_x]
    type = WallEnergyDerivative
    variable = polar_x
    component = 0
    block = 0
  []
  [walled_y]
    type = WallEnergyDerivative
    variable = polar_y
    component = 1
    block = 0
  []
  [walled_z]
    type = WallEnergyDerivative
    variable = polar_z
    component = 2
    block = 0
  []

  [electrostr_x]
    type = CubicParentElasticPDerivative
    variable = polar_x
    component = 0
    block = 0
  []
  [electrostr_y]
    type = CubicParentElasticPDerivative
    variable = polar_y
    component = 1
    block = 0
  []
  [electrostr_z]
    type = CubicParentElasticPDerivative
    variable = polar_z
    component = 2
    block = 0
  []
[]

[BCs]

  [Periodic]
    [xy]
      auto_direction = 'x y'
      variable = 'u_x u_y u_z polar_x polar_y polar_z'
    []
  []

  [sub_fix_x]
    type = DirichletBC
    variable = u_x
    boundary = 'back'
    value = 0
  []
  [sub_fix_y]
    type = DirichletBC
    variable = u_y
    boundary = 'back'
    value = 0
  []
  [sub_fix_z]
    type = DirichletBC
    variable = u_z
    boundary = 'back'
    value = 0
  []
[]

[Postprocessors]

  [Fbulk]
    type = BulkEnergyEighth
    block = 0
    execute_on = 'initial timestep_end'
  []
  [Fwall]
    type = WallEnergy
    block = 0
    execute_on = 'initial timestep_end'
  []
  [Felastic]
    type = CubicParentElasticEnergy
    block = 0
    execute_on = 'initial timestep_end'
  []
  [Ftotal]
    type = LinearCombinationPostprocessor
    pp_names = 'Fbulk Fwall Felastic'
    pp_coefs = '1 1 1'
    execute_on = 'timestep_end'
  []

  [exx_film]
    type = ElementAverageValue
    variable = strain_xx
    block = 0
    execute_on = 'initial timestep_end'
  []
  [ezz_film]
    type = ElementAverageValue
    variable = strain_zz
    block = 0
    execute_on = 'initial timestep_end'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -ksp_type -ksp_rtol -ksp_gmres_restart'
    petsc_options_value = 'bjacobi  ilu          gmres     1e-6      100'
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  l_max_its = 200
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5
    growth_factor = 1.3
    cutback_factor = 0.8
    optimal_iterations = 8
    linear_iteration_ratio = 1000
  []
  dtmax = 5.0

  num_steps = 3
[]

[Outputs]
  print_linear_residuals = false
  [exo]
    type = Exodus
    file_base = PTO_film_substrate_full_exo
  []
[]
