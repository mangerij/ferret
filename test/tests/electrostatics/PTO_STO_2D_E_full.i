TC   = 25.0

Lx     = 6.0
t_pto  = 2.0
t_sto  = 2.0
Ly     = ${fparse 2.0*t_pto + t_sto}
dx     = 0.25
nx     = ${fparse int(Lx/dx + 0.5)}
ny     = ${fparse int(Ly/dx + 0.5)}
y0     = ${t_pto}
y1     = ${fparse t_pto + t_sto}

lam     = 6.0
P_seed  = 0.3
Px_seed = 0.01

eps0    = 0.0088542
eps_r   = 10.0
eps_sto = 300.0

a_pto = 3.957
a_sto = 3.905
a_sub = 3.944
um_pto = ${fparse (a_sub - a_pto)/a_pto}
um_sto = ${fparse (a_sub - a_sto)/a_sto}

l0     = 1.0
alpha1 = ${fparse 3.8e-4*(TC - 479.0)}
G110   = ${fparse l0*l0*abs(alpha1)}
alpha1_sto = ${fparse 1.0/(2.0*(eps_sto - eps_r)*eps0)}

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${nx}
    ny = ${ny}
    xmin = 0.0
    xmax = ${Lx}
    ymin = 0.0
    ymax = ${Ly}
    elem_type = QUAD4
  []
  [sto_layer]
    type = SubdomainBoundingBoxGenerator
    input = gen
    bottom_left = '0.0 ${y0} 0.0'
    top_right = '${Lx} ${y1} 0.0'
    block_id = 1
  []
  [pin]
    type = ExtraNodesetGenerator
    input = sto_layer
    new_boundary = 'pin_node'
    coord = '${fparse 0.5*Lx} 0.0 0.0'
    use_closest_node = true
  []
[]

[GlobalParams]
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_E_int = potential_E_int
[]

[Functions]
  [ic_Py]
    type = ParsedFunction
    expression = '${P_seed}*cos(2*pi*x/${lam})'
  []
  [ic_Px]
    type = ParsedFunction
    expression = '${Px_seed}*sin(2*pi*x/${lam})'
  []
[]

[Variables]
  [polar_x]
    order = FIRST
    family = LAGRANGE
    [InitialCondition]
      type = FunctionIC
      function = ic_Px
    []
  []
  [polar_y]
    order = FIRST
    family = LAGRANGE
    [InitialCondition]
      type = FunctionIC
      function = ic_Py
    []
  []
  [polar_z]
    order = FIRST
    family = LAGRANGE
    [InitialCondition]
      type = ConstantIC
      value = 0.0
    []
  []
  [potential_E_int]
    order = FIRST
    family = LAGRANGE
  []
  [u_x]
    order = FIRST
    family = LAGRANGE
  []
  [u_y]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [Pmag]
    order = FIRST
    family = LAGRANGE
  []
  [wE]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_yy]
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
    execute_on = 'initial timestep_end'
  []
  [wE]
    type = WallEnergyDensity
    variable = wE
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
  [strain_yy]
    type = RankTwoAux
    variable = strain_yy
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 1
    execute_on = 'initial timestep_end'
  []
[]

[Materials]
  [Landau_P_pto]
    type = GenericConstantMaterial
    prop_names  = 'alpha1     alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '${alpha1}  -0.073  0.75    0.26     0.61     -3.7     0.0       0.0       0.0       0.0'
    block = '0'
  []
  [Landau_G_pto]
    type = GenericConstantMaterial
    prop_names  = 'G110     G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '${G110}  0.6      0.0      0.3      0.3'
    block = '0'
  []
  [mat_Q_pto]
    type = GenericConstantMaterial
    prop_names  = 'Q11   Q12    Q44'
    prop_values = '0.089 -0.026 0.03375'
    block = '0'
  []
  [mat_C_pto]
    type = GenericConstantMaterial
    prop_names  = 'C11     C12    C44'
    prop_values = '174.269 79.029 47.62'
    block = '0'
  []
  [elasticity_tensor_pto]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '174.269 79.029 79.029 174.269 79.029 174.269 47.62 47.62 47.62'
    block = '0'
  []
  [misfit_pto]
    type = GenericConstantRankTwoTensor
    tensor_name = global_strain
    tensor_values = '${um_pto} 0 0   0 0 0   0 0 ${um_pto}'
    block = '0'
  []
  [ferro_pto]
    type = ComputeCubicParentElectrostrictiveStrain
    eigenstrain_name = ferro
    block = '0'
  []

  [Landau_P_sto]
    type = GenericConstantMaterial
    prop_names  = 'alpha1          alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '${alpha1_sto}   0.0     0.0     0.0      0.0      0.0      0.0       0.0       0.0       0.0'
    block = '1'
  []
  [Landau_G_sto]
    type = GenericConstantMaterial
    prop_names  = 'G110     G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '${G110}  0.6      0.0      0.3      0.3'
    block = '1'
  []
  [mat_Q_sto]
    type = GenericConstantMaterial
    prop_names  = 'Q11     Q12      Q44'
    prop_values = '0.04575 -0.01350 0.00957'
    block = '1'
  []
  [mat_C_sto]
    type = GenericConstantMaterial
    prop_names  = 'C11   C12   C44'
    prop_values = '318.1 102.5 123.5'
    block = '1'
  []
  [elasticity_tensor_sto]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '318.1 102.5 102.5 318.1 102.5 318.1 123.5 123.5 123.5'
    block = '1'
  []
  [misfit_sto]
    type = GenericConstantRankTwoTensor
    tensor_name = global_strain
    tensor_values = '${um_sto} 0 0   0 0 0   0 0 ${um_sto}'
    block = '1'
  []
  [ferro_sto]
    type = ComputeCubicParentElectrostrictiveStrain
    eigenstrain_name = ferro
    block = '1'
  []

  [permittivity]
    type = GenericConstantMaterial
    prop_names  = 'permittivity'
    prop_values = '${fparse eps_r*eps0}'
  []
  [strain]
    type = ComputeSmallStrain
    displacements = 'u_x u_y'
    eigenstrain_names = 'ferro'
    global_strain = global_strain
  []
  [stress]
    type = ComputeLinearElasticStress
  []
[]

[Kernels]
  [time_x]
    type = TimeDerivative
    variable = polar_x
  []
  [time_y]
    type = TimeDerivative
    variable = polar_y
  []
  [time_z]
    type = TimeDerivative
    variable = polar_z
  []

  [poisson]
    type = Electrostatics
    variable = potential_E_int
  []
  [divP]
    type = PolarElectricEStrong
    variable = potential_E_int
  []
  [Ephi_x]
    type = PolarElectricPStrong
    variable = polar_x
    component = 0
  []
  [Ephi_y]
    type = PolarElectricPStrong
    variable = polar_y
    component = 1
  []
  [Ephi_z]
    type = PolarElectricPStrong
    variable = polar_z
    component = 2
  []

  [div_x]
    type = StressDivergenceTensors
    variable = u_x
    displacements = 'u_x u_y'
    component = 0
  []
  [div_y]
    type = StressDivergenceTensors
    variable = u_y
    displacements = 'u_x u_y'
    component = 1
  []

  [bed_x]
    type = BulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
  []
  [bed_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
  []
  [bed_z]
    type = BulkEnergyDerivativeEighth
    variable = polar_z
    component = 2
  []

  [walled_x]
    type = WallEnergyDerivative
    variable = polar_x
    component = 0
  []
  [walled_y]
    type = WallEnergyDerivative
    variable = polar_y
    component = 1
  []
  [walled_z]
    type = WallEnergyDerivative
    variable = polar_z
    component = 2
  []

  [electrostr_x]
    type = CubicParentElasticPDerivative
    variable = polar_x
    component = 0
    displacements = 'u_x u_y'
  []
  [electrostr_y]
    type = CubicParentElasticPDerivative
    variable = polar_y
    component = 1
    displacements = 'u_x u_y'
  []
  [electrostr_z]
    type = CubicParentElasticPDerivative
    variable = polar_z
    component = 2
    displacements = 'u_x u_y'
  []
[]

[BCs]
  [Periodic]
    [x]
      auto_direction = 'x'
      variable = 'polar_x polar_y polar_z potential_E_int u_x u_y'
    []
  []
  [phi_pin]
    type = DirichletBC
    variable = potential_E_int
    boundary = 'pin_node'
    value = 0
  []
  [clamp_x]
    type = DirichletBC
    variable = u_x
    boundary = 'bottom'
    value = 0
  []
  [clamp_y]
    type = DirichletBC
    variable = u_y
    boundary = 'bottom'
    value = 0
  []
[]

[Postprocessors]
  [Fbulk]
    type = BulkEnergyEighth
    execute_on = 'initial timestep_end'
  []
  [Fwall]
    type = WallEnergy
    execute_on = 'initial timestep_end'
  []
  [Felec]
    type = ElectrostaticEnergy
    execute_on = 'initial timestep_end'
  []
  [Felastic]
    type = CubicParentElasticEnergy
    execute_on = 'initial timestep_end'
  []
  [Ftotal]
    type = LinearCombinationPostprocessor
    pp_names = 'Fbulk Fwall Felastic Felec'
    pp_coefs = '1 1 1 1'
    execute_on = 'timestep_end'
  []
  [Pmag_pto]
    type = ElementAverageValue
    variable = Pmag
    block = '0'
    execute_on = 'initial timestep_end'
  []
  [Pmag_sto]
    type = ElementAverageValue
    variable = Pmag
    block = '1'
    execute_on = 'initial timestep_end'
  []
  [exx_pto]
    type = ElementAverageValue
    variable = strain_xx
    block = '0'
    execute_on = 'initial timestep_end'
  []
  [exx_sto]
    type = ElementAverageValue
    variable = strain_xx
    block = '1'
    execute_on = 'initial timestep_end'
  []
  [phi_max]
    type = NodalExtremeValue
    variable = potential_E_int
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [phi_min]
    type = NodalExtremeValue
    variable = potential_E_int
    value_type = min
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
  end_time = 5.0
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5
    growth_factor = 1.3
    cutback_factor = 0.8
    optimal_iterations = 8
    linear_iteration_ratio = 1000
  []
  dtmax = 5.0

  num_steps = 4
[]

[Outputs]
  print_linear_residuals = false
  [exo]
    type = Exodus
  []
[]
