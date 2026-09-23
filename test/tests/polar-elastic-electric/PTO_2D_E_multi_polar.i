TC   = 25.0
um   = -0.01
l0   = 1.0
Lx   = 8.0
tf   = 2.0
nx   = 40
ny   = 10
lam  = 8.0
P_seed = 0.3
Px_seed = 0.01

eps_r = 10.0
eps0  = 0.0088542

dt_polar = 0.5

alpha1 = ${fparse 3.8e-4*(TC - 479.0)}
G110   = ${fparse l0*l0*abs(alpha1)}

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${nx}
    ny = ${ny}
    xmin = 0.0
    xmax = ${Lx}
    ymin = 0.0
    ymax = ${tf}
    elem_type = QUAD4
  []
  [pin]
    type = ExtraNodesetGenerator
    input = gen
    new_boundary = 'pin_node'
    coord = '8.0 0.0 0.0'
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
[]

[AuxVariables]
  [u_x]
    order = FIRST
    family = LAGRANGE
  []
  [u_y]
    order = FIRST
    family = LAGRANGE
  []
  [Px_old]
    order = FIRST
    family = LAGRANGE
  []
  [Py_old]
    order = FIRST
    family = LAGRANGE
  []
  [Pz_old]
    order = FIRST
    family = LAGRANGE
  []
  [Pmag]
    order = FIRST
    family = LAGRANGE
  []
  [wE]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [Px_old]
    type = ParsedAux
    variable = Px_old
    coupled_variables = 'polar_x'
    expression = 'polar_x'
    execute_on = 'initial timestep_begin'
  []
  [Py_old]
    type = ParsedAux
    variable = Py_old
    coupled_variables = 'polar_y'
    expression = 'polar_y'
    execute_on = 'initial timestep_begin'
  []
  [Pz_old]
    type = ParsedAux
    variable = Pz_old
    coupled_variables = 'polar_z'
    expression = 'polar_z'
    execute_on = 'initial timestep_begin'
  []
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
[]

[Materials]
  [Landau_P]
    type = GenericConstantMaterial
    prop_names  = 'alpha1     alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '${alpha1}  -0.073  0.75    0.26     0.61     -3.7     0.0       0.0       0.0       0.0'
  []
  [Landau_G]
    type = GenericConstantMaterial
    prop_names  = 'G110     G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '${G110}  0.6      0.0      0.3      0.3'
  []
  [mat_C]
    type = GenericConstantMaterial
    prop_names  = 'C11     C12    C44'
    prop_values = '174.269 79.029 47.62'
  []
  [mat_Q]
    type = GenericConstantMaterial
    prop_names  = 'Q11   Q12    Q44'
    prop_values = '0.089 -0.026 0.03375'
  []
  [permittivity]
    type = GenericConstantMaterial
    prop_names  = 'permittivity'
    prop_values = '${fparse eps_r*eps0}'
  []
  [misfit]
    type = GenericConstantRankTwoTensor
    tensor_name = global_strain
    tensor_values = '${um} 0 0   0 0 0   0 0 ${um}'
  []
  [strain]
    type = ComputeSmallStrain
    displacements = 'u_x u_y'
    global_strain = global_strain
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

  [bed_x]
    type = BulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
    polar_x = Px_old
    polar_y = Py_old
    polar_z = Pz_old
  []
  [bed_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
    polar_x = Px_old
    polar_y = Py_old
    polar_z = Pz_old
  []
  [bed_z]
    type = BulkEnergyDerivativeEighth
    variable = polar_z
    component = 2
    polar_x = Px_old
    polar_y = Py_old
    polar_z = Pz_old
  []

  [electrostr_x]
    type = CubicParentElasticPDerivative
    variable = polar_x
    component = 0
    polar_x = Px_old
    polar_y = Py_old
    polar_z = Pz_old
  []
  [electrostr_y]
    type = CubicParentElasticPDerivative
    variable = polar_y
    component = 1
    polar_x = Px_old
    polar_y = Py_old
    polar_z = Pz_old
  []
  [electrostr_z]
    type = CubicParentElasticPDerivative
    variable = polar_z
    component = 2
    polar_x = Px_old
    polar_y = Py_old
    polar_z = Pz_old
  []
[]

[BCs]
  [Periodic]
    [x]
      auto_direction = 'x'
      variable = 'polar_x polar_y polar_z potential_E_int'
    []
  []
  [phi_pin]
    type = DirichletBC
    variable = potential_E_int
    boundary = 'pin_node'
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
  scheme = implicit-euler
  solve_type = LINEAR
  l_max_its = 300
  dt = ${dt_polar}
  end_time = 1e9
[]

[Outputs]
  print_linear_residuals = false
  [exo]
    type = Exodus
    file_base = PTO_2D_E_multi_exo
  []
[]
