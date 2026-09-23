TC   = 25.0
um   = -0.002
C11 = 175.549
C12 = 84.639
C44 = 108.225
umzz = ${fparse -2.0*C12/C11*um}
alpha1   = ${fparse 3.3e-4*(TC - 110.0)}
alpha11  = ${fparse 3.6e-3*(TC - 175.0)}
alpha123 = ${fparse 7.6e-2*(TC - 120.0) + 44.0}
G110 = 1.0
Lx = 2.0
Ly = 2.0
tf = 2.0
nx = 2
ny = 2
nz = 8
dt_polar = 0.05
p0x = 0.0
p0y = 0.0
p0z = 0.0
noise = 0.01

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${ny}
    nz = ${nz}
    xmin = 0.0
    xmax = ${Lx}
    ymin = 0.0
    ymax = ${Ly}
    zmin = 0.0
    zmax = ${tf}
    elem_type = HEX8
  []
[]

[GlobalParams]
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
[]

[Variables]
  [polar_x]
    order = FIRST
    family = LAGRANGE
    [InitialCondition]
      type = RandomIC
      min = ${fparse p0x - noise}
      max = ${fparse p0x + noise}
      seed = 51
    []
  []
  [polar_y]
    order = FIRST
    family = LAGRANGE
    [InitialCondition]
      type = RandomIC
      min = ${fparse p0y - noise}
      max = ${fparse p0y + noise}
      seed = 46
    []
  []
  [polar_z]
    order = FIRST
    family = LAGRANGE
    [InitialCondition]
      type = RandomIC
      min = ${fparse p0z - noise}
      max = ${fparse p0z + noise}
      seed = 124
    []
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
  [u_z]
    order = FIRST
    family = LAGRANGE
  []
  [Pmag]
    order = FIRST
    family = LAGRANGE
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
[]

[Materials]
  [Landau_P]
    type = GenericConstantMaterial
    prop_names  = 'alpha1     alpha11     alpha12 alpha111 alpha112 alpha123     alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '${alpha1}  ${alpha11}  0.49    6.6      2.9      ${alpha123}  0.0       0.0       0.0       0.0'
  []
  [Landau_G]
    type = GenericConstantMaterial
    prop_names  = 'G110     G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '${G110}  0.51     -0.02    0.02     0.0'
  []
  [mat_C]
    type = GenericConstantMaterial
    prop_names  = 'C11      C12      C44'
    prop_values = '${C11}   ${C12}   ${C44}'
  []
  [mat_Q]
    type = GenericConstantMaterial
    prop_names  = 'Q11   Q12     Q44'
    prop_values = '0.11  -0.043  0.0295'
  []
  [misfit]
    type = GenericConstantRankTwoTensor
    tensor_name = global_strain
    tensor_values = '${um} 0 0   0 ${um} 0   0 0 ${umzz}'
  []
  [strain]
    type = ComputeSmallStrain
    displacements = 'u_x u_y u_z'
    global_strain = global_strain
  []
[]

[Kernels]
  [bed_x]
    type = BulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
  []
  [walled_x]
    type = WallEnergyDerivative
    variable = polar_x
    component = 0
  []
  [electrostr_x]
    type = CubicParentElasticPDerivative
    variable = polar_x
    component = 0
  []
  [time_x]
    type = TimeDerivative
    variable = polar_x
  []
  [bed_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
  []
  [walled_y]
    type = WallEnergyDerivative
    variable = polar_y
    component = 1
  []
  [electrostr_y]
    type = CubicParentElasticPDerivative
    variable = polar_y
    component = 1
  []
  [time_y]
    type = TimeDerivative
    variable = polar_y
  []
  [bed_z]
    type = BulkEnergyDerivativeEighth
    variable = polar_z
    component = 2
  []
  [walled_z]
    type = WallEnergyDerivative
    variable = polar_z
    component = 2
  []
  [electrostr_z]
    type = CubicParentElasticPDerivative
    variable = polar_z
    component = 2
  []
  [time_z]
    type = TimeDerivative
    variable = polar_z
  []
[]

[BCs]
  [Periodic]
    [xy]
      auto_direction = 'x y'
      variable = 'polar_x polar_y polar_z'
    []
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
  [Px]
    type = ElementAverageValue
    variable = polar_x
    execute_on = 'initial timestep_end'
  []
  [Py]
    type = ElementAverageValue
    variable = polar_y
    execute_on = 'initial timestep_end'
  []
  [Pz]
    type = ElementAverageValue
    variable = polar_z
    execute_on = 'initial timestep_end'
  []
  [Pmag_max]
    type = ElementExtremeValue
    variable = Pmag
    value_type = max
    execute_on = 'initial timestep_end'
  []
[]

[Executioner]
  type = Transient
  solve_type = LINEAR
  [TimeIntegrator]
    type = ActuallyExplicitEuler
    solve_type = lumped
    use_constant_mass = true
  []
  dt = ${dt_polar}
  end_time = 1e9
[]

[Outputs]
  print_linear_residuals = false
[]
