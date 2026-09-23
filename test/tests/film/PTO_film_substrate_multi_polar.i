TC   = 25.0
l0   = 1.0
L    = 4.0
dx   = 1.0
nx   = ${fparse int(L/dx + 0.5)}

h_film = 4.0
h_sub  = 2.0
nz     = ${fparse int((h_film + h_sub)/dx + 0.5)}

dt_polar = 0.5

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
[]

[Problem]

  kernel_coverage_check = false
  material_coverage_check = false
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
  [Px_old]
    order = FIRST
    family = LAGRANGE
    block = 0
  []
  [Py_old]
    order = FIRST
    family = LAGRANGE
    block = 0
  []
  [Pz_old]
    order = FIRST
    family = LAGRANGE
    block = 0
  []
  [Pmag]
    order = FIRST
    family = LAGRANGE
    block = 0
  []
[]

[AuxKernels]
  [Px_old]
    type = ParsedAux
    variable = Px_old
    coupled_variables = 'polar_x'
    expression = 'polar_x'
    block = 0
    execute_on = 'initial timestep_begin'
  []
  [Py_old]
    type = ParsedAux
    variable = Py_old
    coupled_variables = 'polar_y'
    expression = 'polar_y'
    block = 0
    execute_on = 'initial timestep_begin'
  []
  [Pz_old]
    type = ParsedAux
    variable = Pz_old
    coupled_variables = 'polar_z'
    expression = 'polar_z'
    block = 0
    execute_on = 'initial timestep_begin'
  []
  [Pmag]
    type = ParsedAux
    variable = Pmag
    coupled_variables = 'polar_x polar_y polar_z'
    expression = 'sqrt(polar_x^2 + polar_y^2 + polar_z^2)'
    block = 0
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
  [strain]
    type = ComputeSmallStrain
    block = 0
  []
[]

[Kernels]
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
    polar_x = Px_old
    polar_y = Py_old
    polar_z = Pz_old
  []
  [bed_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
    block = 0
    polar_x = Px_old
    polar_y = Py_old
    polar_z = Pz_old
  []
  [bed_z]
    type = BulkEnergyDerivativeEighth
    variable = polar_z
    component = 2
    block = 0
    polar_x = Px_old
    polar_y = Py_old
    polar_z = Pz_old
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
    polar_x = Px_old
    polar_y = Py_old
    polar_z = Pz_old
  []
  [electrostr_y]
    type = CubicParentElasticPDerivative
    variable = polar_y
    component = 1
    block = 0
    polar_x = Px_old
    polar_y = Py_old
    polar_z = Pz_old
  []
  [electrostr_z]
    type = CubicParentElasticPDerivative
    variable = polar_z
    component = 2
    block = 0
    polar_x = Px_old
    polar_y = Py_old
    polar_z = Pz_old
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
    block = 0
    execute_on = 'initial timestep_end'
  []
  [Fwall]
    type = WallEnergy
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
    file_base = PTO_film_substrate_multi_exo
  []
[]
