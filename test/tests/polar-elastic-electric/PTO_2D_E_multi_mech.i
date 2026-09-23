TC   = 25.0
um   = -0.01
l0   = 1.0
Lx   = 8.0
tf   = 2.0
nx   = 40
ny   = 10

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
[]

[GlobalParams]
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
[]

[Variables]
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
  [polar_x]
    order = FIRST
    family = LAGRANGE
  []
  [polar_y]
    order = FIRST
    family = LAGRANGE
  []
  [polar_z]
    order = FIRST
    family = LAGRANGE
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
  [elasticity_tensor]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '174.269 79.029 79.029 174.269 79.029 174.269 47.62 47.62 47.62'
  []
  [misfit]
    type = GenericConstantRankTwoTensor
    tensor_name = global_strain
    tensor_values = '${um} 0 0   0 0 0   0 0 ${um}'
  []
  [strain]
    type = ComputeSmallStrain
    displacements = 'u_x u_y'
    eigenstrain_names = 'ferro'
    global_strain = global_strain
  []
  [ferro]
    type = ComputeCubicParentElectrostrictiveStrain
    eigenstrain_name = ferro
  []
  [stress]
    type = ComputeLinearElasticStress
  []
[]

[Kernels]
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
[]

[BCs]
  [Periodic]
    [x]
      auto_direction = 'x'
      variable = 'u_x u_y'
    []
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
  [Felastic]
    type = CubicParentElasticEnergy
    execute_on = 'initial timestep_end'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -ksp_type -ksp_rtol -snes_type'
    petsc_options_value = 'bjacobi  ilu          cg        1e-6      ksponly'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  l_max_its = 500
  dt = 1.0
  end_time = 1e9
[]

[Outputs]
  print_linear_residuals = false
[]
