TC   = 25.0
um   = -0.002
C11 = 174.6032
C12 = 79.3651
C44 = 111.1111
umzz = ${fparse -2.0*C12/C11*um}
alpha1   = ${fparse 3.8e-4*(TC - 479.0)}
alpha11  = -0.073
alpha123 = -3.7
G110 = 1.0
Lx = 2.0
Ly = 2.0
tf = 2.0
nx = 2
ny = 2
nz = 8

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

  [f_el]
    order = CONSTANT
    family = MONOMIAL
  []
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
  [strain_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_xz]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [f_el]
    type = ElasticEnergyAux
    variable = f_el
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
  [strain_xz]
    type = RankTwoAux
    variable = strain_xz
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 2
    execute_on = 'initial timestep_end'
  []
[]

[Materials]
  [Landau_P]
    type = GenericConstantMaterial
    prop_names  = 'alpha1     alpha11     alpha12 alpha111 alpha112 alpha123     alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '${alpha1}  ${alpha11}  0.75    0.26     0.61     ${alpha123}  0.0       0.0       0.0       0.0'
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
    prop_values = '0.089 -0.026  0.03375'
  []
  [misfit]
    type = GenericConstantRankTwoTensor
    tensor_name = global_strain
    tensor_values = '${um} 0 0   0 ${um} 0   0 0 ${umzz}'
  []
  [elasticity_tensor]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '${C11} ${C12} ${C12} ${C11} ${C12} ${C11} ${C44} ${C44} ${C44}'
  []
  [strain]
    type = ComputeSmallStrain
    displacements = 'u_x u_y u_z'
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
    displacements = 'u_x u_y u_z'
    component = 0
  []
  [div_y]
    type = StressDivergenceTensors
    variable = u_y
    displacements = 'u_x u_y u_z'
    component = 1
  []
  [div_z]
    type = StressDivergenceTensors
    variable = u_z
    displacements = 'u_x u_y u_z'
    component = 2
  []
[]

[BCs]
  [Periodic]
    [xy]
      auto_direction = 'x y'
      variable = 'u_x u_y u_z'
    []
  []
  [clamp_x]
    type = DirichletBC
    variable = u_x
    boundary = 'back'
    value = 0
  []
  [clamp_y]
    type = DirichletBC
    variable = u_y
    boundary = 'back'
    value = 0
  []
  [clamp_z]
    type = DirichletBC
    variable = u_z
    boundary = 'back'
    value = 0
  []
[]

[Postprocessors]
  [Felastic_true]
    type = ElementIntegralVariablePostprocessor
    variable = f_el
    execute_on = 'initial timestep_end'
  []
  [Felastic]
    type = CubicParentElasticEnergy
    execute_on = 'initial timestep_end'
  []
  [exx]
    type = ElementAverageValue
    variable = strain_xx
    execute_on = 'initial timestep_end'
  []
  [ezz]
    type = ElementAverageValue
    variable = strain_zz
    execute_on = 'initial timestep_end'
  []
  [exz]
    type = ElementAverageValue
    variable = strain_xz
    execute_on = 'initial timestep_end'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -ksp_type -ksp_rtol -snes_type'
    petsc_options_value = 'bjacobi  ilu          cg        1e-8      ksponly'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-12
  l_max_its = 500
  dt = 1.0
  end_time = 1e9
[]

[Outputs]
  print_linear_residuals = false
[]
