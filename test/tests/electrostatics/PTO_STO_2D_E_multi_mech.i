Lx     = 6.0
t_pto  = 2.0
t_sto  = 2.0
Ly     = ${fparse 2.0*t_pto + t_sto}
nx     = 24
ny     = 24
y0     = ${t_pto}
y1     = ${fparse t_pto + t_sto}

a_pto = 3.957
a_sto = 3.905
a_sub = 3.944
um_pto = ${fparse (a_sub - a_pto)/a_pto}
um_sto = ${fparse (a_sub - a_sto)/a_sto}

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
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -ksp_type -ksp_rtol -snes_type'
    petsc_options_value = 'bjacobi  ilu          cg        1e-12     ksponly'
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
