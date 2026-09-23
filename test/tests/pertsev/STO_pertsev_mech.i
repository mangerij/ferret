um = 0
nx = 1
ny = 1
nz = 1
C11 = 336
C12 = 107
C44 = 127
Q11 = 0.0457466385977
Q12 = -0.0134813276811
Q44 = 0.00957174265355
R11 = 8.7e-06
R12 = -7.8e-06
R44 = -9.2e-06

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${ny}
    nz = ${nz}
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    zmin = 0
    zmax = 1
    elem_type = HEX8
  []
[]

[GlobalParams]
  displacements = 'u_x u_y u_z'
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  antiphase_A_x = antiphase_A_x
  antiphase_A_y = antiphase_A_y
  antiphase_A_z = antiphase_A_z
[]

[Variables]
  [u_x] [] [u_y] [] [u_z] []
[]

[AuxVariables]

  [f_el]
    order = CONSTANT
    family = MONOMIAL
  []
  [polar_x] [] [polar_y] [] [polar_z] []
  [antiphase_A_x] [] [antiphase_A_y] [] [antiphase_A_z] []
[]

[Materials]
  [mat_C]
    type = GenericConstantMaterial
    prop_names  = 'C11 C12 C44'
    prop_values = '${C11} ${C12} ${C44}'
  []
  [mat_Q]
    type = GenericConstantMaterial
    prop_names  = 'Q11 Q12 Q44'
    prop_values = '${Q11} ${Q12} ${Q44}'
  []
  [mat_R]
    type = GenericConstantMaterial
    prop_names  = 'R11 R12 R44'
    prop_values = '${R11} ${R12} ${R44}'
  []
  [misfit]
    type = GenericConstantRankTwoTensor
    tensor_name = global_strain
    tensor_values = '${um} 0 0  0 ${um} 0  0 0 0'
  []
  [elasticity_tensor]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '${C11} ${C12} ${C12} ${C11} ${C12} ${C11} ${C44} ${C44} ${C44}'
  []
  [ferro]
    type = ComputeCubicParentElectrostrictiveStrain
    eigenstrain_name = ferro
  []
  [roto]
    type = ComputeSpontaneousRotostrictiveStrain
    eigenstrain_name = roto
  []
  [strain]
    type = ComputeSmallStrain
    displacements = 'u_x u_y u_z'
    eigenstrain_names = 'ferro roto'
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
[]

[BCs]
  [Periodic]
    [xy]
      auto_direction = 'x y'
      variable = 'u_x u_y u_z'
    []
  []
  [fix_x]
    type = DirichletBC
    variable = u_x
    boundary = back
    value = 0
  []
  [fix_y]
    type = DirichletBC
    variable = u_y
    boundary = back
    value = 0
  []
  [fix_z]
    type = DirichletBC
    variable = u_z
    boundary = back
    value = 0
  []
[]

[AuxVariables]
  [exx_a]
    order = CONSTANT
    family = MONOMIAL
  []
  [ezz_a]
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
  [kexx]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = exx_a
    index_i = 0
    index_j = 0
  []
  [kezz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = ezz_a
    index_i = 2
    index_j = 2
  []
[]

[Postprocessors]
  [Felastic_true]
    type = ElementIntegralVariablePostprocessor
    variable = f_el
    execute_on = 'initial timestep_end'
  []
  [exx]
    type = ElementAverageValue
    variable = exx_a
  []
  [ezz]
    type = ElementAverageValue
    variable = ezz_a
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -snes_type'
    petsc_options_value = 'lu       ksponly'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-14
  dt = 1.0
  end_time = 1e9
[]

[Outputs]
  console = false
[]
