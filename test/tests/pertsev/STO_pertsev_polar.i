um = 0
nx = 1
ny = 1
nz = 1
dt_polar = 0.1
alpha1 = 0.0393222495738
alpha11 = 1.69629782974
alpha12 = 4.44268479218
beta1 = -2.42984938579e-06
beta11 = 1.69e-07
beta12 = 3.88e-07
t1111 = 0.000299792012938
t1122 = -2.28639273072e-05
t1212 = -0.000629824117823
C11 = 336
C12 = 107
C44 = 127
Q11 = 0.0457466385977
Q12 = -0.0134813276811
Q44 = 0.00957174265355
R11 = 8.7e-06
R12 = -7.8e-06
R44 = -9.2e-06
p0x = 0
p0y = 0
p0z = 0.15
a0x = 0
a0y = 0
a0z = 7
ts_A = 1e-4

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
  [polar_x]
    [InitialCondition]
      type = ConstantIC
      value = ${p0x}
    []
  []
  [polar_y]
    [InitialCondition]
      type = ConstantIC
      value = ${p0y}
    []
  []
  [polar_z]
    [InitialCondition]
      type = ConstantIC
      value = ${p0z}
    []
  []
  [antiphase_A_x]
    [InitialCondition]
      type = ConstantIC
      value = ${a0x}
    []
  []
  [antiphase_A_y]
    [InitialCondition]
      type = ConstantIC
      value = ${a0y}
    []
  []
  [antiphase_A_z]
    [InitialCondition]
      type = ConstantIC
      value = ${a0z}
    []
  []
[]

[AuxVariables]
  [u_x] [] [u_y] [] [u_z] []
[]

[Materials]
  [Landau_P]
    type = GenericConstantMaterial
    prop_names  = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '${alpha1} ${alpha11} ${alpha12} 0 0 0 0 0 0 0'
  []
  [Landau_A]
    type = GenericConstantMaterial
    prop_names  = 'beta1 beta11 beta12 beta111 beta112 beta123 beta1111 beta1112 beta1122 beta1123'
    prop_values = '${beta1} ${beta11} ${beta12} 0 0 0 0 0 0 0'
  []
  [P_A_couple]
    type = GenericConstantMaterial
    prop_names  = 't1111 t1122 t1212 t42111111 t24111111 t42111122 t24112222 t42112233 t24112233 t42112211 t24111122 t42111212 t42123312 t24121112 t24121233 t6211111111 t2611111111 t6211111122 t2611222222 t4411111111 t4411112222'
    prop_values = '${t1111} ${t1122} ${t1212} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0'
  []
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
  [electrostr_x]
    type = CubicParentElasticPDerivative
    variable = polar_x
    component = 0
  []
  [electrostr_y]
    type = CubicParentElasticPDerivative
    variable = polar_y
    component = 1
  []
  [electrostr_z]
    type = CubicParentElasticPDerivative
    variable = polar_z
    component = 2
  []
  [roto_polar_coupled_x]
    type = RotoPolarCoupledEnergyPolarDerivativeAlt
    variable = polar_x
    component = 0
  []
  [roto_polar_coupled_y]
    type = RotoPolarCoupledEnergyPolarDerivativeAlt
    variable = polar_y
    component = 1
  []
  [roto_polar_coupled_z]
    type = RotoPolarCoupledEnergyPolarDerivativeAlt
    variable = polar_z
    component = 2
  []

  [rbed_x]
    type = RotoBulkEnergyDerivativeEighthAlt
    variable = antiphase_A_x
    component = 0
  []
  [rbed_y]
    type = RotoBulkEnergyDerivativeEighthAlt
    variable = antiphase_A_y
    component = 1
  []
  [rbed_z]
    type = RotoBulkEnergyDerivativeEighthAlt
    variable = antiphase_A_z
    component = 2
  []
  [rotostr_x]
    type = CubicParentElasticADerivative
    variable = antiphase_A_x
    component = 0
  []
  [rotostr_y]
    type = CubicParentElasticADerivative
    variable = antiphase_A_y
    component = 1
  []
  [rotostr_z]
    type = CubicParentElasticADerivative
    variable = antiphase_A_z
    component = 2
  []
  [roto_dis_coupled_x]
    type = RotoPolarCoupledEnergyDistortDerivativeAlt
    variable = antiphase_A_x
    component = 0
  []
  [roto_dis_coupled_y]
    type = RotoPolarCoupledEnergyDistortDerivativeAlt
    variable = antiphase_A_y
    component = 1
  []
  [roto_dis_coupled_z]
    type = RotoPolarCoupledEnergyDistortDerivativeAlt
    variable = antiphase_A_z
    component = 2
  []

  [time_px]
    type = TimeDerivative
    variable = polar_x
  []
  [time_py]
    type = TimeDerivative
    variable = polar_y
  []
  [time_pz]
    type = TimeDerivative
    variable = polar_z
  []
  [time_ax]
    type = TimeDerivativeScaled
    variable = antiphase_A_x
    time_scale = ${ts_A}
  []
  [time_ay]
    type = TimeDerivativeScaled
    variable = antiphase_A_y
    time_scale = ${ts_A}
  []
  [time_az]
    type = TimeDerivativeScaled
    variable = antiphase_A_z
    time_scale = ${ts_A}
  []
[]

[Postprocessors]
  [Px]
    type = ElementAverageValue
    variable = polar_x
  []
  [Py]
    type = ElementAverageValue
    variable = polar_y
  []
  [Pz]
    type = ElementAverageValue
    variable = polar_z
  []
  [Ax]
    type = ElementAverageValue
    variable = antiphase_A_x
  []
  [Ay]
    type = ElementAverageValue
    variable = antiphase_A_y
  []
  [Az]
    type = ElementAverageValue
    variable = antiphase_A_z
  []

  [dPx]
    type = AverageVariableChange
    variable = polar_x
    change_over = time_step
    norm = L1
  []
  [dPy]
    type = AverageVariableChange
    variable = polar_y
    change_over = time_step
    norm = L1
  []
  [dPz]
    type = AverageVariableChange
    variable = polar_z
    change_over = time_step
    norm = L1
  []
  [dAx]
    type = AverageVariableChange
    variable = antiphase_A_x
    change_over = time_step
    norm = L1
  []
  [dAy]
    type = AverageVariableChange
    variable = antiphase_A_y
    change_over = time_step
    norm = L1
  []
  [dAz]
    type = AverageVariableChange
    variable = antiphase_A_z
    change_over = time_step
    norm = L1
  []
  [drift]
    type = LinearCombinationPostprocessor
    pp_names = 'dPx dPy dPz dAx dAy dAz'
    pp_coefs = '${fparse 10.0/dt_polar} ${fparse 10.0/dt_polar} ${fparse 10.0/dt_polar} ${fparse 0.1/dt_polar} ${fparse 0.1/dt_polar} ${fparse 0.1/dt_polar}'
  []
  [Fbulk]
    type = BulkEnergyEighth
  []
  [Froto]
    type = RotoBulkEnergyEighth
  []
  [Fcouple]
    type = RotoPolarCoupledEnergyEighth
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
  console = false
[]
