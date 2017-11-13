[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 8
  ny = 8
  nz = 22
  xmin = -1
  xmax = 1
  ymin = -1
  ymax = 1
  zmin = -4
  zmax = 4
  elem_type = HEX8
[]

[GlobalParams]
  len_scale = 1.0

  ##################################
  ##--Landau/coupling parameters--##
  ##################################

  alpha1 = -3.362
  alpha11 = -2.646
  alpha12 = 3.274
  alpha111 = -0.596
  alpha112 = 0.2634
  alpha123 = -7.132
  alpha1111 = 0.09043
  alpha1112 = -0.2284
  alpha1122 = 0.4636
  alpha1123 = 1.493

  beta1 = -0.03362
  beta11 = 0.00005396
  beta12 = 0.00006314
  beta111 = -0.596
  beta112 = -5.203e-8
  beta123 = -5.91e-8
  beta1111 = 4.89e-11
  beta1112 = 1.598e-11
  beta1122 = 1.194e-10
  beta1123 = 1.002e-10

  t1111 = 0.0172
  t1122 = 0.02273
  t1212 = -0.02844
  t42111111 = 0.002371
  t24111111 = -0.00005689
  t42111122 = -0.009069
  t24112222 = -0.00004608
  t42112233 = -0.008438
  t24112233 = -0.00002421
  t42112211 = 0.006805
  t24111122 = 0.00002594
  t42111212 = -0.01954
  t42123312 = 0.0266
  t24121112 = -8.314e-6
  t24121233 = -0.00005457
  t6211111111 = -0.002573
  t2611111111 = 1.132e-7
  t6211111122 = 0.00137
  t2611222222 = 5.255e-8
  t4411111111 = 0.00004747
  t4411112222 = 0.00001888

  H110 = 0.08
  H11_H110 = 0.6
  H12_H110 = 0
  H44_H110 = 0.3
  H44P_H110 = 0.3

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  potential_int = potential_int

  antiferrodis_A_x = antiferrodis_A_x
  antiferrodis_A_y = antiferrodis_A_y
  antiferrodis_A_z = antiferrodis_A_z
[]



[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-7
      max = 0.1e-7
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-7
      max = 0.1e-7
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-7
      max = 0.1e-7
    [../]
  [../]

  [./potential_int]
    order = FIRST
    family = LAGRANGE
  [../]

  [./antiferrodis_A_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-7
      max = 0.1e-7
    [../]
  [../]
  [./antiferrodis_A_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-7
      max = 0.1e-7
    [../]
  [../]
  [./antiferrodis_A_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-7
      max = 0.1e-7
    [../]
  [../]
[]

[Kernels]

  ### Operators for the polar field: ###
  [./bed_x]
    type = BulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = BulkEnergyDerivativeEighth
    variable = polar_z
    component = 2
  [../]

  [./rpc_P_x]
    type = RotoPolarCoupledEnergyPolarDerivative
    variable = polar_x
    component = 0
  [../]
  [./rpc_P_y]
    type = RotoPolarCoupledEnergyPolarDerivative
    variable = polar_y
    component = 1
  [../]
  [./rpc_P_z]
    type = RotoPolarCoupledEnergyPolarDerivative
    variable = polar_z
    component = 2
  [../]

  ###Time dependence
  [./polar_x_time]
    type = TimeDerivativeScaled
    variable = polar_x
    time_scale = 1.0
  [../]
  [./polar_y_time]
     type = TimeDerivativeScaled
     variable=polar_y
    time_scale = 1.0
  [../]
  [./polar_z_time]
     type = TimeDerivativeScaled
     variable = polar_z
    time_scale = 1.0
  [../]

  ####Electrostatics
  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
  [../]
  [./FE_E_int]
     type = Electrostatics
     permittivity = 0.08854187
     variable = potential_int
  [../]
  [./polar_electric_px]
     type = PolarElectricPStrong
     variable = polar_x
     component = 0
  [../]
  [./polar_electric_py]
     type = PolarElectricPStrong
     variable = polar_y
     component = 1
  [../]
  [./polar_electric_pz]
     type = PolarElectricPStrong
     variable = polar_z
     component = 2
  [../]

  ### Operators for the antiferrodistortive field: ###

  [./bed_A_x]
    type = RotoBulkEnergyDerivativeEighth
    variable = antiferrodis_A_x
    component = 0
  [../]
  [./bed_A_y]
    type = RotoBulkEnergyDerivativeEighth
    variable = antiferrodis_A_y
    component = 1
  [../]
  [./bed_A_z]
    type = RotoBulkEnergyDerivativeEighth
    variable = antiferrodis_A_z
    component = 2
  [../]

  [./rpc_A_x]
    type = RotoPolarCoupledEnergyDistortDerivative
    variable = antiferrodis_A_x
    component = 0
  [../]
  [./rpc_A_y]
    type = RotoPolarCoupledEnergyDistortDerivative
    variable = antiferrodis_A_y
    component = 1
  [../]
  [./rpc_A_z]
    type = RotoPolarCoupledEnergyDistortDerivative
    variable = antiferrodis_A_z
    component = 2
  [../]

  [./Aant_x]
    type = AFDAntiphaseEnergyDerivative
    variable = antiferrodis_A_x
    component = 0
  [../]
  [./Aant_y]
    type = AFDAntiphaseEnergyDerivative
    variable = antiferrodis_A_y
    component = 1
  [../]
  [./Aant_z]
    type = AFDAntiphaseEnergyDerivative
    variable = antiferrodis_A_z
    component = 2
  [../]

  ###Time dependence
  [./antiferrodis_A_x_time]
    type = TimeDerivativeScaled
    variable = antiferrodis_A_x
    time_scale = 1.0
  [../]
  [./antiferrodis_A_y_time]
     type = TimeDerivativeScaled
     variable = antiferrodis_A_y
    time_scale = 1.0
  [../]
  [./antiferrodis_A_z_time]
     type = TimeDerivativeScaled
     variable = antiferrodis_A_z
    time_scale = 1.0
  [../]

[]

[BCs]
  [./potential_cube5]
    type = DirichletBC
    boundary = 'front'
    value = 0.0
    variable = potential_int
  [../]
  [./potential_cube6]
    type = DirichletBC
    boundary = 'back'
    value = 0.01
    variable = potential_int
  [../]

  [./Periodic]
    [./z]
      auto_direction = 'z'
      variable = 'antiferrodis_A_x antiferrodis_A_y antiferrodis_A_z'
    [../]
  [../]
[]

[Debug]
  show_var_residual_norms = false
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '     121              1e-10    1e-8      1e-6          bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    optimal_iterations = 4
    growth_factor = 1.4
    linear_iteration_ratio = 100
    cutback_factor =  0.55
  [../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 0.05
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_BFO_RPcoupled_test
    elemental_as_nodal = true
  [../]
[]
