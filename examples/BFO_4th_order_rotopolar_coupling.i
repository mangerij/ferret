[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 9
  ny = 9
  nz = 25
  xmin = -3
  xmax = 3
  ymin = -3
  ymax = 3
  zmin = -10
  zmax = 10
  elem_type = HEX8
[]

[GlobalParams]
  len_scale = 1.0

  ##################################
  ##--Landau/coupling parameters--##
  ##################################

  alpha1 = -0.0805
  alpha11 = -0.0522
  alpha12 = 0.0687
  alpha111 = 0.0
  alpha112 = 0.0
  alpha123 = 0.0

  beta1 = -0.0805
  beta11 = -0.0522
  beta12 = 0.0687
  beta111 = 0.0
  beta112 = 0.0
  beta123 = 0.0

  t11 = -0.26
  t12 = -0.25
  t44 = 0.05

  H110 = 0.253
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
      min = -0.1e-6
      max = 0.1e-6
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-6
      max = 0.1e-6
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-6
      max = 0.1e-6
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
      min = -0.1e-6
      max = 0.1e-6
    [../]
  [../]
  [./antiferrodis_A_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-6
      max = 0.1e-6
    [../]
  [../]
  [./antiferrodis_A_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-6
      max = 0.1e-6
    [../]
  [../]
[]

[Kernels]

  ### Operators for the polar field: ###
  [./bed_x]
    type = BulkEnergyDerivativeSixth
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = BulkEnergyDerivativeSixth
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = BulkEnergyDerivativeSixth
    variable = polar_z
    component = 2
  [../]

  [./rpc_P_x]
    type = RotopolarCoupledPolarDerivativeFourth
    variable = polar_x
    component = 0
  [../]
  [./rpc_P_y]
    type = RotopolarCoupledPolarDerivativeFourth
    variable = polar_y
    component = 1
  [../]
  [./rpc_P_z]
    type = RotopolarCoupledPolarDerivativeFourth
    variable = polar_z
    component = 2
  [../]

  ###Time dependence
  [./polar_x_time]
     type=TimeDerivativeScaled
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
    type = BulkAntiferrodistortEnergyDerivativeSixth
    variable = antiferrodis_A_x
    component = 0
  [../]
  [./bed_A_y]
    type = BulkAntiferrodistortEnergyDerivativeSixth
    variable = antiferrodis_A_y
    component = 1
  [../]
  [./bed_A_z]
    type = BulkAntiferrodistortEnergyDerivativeSixth
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


  [./rpc_A_x]
    type = RotopolarCoupledDistortDerivativeFourth
    variable = antiferrodis_A_x
    component = 0
  [../]
  [./rpc_A_y]
    type = RotopolarCoupledDistortDerivativeFourth
    variable = antiferrodis_A_y
    component = 1
  [../]
  [./rpc_A_z]
    type = RotopolarCoupledDistortDerivativeFourth
    variable = antiferrodis_A_z
    component = 2
  [../]

  [./A_x_time]
     type=TimeDerivativeScaled
     variable = antiferrodis_A_x
    time_scale = 1.0
  [../]
  [./A_y_time]
     type = TimeDerivativeScaled
     variable = antiferrodis_A_y
    time_scale = 1.0
  [../]
  [./A_z_time]
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
    value = 0.0001
    variable = potential_int
  [../]
  [./Periodic]
    [./z]
      auto_direction = 'z'
      variable = 'antiferrodis_A_x antiferrodis_A_y antiferrodis_A_z'
    [../]
  [../]

  ################################
  ## - - - - BC Comment - - - - ##
  ################################
  # No electric field means that #
  # there should be P || <111>   #
  ################################
[]

[Postprocessors]
  [./bulk_energy]
   type = BulkEnergy
   execute_on = 'initial timestep_end'
  [../]
  [./electrostatic_energy]
   type = ElectrostaticEnergy
   execute_on = 'initial timestep_end'
  [../]
  [./bulk_antiferrodistort_energy]
   type = BulkAntiferrodistortEnergy
   execute_on = 'initial timestep_end'
  [../]
  [./rotopolar_couple_energy]
   type = RotopolarCouplingEnergy
   execute_on = 'initial timestep_end'
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
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_type -snes_rtol -ksp_rtol -ksp_atol -pc_type  -pc_factor_shift_type -pc_factor_shift_amount'
    petsc_options_value = '     121              1e-10   newtonls    1e-8      1e-8    1e-8            lu            NONZERO                1e-10'
  [../]
[]

[Executioner]
  type = Transient
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.001
    optimal_iterations = 4
    growth_factor = 1.4
    linear_iteration_ratio = 100
    cutback_factor =  0.55
  [../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 0.1
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_BFO_cycloid_test
    elemental_as_nodal = true
  [../]
[]
