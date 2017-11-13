[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 2
  ny = 2
  nz = 1
  xmin = -0.5
  xmax = 0.5
  ymin = -1
  ymax = 1
  zmin = -0.1
  zmax = 0.1
  elem_type = HEX8
[]

[GlobalParams]
  len_scale = 1.0

  potential_int = potential_int

  antiferromag_L_x = antiferromag_L_x
  antiferromag_L_y = antiferromag_L_y
  antiferromag_L_z = antiferromag_L_z

  epsilon = 1e-7
  lambda = lambda

  M0 = 1.0

  alphaLL = 0.08

  nx = 1.0
  ny = 1.0
  nz = 0.0

  Ku = 2.5e-1

  A = 4.0e-4

  phi = phi
  theta = theta
[]



[Variables]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
  [../]

  [./antiferromag_L_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      component = 0
    [../]
  [../]
  [./antiferromag_L_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      component = 1
    [../]
  [../]
  [./antiferromag_L_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      component = 2
    [../]
  [../]
  [./lambda]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./phi]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -10.0
      max = 10.0
    [../]
  [../]
  [./theta]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -10.0
      max = 10.0
    [../]
  [../]
[]

[Kernels]
  [./antiferromag_L_x_time]
    type = TimeDerivativeScaled
    variable = antiferromag_L_x
    time_scale = 1.0
  [../]
  [./antiferromag_L_y_time]
    type = TimeDerivativeScaled
    variable = antiferromag_L_y
    time_scale = 1.0
  [../]
  [./antiferromag_L_z_time]
    type = TimeDerivativeScaled
    variable = antiferromag_L_z
    time_scale = 1.0
  [../]

  [./mag_h]
     type = MagHStrong
     variable = potential_int
  [../]
  [./M_int]
     type = Electrostatics
     variable = potential_int
     permittivity = 1.0
  [../]

  [./mag_x]
     type = MagMStrong
     variable = antiferromag_L_x
     component = 0
  [../]
  [./mag_y]
     type = MagMStrong
     variable = antiferromag_L_y
     component = 1
  [../]
  [./mag_z]
     type = MagMStrong
     variable = antiferromag_L_z
     component = 2
  [../]

  [./mag_exchange_x]
    type = MagneticExchangeDerivative
    variable = antiferromag_L_x
    component = 0
  [../]
  [./mag_exchange_y]
    type = MagneticExchangeDerivative
    variable = antiferromag_L_y
    component = 1
  [../]
  [./mag_exchange_z]
    type = MagneticExchangeDerivative
    variable = antiferromag_L_z
    component = 2
  [../]

  [./mag_aniso_x]
    type = MagneticAnisotropyDerivative
    variable = antiferromag_L_x
    component = 0
  [../]
  [./mag_aniso_y]
    type = MagneticAnisotropyDerivative
    variable = antiferromag_L_y
    component = 1
  [../]
  [./mag_aniso_z]
    type = MagneticAnisotropyDerivative
    variable = antiferromag_L_z
    component = 2
  [../]


  [./d_mag_exchange_x]
    type = DampingMagneticExchangeDerivative
    variable = antiferromag_L_x
    component = 0
  [../]
  [./d_mag_exchange_y]
    type = DampingMagneticExchangeDerivative
    variable = antiferromag_L_y
    component = 1
  [../]
  [./d_mag_exchange_z]
    type = DampingMagneticExchangeDerivative
    variable = antiferromag_L_z
    component = 2
  [../]

  [./d_mag_aniso_x]
    type = DampingMagneticAnisotropyDerivative
    variable = antiferromag_L_x
    component = 0
  [../]
  [./d_mag_aniso_y]
    type = DampingMagneticAnisotropyDerivative
    variable = antiferromag_L_y
    component = 1
  [../]
  [./d_mag_aniso_z]
    type = DampingMagneticAnisotropyDerivative
    variable = antiferromag_L_z
    component = 2
  [../]


  [./const_x]
    type = LagrangianMultiplierAntiferromagConstraint
    variable = antiferromag_L_x
    component = 0
  [../]
  [./const_y]
    type = LagrangianMultiplierAntiferromagConstraint
    variable = antiferromag_L_y
    component = 1
  [../]
  [./const_z]
    type = LagrangianMultiplierAntiferromagConstraint
    variable = antiferromag_L_z
    component = 2
  [../]
  [./const_lambda]
    type = LagrangianMultiplierAntiferromagConstraint
    variable = lambda
    component = 3
  [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = ' -snes_check_jacobian'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol  -pc_type -pc_factor_shift_type -pc_factor_shift_amount'
    petsc_options_value = '     121              1e-10      1e-8      1e-6  lu       NONZERO               1e-10 '
  [../]
[]

[Executioner]
  type = Transient
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
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
    file_base = out_mag_test
    elemental_as_nodal = true
  [../]
[]
