
[Mesh]
  file = exodus_brick.e
[]

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z
 
  potential_H_int = potential_H_int
  potential_H_ext = potential_H_ext

  magnetic_x = magnetic_x
  magnetic_y = magnetic_y
  magnetic_z = magnetic_z

  alpha = 0.85

  Ae = 0.013
  Ms = 0.8
  g0 = 176.1

  permittivity = 1.0
  mu0 = 1256.64
[]

#theta = 0  not a problem!!!

[Variables]
  [./mag_x]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 0.79999
      max = 0.80001
      seed = 3 
    [../]
  [../]
  [./mag_y]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 0.79999
      max = 0.80001
      seed = 3 
    [../]
  [../]
  [./mag_z]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 0.79999
      max = 0.80001
      seed = 3 
    [../]
  [../]
  [./potential_H_int]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
  [../]
  [./potential_H_ext]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
  [../]
[]

[AuxVariables]
 
[]

[AuxKernels]

[]

[Kernels]
  ## Time dependence
  [./mag_x_time]
    type = TimeDerivativeScaled
    variable = mag_x
    time_scale = 1.0
    block = '1'
  [../]
  [./mag_y_time]
    type = TimeDerivativeScaled
    variable = mag_y
    time_scale = 1.0
    block = '1'
  [../]

  [./mag_z_time]
    type = TimeDerivativeScaled
    variable = mag_z
    time_scale = 1.0
    block = '1'
  [../]


   #LLG simple

  # Exchange term

  [./d_llg_exch_x]
    type = ExchangeCartLLG
    variable = mag_x
    component = 0
  [../]
  [./d_llg_exch_y]
    type = ExchangeCartLLG
    variable = mag_y
    component = 1
  [../]
  [./d_llg_exch_z]
    type = ExchangeCartLLG
    variable = mag_z
    component = 2
  [../]

  # Anisotropy term. TURNED OFF

  # Magnetic interaction term

  [./d_HM_x]
    type = InteractionCartLLG
    variable = mag_x
    component = 0
  [../]
  [./d_HM_y]
    type = InteractionCartLLG
    variable = mag_y
    component = 1
  [../]
  [./d_HM_z]
    type = InteractionCartLLG
    variable = mag_z
    component = 2
  [../]

  # Magnetostatic Poisson equation

  [./ext_pot_lap]
    type = Electrostatics
    variable = potential_H_ext
    block = '1 2'
  [../]
  [./int_pot_lap]
    type = Electrostatics
    variable = potential_H_int
    block = '1 2'
  [../]
  [./int_bc_pot_lap]
    type = MagHStrongCart
    variable = potential_H_int
    block = '1'
  [../]
[]

[BCs]
   [./bc_int_pot_R]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = '1'
  [../]

  [./bc_int_pot_L]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = '2'
  [../]

  [./bc_int_pot_T]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = '3'
  [../]
  [./bc_int_pot_Bo]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = '4'
  [../]

  [./bc_int_pot_F]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = '5'
  [../]
  [./bc_int_pot_Ba]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = '6'
  [../]

  [./bc_ext_pot_front]
    type = PresetBC
    variable = potential_H_ext
    value = 0.0
    boundary = '1'
  [../]
  [./bc_ext_pot_back]
    type = PresetBC
    variable = potential_H_ext
    value = 0.0
    boundary = '2'
  [../]
[]

[Postprocessors]
  [./Fexchange]
    type = MagneticExchangeEnergy
    execute_on = 'initial timestep_end'
    block = '1'
    magnetic_x = mag_x
    magnetic_y = mag_y
    magnetic_z = mag_z
  [../]
  [./Fdemag]
    type = MagnetostaticEnergyCart
    execute_on = 'initial timestep_end'
    magnetic_x = mag_x
    magnetic_y = mag_y
    magnetic_z = mag_z
    block = '1'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor 
    pp_names = 'Fexchange Fdemag'
    pp_coefs = ' 1 1 ' 
    execute_on = 'initial timestep_end'
  [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type '
    petsc_options_value = '    121               1e-10      1e-8      1e-8       bjacobi'
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[Executioner]
  type = Transient            
  solve_type = 'NEWTON'
  scheme = 'implicit-euler'   #, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-16
  dtmax = 5.0e-7
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 12
    growth_factor = 1.2
    cutback_factor = 0.4
    dt = 1.0e-7
  [../]
  verbose = true
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = outCartLLG_test_VB
    interval = 1
    elemental_as_nodal = true
  [../]
[]
