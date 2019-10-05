
[Mesh]
  file = exodus_brick.e
[]

[GlobalParams]
  polar_theta = polar_theta
  azimuth_phi = azimuth_phi
 
  potential_H_int = potential_H_int
  potential_H_ext = potential_H_ext

  magnetic_x = magnetic_x
  magnetic_y = magnetic_y
  magnetic_z = magnetic_z

  alpha = 1.0 #what is the sign of alpha?? Flipping this makes the exchange energy decrease.

  Ae = 0.013
  Ms = 0.8
  g0 = 175.88 

  permittivity = 1.0 # this scalar is in the Electrostatics.C kernel. It is grad * permitivitty * grad * potential = ...
[]

[Variables]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 3.1459
      seed = 3 
    [../]
  [../]
  [./azimuth_phi]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 6.283185307178
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
  [./magnetic_x]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./magnetic_y]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./magnetic_z]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]

  [./Hexch_x]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  [./Hexch_y]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  [./Hexch_z]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./Hdemag_x]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  [./Hdemag_y]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  [./Hdemag_z]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
[]

[AuxKernels]
 [./mag_x_c]
   type = MagFieldAux
   component = 0
   variable = magnetic_x
   execute_on = 'initial linear nonlinear timestep_end'
   block = '1'
 [../]
 [./mag_y_c]
   type = MagFieldAux
   component = 1
   variable = magnetic_y
   execute_on = 'initial linear nonlinear timestep_end'
   block = '1'
 [../]
 [./mag_z_c]
   type = MagFieldAux
   component = 2
   variable = magnetic_z
   execute_on = 'initial linear nonlinear timestep_end'
   block = '1'
 [../]

 [./Hexch_x_c]
   type = ExchangeFieldAux
   component = 0
   variable = Hexch_x
   execute_on = 'initial timestep_end'
   block = '1'
 [../]
 [./Hexch_y_c]
   type =  ExchangeFieldAux
   component = 1
   variable = Hexch_y
   execute_on = 'initial timestep_end'
   block = '1'
 [../]
 [./Hexch_z_c]
   type =  ExchangeFieldAux
   component = 2
   variable = Hexch_z
   execute_on = 'initial timestep_end'
   block = '1'
 [../]

 [./Hdemag_x_c]
   type = DemagFieldAux
   component = 0
   variable = Hexch_x
   execute_on = 'initial timestep_end'
   block = '1'
 [../]
 [./Hdemag_y_c]
   type =  DemagFieldAux
   component = 1
   variable = Hexch_y
   execute_on = 'initial timestep_end'
   block = '1'
 [../]
 [./Hdemag_z_c]
   type =  DemagFieldAux
   component = 2
   variable = Hexch_z
   execute_on = 'initial timestep_end'
   block = '1'
 [../]
[]

[Kernels]
  ## Time dependence
  [./polar_time]
    type = TimeDerivativeScaled
    variable = polar_theta
    time_scale = 1.0
    block = '1'
  [../]
  [./azimuthal_time]
    type = TimeDerivativeScaled
    variable = azimuth_phi
    time_scale = 1.0
    block = '1'
  [../]

   #LLG simple

  # Exchange term

  [./d_llg_exch_th]
    type = ExchangeUSLLG
    variable = polar_theta
    component = 0
  [../]
  [./d_llg_exch_phi]
    type = ExchangeUSLLG
    variable = azimuth_phi
    component = 1
  [../]

  # Anisotropy term. TURNED OFF

  #[./d_llg_anis_th]
  #  type = AnisotropyUSLLG
  #  variable = polar_theta
  #  component = 0
  #[../]
  #[./d_llg_anis_phi]
  #  type = AnisotropyUSLLG
  #  variable = azimuth_phi
  #  component = 1
  #[../]

  # Magnetic interaction term

  [./d_HM_1]
    type = InteractionUSLLG
    variable = azimuth_phi
    component = 1
  [../]
  [./d_HM_0]
    type = InteractionUSLLG
    variable = polar_theta
    component = 0
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
    type = MagHStrong
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
  [../]
  [./Fdemag]
    type = MagnetostaticEnergy
    execute_on = 'initial timestep_end'
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
  [../]
[]

[Debug]
  show_var_residual_norms = false
[]

[Executioner]
  type = Transient
  #dt = 1.0e-7                
  solve_type = 'NEWTON'
  scheme = 'implicit-euler'   #, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-16
  dtmax = 1.0e-8
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 12
    growth_factor = 1.2
    cutback_factor = 0.85
    dt = 1.0e-10
  [../]
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = outUSLLG_test_0
    interval = 1
    elemental_as_nodal = true
  [../]
[]
