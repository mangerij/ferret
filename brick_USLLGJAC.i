
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 5
  ny = 5
  nz = 5
  xmin = -0.1
  xmax = 0.1
  ymin = -0.1
  ymax = 0.1
  zmin = -0.1
  zmax = 0.1
  elem_type = HEX8
[]

[GlobalParams]
  polar_theta = polar_theta
  azimuth_phi = azimuth_phi
 
  potential_H_int = potential_H_int
  potential_H_ext = potential_H_ext

  magnetic_x = magnetic_x
  magnetic_y = magnetic_y
  magnetic_z = magnetic_z

  alpha = 0.85

  Ae = 0.013
  Ms = 0.8
  g0 = 175.88

  permittivity = 1.0 # this scalar is in the Electrostatics.C kernel. It is grad * permitivitty * grad * potential = ...
  mu0 = 1256.64
[]

[Variables]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    block = '0'
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
    block = '0'
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
    block = '0'
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 6.283185307178
      seed = 3
    [../]
  [../]
  [./potential_H_ext]
    order = FIRST
    family = LAGRANGE
    block = '0'
  [../]
[]

[AuxVariables]
  [./magnetic_x]
    order = FIRST
    family = LAGRANGE
    block = '0'
  [../]
  [./magnetic_y]
    order = FIRST
    family = LAGRANGE
    block = '0'
  [../]
  [./magnetic_z]
    order = FIRST
    family = LAGRANGE
    block = '0'
  [../]

  [./Hexch_x]
    order = CONSTANT
    family = MONOMIAL
    block = '0'
  [../]
  [./Hexch_y]
    order = CONSTANT
    family = MONOMIAL
    block = '0'
  [../]
  [./Hexch_z]
    order = CONSTANT
    family = MONOMIAL
    block = '0'
  [../]

  [./Hdemag_x]
    order = CONSTANT
    family = MONOMIAL
    block = '0'
  [../]
  [./Hdemag_y]
    order = CONSTANT
    family = MONOMIAL
    block = '0'
  [../]
  [./Hdemag_z]
    order = CONSTANT
    family = MONOMIAL
    block = '0'
  [../]
[]

[AuxKernels]
 [./mag_x_c]
   type = MagFieldAux
   component = 0
   variable = magnetic_x
   execute_on = 'initial linear nonlinear timestep_end'
   block = '0'
 [../]
 [./mag_y_c]
   type = MagFieldAux
   component = 1
   variable = magnetic_y
   execute_on = 'initial linear nonlinear timestep_end'
   block = '0'
 [../]
 [./mag_z_c]
   type = MagFieldAux
   component = 2
   variable = magnetic_z
   execute_on = 'initial linear nonlinear timestep_end'
   block = '0'
 [../]

 [./Hexch_x_c]
   type = ExchangeFieldAux
   component = 0
   variable = Hexch_x
   execute_on = 'initial timestep_end'
   block = '0'
 [../]
 [./Hexch_y_c]
   type =  ExchangeFieldAux
   component = 1
   variable = Hexch_y
   execute_on = 'initial timestep_end'
   block = '0'
 [../]
 [./Hexch_z_c]
   type =  ExchangeFieldAux
   component = 2
   variable = Hexch_z
   execute_on = 'initial timestep_end'
   block = '0'
 [../]

 [./Hdemag_x_c]
   type = DemagFieldAux
   component = 0
   variable = Hexch_x
   execute_on = 'initial timestep_end'
   block = '0'
 [../]
 [./Hdemag_y_c]
   type =  DemagFieldAux
   component = 1
   variable = Hexch_y
   execute_on = 'initial timestep_end'
   block = '0'
 [../]
 [./Hdemag_z_c]
   type =  DemagFieldAux
   component = 2
   variable = Hexch_z
   execute_on = 'initial timestep_end'
   block = '0'
 [../]
[]

[Kernels]
  ## Time dependence
  [./polar_time]
    type = TimeDerivativeScaled
    variable = polar_theta
    time_scale = 1.0
    block = '0'
  [../]
  [./azimuthal_time]
    type = TimeDerivativeScaled
    variable = azimuth_phi
    time_scale = 1.0
    block = '0'
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
    block = '0'
  [../]
  [./int_pot_lap]
    type = Electrostatics
    variable = potential_H_int
    block = '0'
  [../]
  [./int_bc_pot_lap]
    type = MagHStrong
    variable = potential_H_int
    block = '0'
  [../]
[]

[BCs]
  ##
  ## NOTE THE JAC DEBUGGER TURNS OFF ALL BCs
  ##
  ##  ##  ##  ##  ##  ##  ##  ##  ##
   [./bc_int_pot_R]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = 'front'
  [../]

  [./bc_int_pot_L]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = 'back'
  [../]

  [./bc_int_pot_T]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = 'left'
  [../]
  [./bc_int_pot_Bo]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = 'right'
  [../]

  [./bc_int_pot_F]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = 'top'
  [../]
  [./bc_int_pot_Ba]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = 'bottom'
  [../]

  [./bc_ext_pot_front]
    type = PresetBC
    variable = potential_H_ext
    value = 0.0
    boundary = 'front'
  [../]
  [./bc_ext_pot_back]
    type = PresetBC
    variable = potential_H_ext
    value = 0.0
    boundary = 'back'
  [../]
[]

[Postprocessors]
  [./Fexchange]
    type = MagneticExchangeEnergy
    execute_on = 'initial timestep_end'
    block = '0'
  [../]
  [./Fdemag]
    type = MagnetostaticEnergy
    execute_on = 'initial timestep_end'
    block = '0'
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
    petsc_options = '-snes_check_jacobian'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type '
    petsc_options_value = '    121               1e-10      1e-8      1e-6       bjacobi'
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[Executioner]
  type = Transient
  #dt = 1.0e-7                
  solve_type = 'NEWTON'
  scheme = 'implicit-euler'   #, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-16
  dtmax = 1.0e-6
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 12
    growth_factor = 1.2
    cutback_factor = 0.85
    dt = 1.0e-8
  [../]
  num_steps = 1
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = outUSLLG_test
    interval = 1
    elemental_as_nodal = true
  [../]
[]
