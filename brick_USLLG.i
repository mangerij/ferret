
[Mesh]
  file = exodus_Lbrick.e
[]



[GlobalParams]
  polar_theta = polar_theta
  azimuth_phi = azimuth_phi
 
  potential_H_int = potential_H_int
  potential_H_ext = potential_H_ext

  alpha = -10000.0

  Ae = -0.0000129313 #can multiply by 1/mu0 and then hit alpha with 100 and get same structure faster (10 time steps of 1e-7)

  Ms = 0.8
  M = 0.8

  K1 = 0.0
  K2 = 0.0


  nx = 1.0
  ny = 0.0
  nz = 0.0

  Ku = 0.0

  g0 = 175929.0 #sign of gamma for conservative term just which way its spinning, damping term => pumping energy in or out

  permittivity = 0.000795775 #0.000795775 This should be 1/mu0... in these units
  mu0 = 1.0
[]

[Variables]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 1.9
      max = 3.14159265357
      seed = 3 
    [../]
  [../]
  [./azimuth_phi]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 6.1
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
  [./bounds_theta_dummy]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./bounds_phi_dummy]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
[]

[AuxKernels]
 [./mag_x_c]
   type = MagFieldAux
   component = 0
   variable = magnetic_x
   execute_on = 'initial timestep_end'
   block = '1'
 [../]
 [./mag_y_c]
   type = MagFieldAux
   component = 1
   variable = magnetic_y
   execute_on = 'initial timestep_end'
   block = '1'
 [../]
 [./mag_z_c]
   type = MagFieldAux
   component = 2
   variable = magnetic_z
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

  # LLG simple
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

  [./d_llg_aniso_th]
    type = AnisotropyUSLLG
    variable = polar_theta
    component = 0
  [../]
  [./d_llg_aniso_phi]
    type = AnisotropyUSLLG
    variable = azimuth_phi
    component = 1
  [../]

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

  # Magnetostatic equations

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

[Postprocessors]
  [./Fexchange]
    type = MagneticConstrainedExchangeEnergy
    execute_on = 'initial timestep_end'
    block = '1'
  [../]
  [./Faniso]
    type = MagneticConstrainedAltAnisotropyEnergy
    execute_on = 'initial timestep_end'
    block = '1'
  [../]
  [./Fmag]
    type = MagnetostaticEnergy
    execute_on = 'initial timestep_end'
    block = '1'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor 
    pp_names = 'Fexchange Faniso Fmag'
    pp_coefs = ' 1 1 1' 
    execute_on = 'initial timestep_end'
  [../]
[]

[BCs]
   [./bc_int_pot_R]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '1'
  [../]

  [./bc_int_pot_L]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '2'
  [../]

  [./bc_int_pot_T]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '3'
  [../]
  [./bc_int_pot_Bo]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '4'
  [../]

  [./bc_int_pot_F]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '5'
  [../]
  [./bc_int_pot_Ba]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '6'
  [../]


  [./bc_ext_pot_front]
    type = DirichletBC
    variable = potential_H_ext
    value = 0.0
    boundary = '1'
  [../]
  [./bc_ext_pot_back]
    type = DirichletBC
    variable = potential_H_ext
    value = 0.0
    boundary = '2'
  [../]
[]

[Bounds]
  [./u_bounds]
    type = BoundsAux
    variable = bounds_theta_dummy
    bounded_variable = polar_theta
    upper = 3.14159265357                  #should this just be \pi?
    lower = -3.14159265357
    execute_on = 'Initial Linear Nonlinear'
  [../]

  [./v_bounds]
    type = BoundsAux
    variable = bounds_phi_dummy
    bounded_variable = azimuth_phi
    upper = 6.283185307177
    lower = 0.0
    execute_on = 'Initial Linear Nonlinear'
  [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason'
    petsc_options_iname = '-snes_type  -snes_mf_operator -ksp_gmres_restart -snes_stol -snes_atol  -snes_rtol -ksp_rtol -pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_shift_amount'
    petsc_options_value = ' vinewtonrsls     0                 121          0   1e-9        1e-8      1e-12      asm        8              lu         nonzero                1e-10'
  [../]
[]

# -snes_mf_operator 0???

[Debug]
  show_var_residual_norms = false
[]

[Executioner]
  type = Transient
  #dt = 1.0e-7                 # NOTE FOR MAGNETS, TIMESTEP AND ALPHA CHOICE ARE INTERTWINED
  solve_type = 'NEWTON'       # PJFNK not supported for rsls (or is it)! 
  scheme = 'bdf2'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-16
  dtmax = 1.5e-6
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 30
    growth_factor = 1.2
    cutback_factor = 0.5
    dt = 1e-9
  [../]
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = USLLG_test3G
    interval = 1
    elemental_as_nodal = true
  [../]
[]
