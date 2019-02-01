
[Mesh]
  file = sphere_medium_exodus.e
[]



[GlobalParams]
  polar_theta = polar_theta
  azimuth_phi = azimuth_phi
 
  potential_H_int = potential_H_int
  potential_H_ext = potential_H_ext

  alpha = -10.0  #can make alpha much larger (more negative) and therefore, then the domains start to form but then gets stiff

  Ae = 0.013  #seems should be POSITIVE

  Ms = 0.8

  K1 = 0.0
  K2 = 0.0

  nx = 1.0
  ny = 0.0
  nz = 0.0

  g0 = 175929.0              #sign of gamma for conservative term just which way its spinning, damping => pumping energy in or out

  permittivity = 0.000795775 # This should be 1/mu0... in these units
  mu0 = 1256.64              # mu0
[]

[Variables]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    block = '2'
    [./InitialCondition]
      type = RandomIC
      min = 1.5
      max = 3.14159265357
      seed = 3 
    [../]
  [../]
  [./azimuth_phi]
    order = FIRST
    family = LAGRANGE
    block = '2'
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
    block = '2'
  [../]
  [./magnetic_y]
    order = FIRST
    family = LAGRANGE
    block = '2'
  [../]
  [./magnetic_z]
    order = FIRST
    family = LAGRANGE
    block = '2'
  [../]
  [./bounds_theta_dummy]
    order = FIRST
    family = LAGRANGE
    block = '2'
  [../]
  [./bounds_phi_dummy]
    order = FIRST
    family = LAGRANGE
    block = '2'
  [../]
[]

[AuxKernels]
 [./mag_x_c]
   type = MagFieldAux
   component = 0
   variable = magnetic_x
   execute_on = 'initial timestep_end'
   block = '2'
 [../]
 [./mag_y_c]
   type = MagFieldAux
   component = 1
   variable = magnetic_y
   execute_on = 'initial timestep_end'
   block = '2'
 [../]
 [./mag_z_c]
   type = MagFieldAux
   component = 2
   variable = magnetic_z
   execute_on = 'initial timestep_end'
   block = '2'
 [../]
[]

[Kernels]
  ## Time dependence
  [./polar_time]
    type = TimeDerivativeScaled
    variable = polar_theta  
    time_scale = 1.0
    block = '2'
  [../]
  [./azimuthal_time]
    type = TimeDerivativeScaled
    variable = azimuth_phi
    time_scale = 1.0
    block = '2'
  [../]

  # LLG simple
  [./d_llg_exch_th]
    type = ExchangeUSLLG
    variable = polar_theta
    component = 0
    block = '2'
  [../]
  [./d_llg_exch_phi]
    type = ExchangeUSLLG
    variable = azimuth_phi
    component = 1
    block = '2'
  [../]

  [./d_llg_aniso_th]
    type = AnisotropyUSLLG
    variable = polar_theta
    component = 0
    block = '2'
  [../]
  [./d_llg_aniso_phi]
    type = AnisotropyUSLLG
    variable = azimuth_phi
    component = 1
    block = '2'
  [../]

  # Magnetic interaction term

  [./d_HM_0]
    type = InteractionUSLLG 
    variable = polar_theta
    component = 0
    block = '2'
  [../]
  [./d_HM_1]
    type = InteractionUSLLG
    variable = azimuth_phi
    component = 1
    block = '2'
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
    block = '2'
  [../]
[]

[Postprocessors]
  [./Fexchange]
    type = MagneticConstrainedExchangeEnergy
    execute_on = 'initial timestep_end'
    block = '2'
  [../]
  [./Faniso]
    type = MagneticConstrainedAltAnisotropyEnergy
    execute_on = 'initial timestep_end'
    block = '2'
  [../]
  [./Fmag]
    type = MagnetostaticEnergy
    execute_on = 'initial timestep_end'
    block = '2'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor 
    pp_names = 'Fexchange Faniso Fmag'
    pp_coefs = ' 1 1 1' 
    execute_on = 'initial timestep_end'
  [../]
[]

[BCs]
   [./bc_int_pot_1]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '1'
  [../]
   [./bc_int_pot_2]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '2'
  [../]
   [./bc_int_pot_3]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '3'
  [../]
   [./bc_int_pot_4]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '4'
  [../]
   [./bc_int_pot_5]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '5'
  [../]
   [./bc_int_pot_6]
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
    upper = 3.14158                  #should this just be \pi?
    lower = 0.001
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
    petsc_options = '-snes_ksp_ew -snes_converged_reason'          #-snes_monitor -snes_linesearch_monitor' 
                                                                    #-ksp_monitor_true_residual -ksp_monitor_singular_value                                                                                 
    petsc_options_iname = '-snes_type   -snes_rtol -ksp_rtol  -snes_atol -snes_stol -pc_type  -pc_factor_mat_solver_package -pc_factor_shift_type   -pc_factor_shift_amount'
    petsc_options_value = ' vinewtonrsls     1e-8      1e-6     1e-10      0          lu         superlu_dist                  NONZERO                          1e-10'
  [../]
[]

[Debug]
  show_var_residual_norms = false
[]

[Executioner]
  type = Transient
                               # NOTE FOR MAGNETS, TIMESTEP AND ALPHA CHOICE ARE INTERTWINED
  solve_type = 'PJFNK'         # PJFNK not supported for rsls (or is it?)! 
  dtmin = 1e-16
  dtmax = 1e-6               #alpha large might make dt too small not allowed!

  scheme = 'implicit-euler' 

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 30
    growth_factor = 1.2
    cutback_factor = 0.5
    dt = 1e-6
  [../]
  num_steps = 2
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_USLLG_test
    interval = 1
    elemental_as_nodal = true
  [../]
[]
