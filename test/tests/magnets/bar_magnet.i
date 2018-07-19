## BENCHMARK PROBLEM #1
## DOES NOT WORK AT THE MOMENT
##
##


[Mesh]
  file = exodus_brick_r.e
[]

[GlobalParams]
  polar_theta = polar_theta
  azimuth_phi = azimuth_phi
 
  potential_H_int = potential_H_int
  potential_H_ext = potential_H_ext


  # in units of picograms, micrometers, microseconds, microcoulombs, and microvolts
  alpha = 100.0

  Ae = 0.1 #should be 0.013
  M = 0.8

  nx = 1.0
  ny = 1.0
  nz = 1.0

  Ku = 0.5

  g0 = 1.0

  #try uniform state and calculate demag field

  # the below "permittivity is just a 
  permittivity = 1.25664
  mu0 = 1.25664
[]

[Variables]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 1.57078
      max = 1.57079632679
      seed = 1
    [../]
    block = '1'
  [../]
  [./azimuth_phi]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 3.14158
      max = 3.14159
      seed = 1
    [../]
    block = '1'
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

  [./demag_x]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./demag_y]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./demag_z]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
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

 [./dmag_x_c]
   type = ElecFieldAux
   component = 0
   variable = demag_x
   execute_on = 'initial timestep_end'
   potential_E_int = potential_H_int
   potential_E_ext = potential_H_ext
 [../]
 [./dmag_y_c]
   type = ElecFieldAux
   component = 1
   variable = demag_y
   execute_on = 'initial timestep_end'
   potential_E_int = potential_H_int
   potential_E_ext = potential_H_ext
   block = '1 2'
 [../]
 [./dmag_z_c]
   type = ElecFieldAux
   component = 2
   variable = demag_z
   execute_on = 'initial timestep_end'
   potential_E_int = potential_H_int
   potential_E_ext = potential_H_ext
   block = '1 2'
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
    type = ConstrainedExchangeLLG
    variable = polar_theta
    component = 0
    block = '1'
  [../]
  [./d_llg_exch_phi]
    type = ConstrainedExchangeLLG
    variable = azimuth_phi
    component = 1
    block = '1'
  [../]
  [./d_llg_aniso_th]
    type = ConstrainedAnisotropyLLG
    variable = polar_theta
    component = 0
  [../]
  [./d_llg_aniso_phi]
    type = ConstrainedAnisotropyLLG
    variable = azimuth_phi
    component = 1
    block = '1'
  [../]

  # Magnetic interaction term

  [./d_HM_1]
    type = ConstrainedInteractionLLG
    variable = azimuth_phi
    component = 1
    block = '1'
  [../]
  [./d_HM_0]
    type = ConstrainedInteractionLLG
    variable = polar_theta
    component = 0
    block = '1'
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
    pp_names = 'Fexchange Faniso'
    pp_coefs = ' 1 1' 
    execute_on = 'initial timestep_end'
    block = '1'
  [../]
[]

[BCs]
  #[./Periodic]
  #  [./xy]
  #    auto_direction = 'x y z'
  #    variable = 'polar_theta azimuth_phi'
  #  [../]
  #[../]

  [./bc_int_pot_front]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '1'
  [../]

  [./bc_int_pot_back]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '2'
  [../]

  [./bc_ext_pot_back]
    type = DirichletBC
    variable = potential_H_ext
    value = 0.0
    boundary = '2'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '     121             1e-10       1e-10      1e-6      lu   '
  [../]
[]

[Executioner]
  type = Transient
  dt = 1e-5               # NOTE FOR MAGNETS, TIMESTEP AND ALPHA CHOICE ARE INTERTWINED
  solve_type = 'NEWTON'   # "PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'         # "implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-10
  dtmax = 1.0e-2
[]

[Debug]
  show_var_residual_norms = false
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_mag_uni_r
    interval = 1
    elemental_as_nodal = true
  [../]
[]
