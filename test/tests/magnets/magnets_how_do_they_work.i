
[Mesh]
  # Mesh details
  type = GeneratedMesh
  dim = 3
  nx = 8
  ny = 8
  nz = 2
  elem_type = HEX8


  # dimensions in units of microns
  xmin = -0.02
  xmax = 0.02
  ymin = -0.02
  ymax = 0.02
  zmin = -0.005
  zmax = 0.005
[]

[GlobalParams]
  polar_theta = polar_theta
  azimuth_phi = azimuth_phi
 
  potential_H_int = potential_H_int
  potential_H_ext = potential_H_ext

  alpha = 1.01

  Ae = -0.13
  M = 0.8

  nx = 0.0
  ny = 1.0
  nz = 0.0

  Ku = 0.5

  g0 = 1.0

  permittivity = 1.25
  mu0 = 1.25
[]

[Variables]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.00000001
      max = 3.14159265357
      seed = 5
    [../]
  [../]
  [./azimuth_phi]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.000000001
      max = 6.283185307178
      seed = 5
    [../]
  [../]
  [./potential_H_int]
    order = FIRST
    family = LAGRANGE
  [../]
  [./potential_H_ext]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./magnetic_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./magnetic_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./magnetic_z]
    order = FIRST
    family = LAGRANGE
  [../]
  [./bounds_dummy]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxKernels]
 [./mag_x_c]
   type = MagFieldAux
   component = 0
   variable = magnetic_x
   execute_on = 'initial timestep_end'
 [../]
 [./mag_y_c]
   type = MagFieldAux
   component = 1
   variable = magnetic_y
   execute_on = 'initial timestep_end'
 [../]
 [./mag_z_c]
   type = MagFieldAux
   component = 2
   variable = magnetic_z
   execute_on = 'initial timestep_end'
 [../]
[]

[Kernels]
  ## Time dependence
  [./polar_time]
    type = TimeDerivativeScaled
    variable = polar_theta
    time_scale = 1.0
  [../]
  [./azimuthal_time]
    type = TimeDerivativeScaled
    variable = azimuth_phi
    time_scale = 1.0
  [../]

  # LLG simple
  [./d_llg_exch_th]
    type = ConstrainedExchangeLLG
    variable = polar_theta
    component = 0
  [../]
  [./d_llg_exch_phi]
    type = ConstrainedExchangeLLG
    variable = azimuth_phi
    component = 1
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
  [../]

  # Magnetic interaction term

  [./d_HM_1]
    type = ConstrainedInteractionLLG
    variable = azimuth_phi
    component = 1
  [../]
  [./d_HM_0]
    type = ConstrainedInteractionLLG
    variable = polar_theta
    component = 0
  [../]

  # Magnetostatic equations

  [./ext_pot_lap]
    type = Electrostatics
    variable = potential_H_ext
  [../]
  [./int_pot_lap]
    type = Electrostatics
    variable = potential_H_int
  [../]
  [./int_bc_pot_lap]
    type = MagHStrong
    variable = potential_H_int
  [../]
[]

[Postprocessors]
  [./Fexchange]
    type = MagneticConstrainedExchangeEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./Faniso]
    type = MagneticConstrainedAltAnisotropyEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./Fmag]
    type = MagnetostaticEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor 
    pp_names = 'Fexchange Faniso'
    pp_coefs = ' 1 1' 
    execute_on = 'initial timestep_end'
  [../]
[]

[BCs]
  [./bc_int_pot_R]
    type = DirichletBC
    variable = potential_H_int
    value = 500.5
    boundary = 'right'
  [../]
  [./bc_int_pot_L]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = 'left'
  [../]

  [./bc_ext_pot_front]
    type = DirichletBC
    variable = potential_H_ext
    value = 0.0
    boundary = 'front'
  [../]
  [./bc_ext_pot_back]
    type = DirichletBC
    variable = potential_H_ext
    value = 0.0
    boundary = 'back'
  [../]
[]

[Bounds]
  [./u_bounds]
    type = BoundsAux
    variable = bounds_dummy
    bounded_variable = polar_theta
    upper = 3.14159265357
    lower = 0.0
    execute_on = 'Initial Linear Nonlinear'
  [../]

  [./v_bounds]
    type = BoundsAux
    variable = bounds_dummy
    bounded_variable = azimuth_phi
    upper = 6.283185307178
    lower = 0.0
    execute_on = 'Initial Linear Nonlinear'
  [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-snes_type  -snes_mf_operator -ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -pc_factor_shift_type -pc_factor_shift_amount'
    petsc_options_value = ' vinewtonrsls     0                 121            1e-10        1e-8      1e-8                lu     nonzero                  1e-10'
  [../]
  #vinewtonssls
[]

[Debug]
  show_var_residual_norms = false
[]

[Executioner]
  type = Transient
  dt = 1.0e-5                 #NOTE FOR MAGNETS, TIMESTEP AND ALPHA CHOICE ARE INTERTWINED
  solve_type = 'NEWTON'       # PJFNK not supported for rsls 
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 1.0
  num_steps = 8
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_LLG_test
    interval = 1
    elemental_as_nodal = true
  [../]
[]
