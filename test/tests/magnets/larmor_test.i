
[Mesh]
  file = exodus_sphere.e
[]

[GlobalParams]
  polar_theta = polar_theta
  azimuth_phi = azimuth_phi
 
  potential_H_int = potential_H_int
  potential_H_ext = potential_H_ext

  alpha = 0.0
  Ms = 0.8

  g0 = 175929.0

  permittivity = 0.000795775 #This is equal to 1/mu0 in these units.
  mu0 = 1256.64
[]

[Variables]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 0.60000000
      max = 0.60000001
      seed = 3 
    [../]
  [../]
  [./azimuth_phi]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 6.283185307177
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

  # Magnetostatic equations (this just solves for a constant field along some direction).

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
[]

[BCs]
   [./bc_int_pot_R]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '3'
  [../]

  [./bc_int_pot_L]
    type = DirichletBC
    variable = potential_H_int
    value = 1000.0
    boundary = '5'
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

[Postprocessors]
  [./Fmag]
    type = MagnetostaticEnergy
    execute_on = 'initial timestep_end'
    outputs = 'none'
    block = '1'
  [../]
[]


[Bounds]
  [./u_bounds]
    type = BoundsAux
    variable = bounds_theta_dummy
    bounded_variable = polar_theta
    upper = 3.14159265357
    lower = 0.00000000001
    execute_on = 'Initial Linear Nonlinear'
  [../]

  [./v_bounds]
    type = BoundsAux
    variable = bounds_phi_dummy
    bounded_variable = azimuth_phi
    upper = 6.283185307177
    lower = 0.000000000001
    execute_on = 'Initial Linear Nonlinear'
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
  solve_type = 'NEWTON'
  scheme = 'implicit-euler'

  dtmin = 1.0e-16
  dtmax = 1.0e-6
  dt = 1.0e-6
  num_steps = 5
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = Larmor_USLLG_test_A0
    interval = 1
    elemental_as_nodal = true
  [../]
  #[./pgraph]
  #  type = PerfGraphOutput
  #  level = 1
  #  heaviest_branch = true
  #  heaviest_sections = 10
  #[]
[]
