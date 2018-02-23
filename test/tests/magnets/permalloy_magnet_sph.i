## BENCHMARK PROBLEM #1
## DOES NOT WORK AT THE MOMENT
##
##


[Mesh]
  # Mesh details
  type = GeneratedMesh
  dim = 3
  nx = 12
  ny = 12
  nz = 2
  elem_type = HEX8


  # dimensions in units of microns
  xmin = -0.0075
  xmax = 0.0075
  ymin = -0.0075
  ymax = 0.0075
  zmin = -0.0005
  zmax = 0.0005
[]

[GlobalParams]
  polar_theta = polar_theta
  azimuth_phi = azimuth_phi
 
  potential_H_int = potential_H_int
  potential_H_ext = potential_H_ext

  alpha = 0.85

  Ae = 1.0e-3
  M = 11.43
  K1 = -1.0e-1
  K2 = 2.25e-1
  g0 = 1.0
  permittivity = 1.0
  mu0 = 1.0
[]

[Variables]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 1.5
      max = 1.55
      seed = 5
    [../]
  [../]
  [./azimuth_phi]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 1.5 #6.283185307179586
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
    type = DampedConstrainedExchangeLLG
    variable = polar_theta
    component = 0
  [../]
  [./d_llg_exch_phi]
    type = DampedConstrainedExchangeLLG
    variable = azimuth_phi
    component = 1
  [../]
  [./d_llg_aniso_th]
    type = DampedConstrainedAnisotropyLLG
    variable = polar_theta
    component = 0
  [../]
  [./d_llg_aniso_phi]
    type = DampedConstrainedAnisotropyLLG
    variable = azimuth_phi
    component = 1
  [../]

  # Magnetic interaction term

  [./d_HM_1]
    type = MagMStrong
    variable = azimuth_phi
    component = 1
  [../]
  [./d_HM_0]
    type = MagMStrong
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
  [./normalizedtotalsaturation]
    type = PolarizationValue
    polar_x = magnetic_x
    polar_y = magnetic_y
    polar_z = magnetic_z
    execute_on = 'initial timestep_end'
  [../]
  [./Fexchange]
    type = MagneticConstrainedExchangeEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./Faniso]
    type = MagneticConstrainedAnisotropyEnergy
    execute_on = 'initial timestep_end'
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[BCs]
 # [./bc_az]
 #   type = NeumannBC
 #   variable = azimuth_phi
 #   value = 0.0
 #   boundary = 'front back left right top bottom'
 # [../]
 # [./bc_pt]
 #   type = NeumannBC
 #   variable = polar_theta
 #   value = 0.0
 #   boundary = 'front back left right top bottom'
 # [../]
 # [./bc_int_pot]
 #   type = NeumannBC
 #   variable = potential_H_int
 #   value = 0.0
 #   boundary = 'front back left right top bottom'
 # [../]

  [./bc_ext_pot_front]
    type = DirichletBC
    variable = potential_H_ext
    value = 1.5
    boundary = 'front'
  [../]
  [./bc_ext_pot_back]
    type = DirichletBC
    variable = potential_H_ext
    value = 0.0
    boundary = 'back'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type '
    petsc_options_value = '     121              1e-10        1e-8      1e-6     lu     '
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.00001               #NOTE FOR MAGNETS, TIMESTEP AND ALPHA CHOICE ARE INTERTWINED
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'bdf2'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 1.0
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_mag_sph_test_high_H2_bdf_damp
    interval = 1
    elemental_as_nodal = true
  [../]
[]
