## BENCHMARK PROBLEM #1
## DOES NOT WORK AT THE MOMENT
##
##


[Mesh]
  # Mesh details
  type = GeneratedMesh
  dim = 3
  nx = 20
  ny = 20
  nz = 10
  elem_type = HEX8


  # dimensions in units of microns
  xmin = -10.0
  xmax = 10.0
  ymin = -10.0
  ymax = 10.0
  zmin = -5.0
  zmax = 5.0
[]

[GlobalParams]
  polar_theta = polar_theta
  azimuth_phi = azimuth_phi

  alpha = 1.0

  Ae = 1.0
  M = 11.348

  K1 = 1.0
  K2 = -2.25

  g0 = -1.0
[]

[Variables]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 360.0
      seed = 5
    [../]
  [../]
  [./azimuth_phi]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 360.0
      seed = 5
    [../]
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
 [../]
 [./mag_y_c]
   type = MagFieldAux
   component = 0
   variable = magnetic_y
 [../]
 [./mag_z_c]
   type = MagFieldAux
   component = 0
   variable = magnetic_z
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
[]

[Postprocessors]
  [./Mx_extreme]
    type = ElementExtremeValue
    variable = magnetic_x
    execute_on = 'initial timestep_end'
  [../]
  [./My_extreme]
    type = ElementExtremeValue
    variable = magnetic_y
    execute_on = 'initial timestep_end'
  [../]
  [./Mz_extreme]
    type = ElementExtremeValue
    variable = magnetic_z
    execute_on = 'initial timestep_end'
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[BCs]

[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol  -pc_type '
    petsc_options_value = '     121              1e-10         1e-8      1e-8     lu     '
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.00001
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'bdf2'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 1.0
  num_steps = 2
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_mag_sph_test
    interval = 5
    elemental_as_nodal = true
  [../]
[]
