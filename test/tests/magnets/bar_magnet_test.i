## BENCHMARK PROBLEM #1
## DOES NOT WORK AT THE MOMENT
##
##

[Problem]
  type = FerretProblem
  polar_var = polar_theta
  azimuth_phi_var = azimuth_phi
  execute_on = 'Initial Linear Nonlinear'
[]

#[Mesh]
#  file = exodus_brick.e
#[]

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
      min = 0.0
      max = 1.0

    [../]
  [../]
  [./azimuth_phi]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 1.0
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

  [./demag_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./demag_y]
    order = CONSTANT
    family = MONOMIAL

  [../]
  [./demag_z]
    order = CONSTANT
    family = MONOMIAL
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
 [../]
 [./dmag_z_c]
   type = ElecFieldAux
   component = 2
   variable = demag_z
   execute_on = 'initial timestep_end'
   potential_E_int = potential_H_int
   potential_E_ext = potential_H_ext
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
    boundary = 'left'
  [../]

  [./bc_int_pot_back]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = 'right'
  [../]

  [./bc_ext_pot_back]
    type = DirichletBC
    variable = potential_H_ext
    value = 0.0
    boundary = 'left'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type '
    petsc_options_value = '     121             1e-10       1e-10      1e-8      lu'
  [../]
[]

[Executioner]
  type = Transient
  dt = 1e-4               # NOTE FOR MAGNETS, TIMESTEP AND ALPHA CHOICE ARE INTERTWINED
  solve_type = 'NEWTON'   # "PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'         # "implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-10
  dtmax = 1.0e-2
[]

[Debug]
  show_var_residual_norms = false
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_mag
    interval = 1
    elemental_as_nodal = true
  [../]
[]
