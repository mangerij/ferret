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
  nz = 6
  elem_type = HEX8


  # dimensions in units of microns
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  zmin = -0.1
  zmax = 0.1
[]

[GlobalParams]
  len_scale = 1.0
  mag_x = magnetic_x
  mag_y = magnetic_y
  mag_z = magnetic_z
  lambda = lambda

  #M0 = 80.0 #microns!

  alphaLL = 1.0

  nx = 1.0
  ny = 0.0
  nz = 0.0

  eps = 1.0e-7

  Ku = 0.00078

  A = 20.31

  phi = phi
  theta = theta
  potential_H_int = potential_H_int

  var_mag = 1e-7
  sat_penalty = 1e2
[]

[Variables]
  [./potential_H_int]
    order = FIRST
    family = LAGRANGE
  [../]

  [./magnetic_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      component = 0
    [../]
      scaling = 1e-5
  [../]
  [./magnetic_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      component = 1
    [../]
      scaling = 1e-5
  [../]
  [./magnetic_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      component = 2
    [../]
      scaling = 1e-5
  [../]
 # [./lambda]
 #   order = FIRST
 #   family = LAGRANGE
 # [../]
[]

[AuxVariables]
  [./phi]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 360.0
    [../]
  [../]
  [./theta]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 360.0
    [../]
  [../]
  #[./magneticsq_x]
  #  order = FIRST
  #  family = LAGRANGE
  #[../]
  #[./magneticsq_y]
  #  order = FIRST
  #  family = LAGRANGE
  #[../]
  #[./magneticsq_z]
  #  order = FIRST
  #  family = LAGRANGE
  #[../]
[]

#[AuxKernels]
# [./mag_x_sq_c]
#   type = PzSq
#   variable = magneticsq_x
#   polar_z = magnetic_x
# [../]
# [./mag_y_sq_c]
#   type = PzSq
#   variable = magneticsq_y
#   polar_z = magnetic_y
# [../]
# [./mag_z_sq_c]
#   type = PzSq
#   variable = magneticsq_z
#   polar_z = magnetic_z
# [../]
#[]

[Kernels]
  ## Time dependence
  [./mag_x_time]
    type = TimeDerivativeScaled
    variable = magnetic_x
    time_scale = 1.0
  [../]
  [./mag_y_time]
    type = TimeDerivativeScaled
    variable = magnetic_y
    time_scale = 1.0
  [../]
  [./mag_z_time]
    type = TimeDerivativeScaled
    variable = magnetic_z
    time_scale = 1.0
  [../]

  ## Magnetostatic terms:
  [./mag_h]
     type = MagHStrong
     variable = potential_H_int
  [../]
  [./M_int]
     type = Electrostatics
     variable = potential_H_int
     permittivity = 1.0
  [../]
  
  [./mag_x]
     type = MagMStrong
     variable = magnetic_x
     component = 0
  [../]
  [./mag_y]
     type = MagMStrong
     variable = magnetic_y
    component = 1
  [../]
  [./mag_z]
     type = MagMStrong
     variable = magnetic_z
     component = 2
  [../]

  ## LLG precessional terms
  #
  # to be implemented.
  #
  #
  ##

  # LLG damping terms
  [./d_mag_aniso_x]
    type = DampedAnisotropyLLG
    variable = magnetic_x
    component = 0
  [../]
  [./d_mag_aniso_y]
    type = DampedAnisotropyLLG
    variable = magnetic_y
    component = 1
  [../]
  [./d_mag_aniso_z]
    type = DampedAnisotropyLLG
    variable = magnetic_z
    component = 2
  [../]

  [./d_mag_exch_x]
    type = DampedExchangeLLG
    variable = magnetic_x
    component = 0
  [../]
  [./d_mag_exch_y]
    type = DampedExchangeLLG
    variable = magnetic_y
    component = 1
  [../]
  [./d_mag_exch_z]
    type = DampedExchangeLLG
    variable = magnetic_z
    component = 2
  [../]

  ## Magnetic constraint kernels
 # [./Lx]
 #   type = LagrangeMagConstraint
 #   component  = 0
 #   variable = magnetic_x
 # [../]
 # [./Ly]
 #   type = LagrangeMagConstraint
 #   component  = 1
 #   variable = magnetic_y
 # [../]
 # [./Lz]
 #   type = LagrangeMagConstraint
 #   component  = 2
 #   variable = magnetic_z
 # [../]
 # [./LL]
 #   type = LagrangeLambdaConstraint
 #   variable = lambda
 # [../]
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
  [./Fexch]
    type = MagneticExchangeEnergy
    execute_on = timestep_end
  [../]
  [./Faniso]
    type = MagneticAnisotropyEnergy
    execute_on = timestep_end
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[BCs]
  #[./lambda_back]
  #  type = DirichletBC
  #  variable = lambda
  #  boundary = 'back'
  #  value = 0.0
  #[../]
  #[./lambda_front]
  #  type = DirichletBC
  #  variable = lambda
  #  boundary = 'front'
  #  value = 0.0
  #[../]
  #[./lambda_top]
  #  type = DirichletBC
  #  variable = lambda
  #  boundary = 'top'
  #  value = 0.0
  #[../]
  #[./lambda_bottom]
  #  type = DirichletBC
  #  variable = lambda
  #  boundary = 'bottom'
  #  value = 0.0
  #[../]
  #[./lambda_left]
  #  type = DirichletBC
  #  variable = lambda
  #  boundary = 'left'
  #  value = 0.0
  #[../]
  #[./lambda_right]
  #  type = DirichletBC
  #  variable = lambda
  #  boundary = 'right'
  #  value = 0.0
  #[../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol  -pc_type -sub_pc_type -pc_asm_overlap'
    petsc_options_value = '     121              1e-10         1e-8      1e-8     asm       lu              2'
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.1
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 1.0
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_mag_test
    elemental_as_nodal = true
  [../]
[]
