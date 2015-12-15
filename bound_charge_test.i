[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 3
  ny = 3
  nz = 3
  ymax = 2
  ymin = -2
  xmin = -2
  xmax = 2
  zmin = -2
  zmax = 2
[]

[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[GlobalParams]
  len_scale = 1.0
  alpha1 = -0.1722883 # (3.766(T-765.1)*10^5) C^{-2} nm^2 (T = 673 K)
  alpha11 = -0.07253
  alpha111 = 0.26
  alpha12 = 0.75
  alpha112 = 0.61
  alpha123 = -3.67
  permittivity = 0.5843763
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
[]

[Kernels]
  [./bed_x]
    type = BulkEnergyDerivativeSixth
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = BulkEnergyDerivativeSixth
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = BulkEnergyDerivativeSixth
    variable = polar_z
    component = 2
  [../]
  [./polar_x_time]
     type=TimeDerivative
     variable = polar_x
  [../]
  [./polar_y_time]
     type=TimeDerivative
     variable = polar_y
  [../]
  [./polar_z_time]
     type=TimeDerivative
     variable = polar_z
  [../]
[]


[AuxVariables]
  [./rho_b]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[AuxKernels]
  type = BoundCharge
  variable = rho_b
[]

[ICs]
  [./polar_x_constic_rand]
     type = RandomIC
     variable = polar_x
     min = 0.0003
     max = 0.0004
  [../]
  [./polar_y_constic_rand]
     type = RandomIC
     variable = polar_y
     min = 0.0003
     max = 0.0004
  [../]
  [./polar_z_constic_rand]
     type = RandomIC
     variable = polar_z
     min = 0.0003
     max = 0.0004
  [../]
[]

[Executioner]
  type=Transient

  #[./TimeStepper]
  #  type = IterationAdaptiveDT
  #  dt = 1.1e-10
  #  optimal_iterations = 3
  #  growth_factor = 1.001
  #  cutback_factor =  0.999
  #[../]
  scheme = 'explicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1.0e-29
  dtmax = 1.70e-9
  dt = 1.1e-3
  num_steps = 2
  petsc_options='-snes_converged_reason'
  petsc_options_iname='-ksp_type -snes_type  -snes_rtol -ksp_rtol -pc_type  -pc_factor_zeropivot'
  petsc_options_value=' gmres     newtonls       1e-8     1e-10      gamg           1e-50  '
[]


[Outputs]
  print_linear_residuals = true
  print_nonlinear_residuals = true
  print_perf_log = true

  [./out]
    type = Exodus
    file_base = out_bound_charge_test
    output_initial = true
    elemental_as_nodal = false
  [../]
[]
