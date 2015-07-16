[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 2
  ny = 2
  nz = 2
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
   len_scale = 1e-9
   alpha1 = -1.8202e8 # 3.8(T-479)*10^5 C^{-2}m^2 (T=0 K)
   alpha11 = -7.3e7
   alpha111 = 2.6e8
   alpha12 = 7.5e8
   alpha112 = 6.1e8
   alpha123 = -3.7e9
   permittivity = 8.85e-12
   polar_x = polar_x
   polar_y = polar_y
   polar_z = polar_z
   potential_int = potential_int
   potential_ext = potential_ext
   time_scale = 1.0e-29
[]

[Kernels]
  [./bed_x]
    type = BulkEnergyDerivative
    variable = polar_x
    component=0
  [../]
  [./bed_y]
    type = BulkEnergyDerivative
    variable = polar_y
    component=1
  [../]
  [./bed_z]
    type = BulkEnergyDerivative
    variable = polar_z
    component = 2
  [../]
  [./polar_x_time]
     type=TimeDerivativeScaled
     variable = polar_x
  [../]
  [./polar_y_time]
     type=TimeDerivativeScaled
     variable = polar_y
  [../]
  [./polar_z_time]
     type=TimeDerivativeScaled
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
     min = 0.6
     max = 0.75
  [../]
  [./polar_y_constic_rand]
     type = RandomIC
     variable = polar_y
     min = 0.6
     max = 0.75
  [../]
  [./polar_z_constic_rand]
     type = RandomIC
     variable = polar_z
     min = 0.6
     max = 0.65
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
  dtmin=1.0e-29
  dtmax=1.70e-9
  dt = 1.1e-12
  num_steps = 250
  petsc_options='-snes_converged_reason'
  petsc_options_iname='-ksp_type -snes_type  -snes_rtol -ksp_rtol -pc_type -pc_hypre_boomeramg_strong_threshold -snes_linesearch_type -pc_factor_zeropivot'
  petsc_options_value=' gmres     newtonls       1e-8     1e-12      hypre     0.5                                 basic         1e-50  '
[]


[Outputs]
  print_linear_residuals = true
  print_perf_log = true

  [./out]
    type = Exodus
    file_base = out_bound_charge_test
    output_initial = true
    elemental_as_nodal = false
    interval = 50
  [../]

[]
