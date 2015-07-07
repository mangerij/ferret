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
  [./potential_int]
    order=FIRST
    family = LAGRANGE
  [../]
[]

[GlobalParams]
   permittivity=8.85e-12
   len_scale = 1e10
   potential_int=potential_int
[]

[Kernels]
  [./FE_E_int]
     type=Electrostatics
     variable=potential_int
  [../]
[]

[BCs]
  [./potential_int_upz]
    type = DirichletBC
    variable = potential_int
    boundary = 'top'
    value = 1.0
  [../]
  [./potential_int_downz]
    type = DirichletBC
    variable = potential_int
    boundary = 'bottom'
    value = 0.0
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'newton'
  petsc_options = '-snes_converged_reason -snes_check_jacobian'
  petsc_options_iname = '-ksp_type -snes_type  -snes_rtol -ksp_rtol -pc_type -pc_hypre_boomeramg_strong_threshold -snes_linesearch_type -pc_factor_zeropivot '
  petsc_options_value = ' gmres     newtonls       1e-8     1e-12      hypre     0.5                                 basic         1e-50  '
[]


[Outputs]
  print_linear_residuals = true
  print_perf_log = true

  [./out]
    type = Exodus
    file_base = out_electrostatic_test
    output_initial = true
    elemental_as_nodal = false
  [../]

[]
