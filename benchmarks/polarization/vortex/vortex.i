[Mesh]
  file = ../../../mesh/spuck.e
  uniform_refine = 2
[]

[Variables]
  [./P_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./P_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./P_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]
[PolarizationVortex]
  P_x = P_x
  P_y = P_y
  P_z = P_z
  a_y = 0.5
  c = 0.5
  R = 5.0
  L = 1.0
  debug = true
  kernel_debug = true
[]

[Executioner]
  type = Transient
  dt = 0.1
  start_time = -0.1
  end_time =  1.0
  petsc_options = '-snes_monitor -ksp_monitor -pc_svd_monitor -options_monitor'
  petsc_options_iname = '-ksp_type       -pc_type -pc_factor_levels'
  petsc_options_value = '         gmres                    ilu                                     3'
[]

[Output]
  file_base = vortex
  interval = 1
  exodus = true
  perf_log = true
[]


