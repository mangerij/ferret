[GlobalParams]
  Debug = 'PolarizationSurfaceCharge::computeQpResidual'
[]

[Mesh]
  file = ../../../mesh/puck.R5.L1.i12.e
  uniform_refine = 0
[]

[Variables]
  [./phi]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./Diffusion]
    type = Diffusion
    variable = phi
  [../]
[]

[BCs]
  [./DoubleLayerCharge]
    type = PolarizationSurfaceCharge
    variable = phi
    P = 'P_x P_y P_z'
    boundary = '3'
  [../]
[]

[AuxVariables]
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
[PolarizationVortexAux]
  P_x = P_x
  P_y = P_y
  P_z = P_z
  a_y = 0.5
  c = 2.0
  R = 5.0
  L = 1.0
[]

[Executioner]
  type = Transient
  dt = 0.02
  start_time = -0.02
  end_time =  1.0

  print_linear_residuals = true
  petsc_options = '-snes_monitor -pc_svd_monitor -options_monitor'
  petsc_options_iname = '-ksp_type   -pc_type -pc_factor_levels'
  petsc_options_value = '    gmres        ilu                 5'
[]

[Output]
  file_base = polarization_vortex.0.R5.L1.c2.out
  interval = 1
  exodus = true
  perf_log = true
[]


