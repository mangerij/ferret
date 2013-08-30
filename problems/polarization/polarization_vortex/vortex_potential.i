[Mesh]
  file = puck.e
  uniform_refine = 0
[]

[Variables]
  active = 'Phi'
  [./Phi]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  active = 'Poisson'
  [./Poisson]
    type = Diffusion
    variable = Phi
  [../]

[]

[BCs]
  active = 'vortex_surface_charge'
  [./vortex_surface_charge]
    type = VortexSurfaceCharge
    variable = Phi
    boundary = '1'
    a        = 1.0
  [../]
[]

[Executioner]
  type = Steady

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options = '-snes_monitor -ksp_monitor -pc_svd_monitor'
  petsc_options_iname = '-pc_type'
  petsc_options_value = '     svd'
[]

[Output]
  file_base = vortex_potential
  interval = 1
  exodus = true
  perf_log = true
[]


