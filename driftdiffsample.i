[Mesh]
    type = GeneratedMesh
    dim = 2
    nx = 300
    ny = 60
    xmin = 0
    xmax = 10
    ymin = 0
    ymax = 2
[]

[Variables]
  [./n]
    order = FIRST
    family = LAGRANGE
    # [./InitialCondition]
    #   type = SmoothCircleIC
    #   x1 = 0
    #   y1 = 10
    #   radius = 10
    #   invalue = 0.5
    #   outvalue = 0
    # [../]
  [../]
  [./p]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
    #   type = SmoothCircleIC
    #   x1 = 50
    #   y1 = 10
    #   radius = 3
    #   invalue = 0.25
    #   outvalue = 0
    # [../]
  [../]
  [./phi]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./Diffusion_e]
    type = ElectronDiffusion
    variable = n
  [../]
  [./Drift_e]
    type = ElectronDrift
    variable = phi
    n = n
  [../]
  # [./Time_e]
  #   type = TimeDerivative
  #   variable = n
  # [../]

  [./Diffusion_p]
    type = HoleDiffusion
    variable = p
  [../]
  [./Drift_p]
    type = HoleDrift
    variable = phi
    p = p
  [../]
  # [./Time_p]
  #   type = TimeDerivative
  #   variable = p
  # [../]

  [./ElectrostaticsCombo]
    type = ElectrostaticsCombo
    variable = phi
    n = n
    p = p
    permittivity = 1
  [../]
[]

[BCs]
  [./bottom_convected]
    type = DirichletBC
    variable = phi
    boundary = 'left'
    value = 1
  [../]
  [./top_convected]
    type = DirichletBC
    variable = phi
    boundary = 'right'
    value = 0
  [../]
  [./bottom_diffused]
    type = DirichletBC
    variable = n
    boundary = 'bottom'
    value = 1
  [../]
  [./top_diffused]
    type = DirichletBC
    variable = n
    boundary = 'top'
    value = 0
  [../]
  [./bottom1_diffused]
    type = DirichletBC
    variable = p
    boundary = 'bottom'
    value = 0
  [../]
  [./top2_diffused]
    type = DirichletBC
    variable = p
    boundary = 'top'
    value = 1
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
   # dt = 0.000001
[]


[Outputs]
  file_base = out
  exodus = true
  [./console]
    type = Console
    perf_log = true
    linear_residuals = false
  [../]
[]
