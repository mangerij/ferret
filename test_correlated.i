
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 50
  xmin = 0.0
  xmax = 50.0
  ymin = 0.0
  ymax = 50.0
[]

[Variables]
  [./test_var]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = CorrelatedRandomFieldIC
      Lcorr = 10.0
      dim = 2
      Nnodes = 100
      seed = 2                       #internal seeding works properly
    [../]
  [../]
[]

[Kernels]
  ## Time dependence
  [./polar_time]
    type = TimeDerivativeScaled
    variable = test_var
    time_scale = 1.0
  [../]
  [./diff]
    type = Diffusion
    variable = test_var
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  dt = 1.0
  solve_type = 'NEWTON'
  dtmin = 1e-5
  dtmax = 10.0
  num_steps = 1
[]

[Outputs]
  print_linear_residuals = true
  [./out]
    type = Exodus
    file_base = test_N100_Lc100
    interval = 1
    elemental_as_nodal = true
  [../]
[]
