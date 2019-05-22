
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100 #number of elements
  xmin = 0.0
  xmax = 100.0
[]

[GlobalParams]
  xmin = 0.0
  xmax = 100.0
  ymin = 0.0
  ymax = 100.0
  zmin = 0.0
  zmax = 100.0
[]

[Variables]
  [./test_var]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = CorrelatedRandomFieldIC
      Lcorr = 10.0
      dim = 3
      Nnodes = 25
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

[BCs]
 [./dc]
   type = DirichletBC
   boundary = right
   variable = test_var
   value = 1.0
 [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason'
    #petsc_options_iname = '-snes_atol'
    #petsc_options_value = '1e-10'
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
    file_base = test_N100_NNod100_Lc5_d3_evaluateLQ
    interval = 1
    elemental_as_nodal = true
  [../]
[]
