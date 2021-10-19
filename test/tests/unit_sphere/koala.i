[Mesh]
  file = sphere_tet_approx_size0_05.e
  dim = 3
[]

[Variables]
  [./phi]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Problem]
  type = FEProblem
  solve = false
[]

[MultiApps]
  [./poisson]
    type = FullSolveMultiApp
    input_files = poisson.i
    #execute_on = 'initial linear nonlinear timestep_begin'
  [../]
  [./laplace]
    type = FullSolveMultiApp
    input_files = laplace.i
    #execute_on = 'initial linear nonlinear timestep_begin'
  [../]
[]

[Transfers]
  [./from_sub]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = phi
    variable = phi
    multi_app = poisson
  [../]
  [./to_sub]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    source_variable = phi
    variable = phi1
    multi_app = laplace
  [../]
  [./from_sub2]
    type = MultiAppAddTransfer
    direction = from_multiapp
    source_variable = phi
    variable = phi
    multi_app = laplace

  [../]
[]

[Executioner]
  type = Steady
  num_steps = 1
[]


[Outputs]
  exodus = true
  execute_on = timestep_end
[]
