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

[Executioner]
  type = Transient
  num_steps = 1
[]

[MultiApps]
  [./poisson]
    type = FullSolveMultiApp
    input_files = poisson.i
    execute_on = initial
  [../]
  [./laplace]
    type = FullSolveMultiApp
    input_files = laplace.i
    execute_on = timestep_begin
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

[Outputs]
  exodus = true
  execute_on = timestep_end
[]
