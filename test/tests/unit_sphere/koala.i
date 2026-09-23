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
  [../]
  [./laplace]
    type = FullSolveMultiApp
    input_files = laplace.i
  [../]
[]

[Transfers]
  [./from_sub]
    type = MultiAppCopyTransfer
    source_variable = phi
    variable = phi
    from_multi_app = poisson
  [../]
  [./to_sub]
    type = MultiAppCopyTransfer
    source_variable = phi
    variable = phi1
    to_multi_app = laplace
  [../]
  [./from_sub2]
    type = MultiAppAddTransfer
    source_variable = phi
    variable = phi
    from_multi_app = laplace
  [../]
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  execute_on = timestep_end
[]
