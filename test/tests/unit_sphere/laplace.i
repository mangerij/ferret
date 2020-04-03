[Mesh]
  file = sphere_tet_approx_size0_05.e
  dim = 3
[]

[Variables]
  [./phi]
    order = FIRST
    family = LAGRANGE
    initial_condition = 1.0
  [../]
[]

[AuxVariables]
  [./phi1]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./diff]
    type = Electrostatics
    variable = phi
    permittivity = 1.0
  [../]
[]

[BCs]
  [./dirichlet]
    type = CoupledDirichletBC
    variable = phi
    boundary = 1

    coupled_var = phi1
  [../]
[]


[UserObjects]
  [./bifmm]
    type = BoundaryIntegralFMM
    cx = 0.0
    cy = 0.0
    cz = 0.0
    boxWidth = 2.1
    TreeHeight = 5
    execute_on = timestep_begin
  [../]
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  execute_on = timestep_end
[]
