[Mesh]
  file = sphere_tet_approx_size0_05.e
  dim = 3
[]

[Variables]
  [./potential_H_int]
    order = FIRST
    family = LAGRANGE
    initial_condition = 1.0
  [../]
[]

[AuxVariables]
  [./potential_H_int1]   #the transfer system just puts this here :)
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Materials]
  [./permitivitty_1]
    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '1.0'
  [../]
[]

[Kernels]
  [./diff]
    type = Electrostatics
    variable = potential_H_int
  [../]
[]

[BCs]
  [./dirichlet]
    type = CoupledDirichletBC
    variable = potential_H_int
    boundary = 1
    coupled_var = potential_H_int1
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
  print_linear_residuals = false
  execute_on = timestep_end
[]
