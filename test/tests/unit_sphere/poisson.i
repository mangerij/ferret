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

[AuxVariables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]
[]

[Kernels]
  [./diffusion]
     type = Electrostatics
     variable = phi
     permittivity = 1.0
  [../]
  [./forcing]
    type = PolarElectricEStrong
    variable = phi
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
  [../]
[]

[Executioner]
  type = Steady
#  solve_type = 'JFNK'
  nl_abs_tol = 1e-8
[]

[Outputs]
  print_linear_residuals = true
  exodus = true
  execute_on = timestep_end
[]
