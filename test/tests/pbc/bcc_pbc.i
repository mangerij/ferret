[Mesh]
  file = supercell.e
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
    block = '1 2 3'
  [../]
[]

[Functions]

  [./ic_func]
    type = ParsedFunction
    expression = 'sin(0.6283185307179586 * x)'
  [../]
  [./exact_func]
    type = ParsedFunction
    expression = 'sin(0.6283185307179586 * x) * exp(-0.3947841760435743 * t)'
  [../]
[]

[ICs]
  [./u_ic]
    type = FunctionIC
    variable = u
    function = ic_func
  [../]
[]

[Kernels]
  [./dudt]
    type = TimeDerivative
    variable = u
  [../]
  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[BCs]

  [./Periodic]
    [./xy]
      auto_direction = 'x y'
      variable = 'u'
    [../]
  [../]
[]

[Postprocessors]
  [./l2_error]
    type = ElementL2Error
    variable = u
    function = exact_func
    execute_on = 'initial timestep_end'
  [../]
  [./total_u]
    type = ElementIntegralVariablePostprocessor
    variable = u
    execute_on = 'initial timestep_end'
  [../]
  [./u_max]
    type = NodalExtremeValue
    variable = u
    value_type = max
    execute_on = 'initial timestep_end'
  [../]
  [./u_min]
    type = NodalExtremeValue
    variable = u
    value_type = min
    execute_on = 'initial timestep_end'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  scheme = 'crank-nicolson'

  petsc_options_iname = '-pc_type -ksp_rtol'
  petsc_options_value = 'lu       1e-12'

  dt = 0.25
  num_steps = 8
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  print_linear_residuals = false
  perf_graph = false
  [./out]
    type = Exodus
    file_base = out_bcc_pbc
  [../]
[]
