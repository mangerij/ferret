[Variables]
  [mag_x]
    block = film
    [InitialCondition]
      type = FunctionIC
      function = ic_mag_x
    []
  []
  [mag_z]
    block = film
    [InitialCondition]
      type = FunctionIC
      function = ic_mag_z
    []
  []
  [potential_H_int]
    block = 'film vacuum'
    [InitialCondition]
      type = FunctionIC
      function = ic_potential
    []
  []
[]

[Kernels]
  [mag_x_time]
    type = TimeDerivative
    variable = mag_x
    block = film
  []
  [mag_z_time]
    type = TimeDerivative
    variable = mag_z
    block = film
  []
  [magnetostatic_source]
    type = MagHStrongCart
    variable = potential_H_int
    mag_y = 0.0
    block = film
  []
[]

[Postprocessors]
  [gx]
    type = PointValue
    variable = mag_x
    point = '0.0 0.0 0.0'
    execute_on = 'initial timestep_end'
  []
  [gz]
    type = PointValue
    variable = mag_z
    point = '0.0 0.0 0.0'
    execute_on = 'initial timestep_end'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  [TimeIntegrator]
    type = CrankNicolson
  []
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-30
  dt = ${dt}
  num_steps = ${num_steps}
[]
