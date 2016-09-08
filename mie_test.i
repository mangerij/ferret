[Mesh]
  file = mie_test.e
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./E_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./E_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./E_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./mie_field_e_x]
    type = MieFieldReals
    variable = E_x
    execute_on = 'timestep_end'
    a = 450
    omega = 1.0
    c = 1.0
    epsilonI = 1.0 
    sigmaI = 0.0
    epsilonII = 1.0
    sigmaII = 0.0
    L = 450
    order = 17
    component = 0
    block = '2'
  [../]
  [./mie_field_e_y]
    type = MieFieldReals
    variable = E_y
    execute_on = 'timestep_end'
    a = 450
    omega = 1.0
    c = 1.0
    epsilonI = 1.0 
    sigmaI = 0.0
    epsilonII = 1.0
    sigmaII = 0.0
    L = 450
    order = 17
    component = 1
    block = '2'
  [../]
  [./mie_field_e_z]
    type = MieFieldReals
    variable = E_z
    execute_on = 'timestep_end'
    a = 450
    omega = 1.0
    c = 1.0
    epsilonI = 1.0 
    sigmaI = 0.0
    epsilonII = 1.0
    sigmaII = 0.0
    L = 450
    order = 17
    component = 2
    block = '2'
  [../]
[]

[Kernels]
  [./u]
    type = Diffusion
    variable = u
  [../]
[]

[BCs]
  [./top]
    type = DirichletBC
    variable = u
    boundary = '1'
    value = 0.0
  [../]
  [./bot]
    type = DirichletBC
    variable = u
    boundary = '2'
    value = 1
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-snes_rtol -ksp_rtol'
    petsc_options_value = '  1e-1      1e-1 '
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
[]

[Outputs]
  print_perf_log = true
  [./out]
    type = Exodus
    elemental_as_nodal = true
    execute_on = 'timestep_end'
  [../]
[]

