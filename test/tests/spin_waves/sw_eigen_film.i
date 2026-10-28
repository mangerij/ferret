!include film_common.i

[Variables]
  [mag_x]
    block = film
  []
  [mag_z]
    block = film
  []
  [potential_H_int]
    block = 'film vacuum'
  []
[]

[Kernels]
  [mag_x_time]
    type = CoefReaction
    variable = mag_x
    coefficient = -1.0
    extra_vector_tags = 'eigen'
    block = film
  []
  [mag_z_time]
    type = CoefReaction
    variable = mag_z
    coefficient = -1.0
    extra_vector_tags = 'eigen'
    block = film
  []
  [magnetostatic_source]
    type = MagHStrongCart
    variable = potential_H_int
    block = film
  []
[]

[Problem]
  active_eigen_index = 2
[]

[VectorPostprocessors]
  [eigenvalues]
    type = Eigenvalues
    execute_on = 'timestep_end'
  []
[]

[Postprocessors]
  [eig_0]
    type = VectorPostprocessorComponent
    vectorpostprocessor = eigenvalues
    vector_name = eigen_values_imag
    index = 0
    execute_on = 'timestep_end'
    outputs = none
  []
  [omega_n0]
    type = ParsedPostprocessor
    expression = 'abs(eig_0)'
    pp_names = 'eig_0'
    execute_on = 'timestep_end'
  []
  [eig_2]
    type = VectorPostprocessorComponent
    vectorpostprocessor = eigenvalues
    vector_name = eigen_values_imag
    index = 2
    execute_on = 'timestep_end'
    outputs = none
  []
  [omega_n1]
    type = ParsedPostprocessor
    expression = 'abs(eig_2)'
    pp_names = 'eig_2'
    execute_on = 'timestep_end'
  []
  [eig_4]
    type = VectorPostprocessorComponent
    vectorpostprocessor = eigenvalues
    vector_name = eigen_values_imag
    index = 4
    execute_on = 'timestep_end'
    outputs = none
  []
  [omega_n2]
    type = ParsedPostprocessor
    expression = 'abs(eig_4)'
    pp_names = 'eig_4'
    execute_on = 'timestep_end'
  []
  [eig_6]
    type = VectorPostprocessorComponent
    vectorpostprocessor = eigenvalues
    vector_name = eigen_values_imag
    index = 6
    execute_on = 'timestep_end'
    outputs = none
  []
  [omega_n3]
    type = ParsedPostprocessor
    expression = 'abs(eig_6)'
    pp_names = 'eig_6'
    execute_on = 'timestep_end'
  []
  [rel_kittel]
    type = ParsedPostprocessor
    expression = 'omega_n0 / (140.206695 * 0.238434) - 1'
    pp_names = 'omega_n0'
    execute_on = 'timestep_end'
  []
[]

[UserObjects]
  [kittel]
    type = Terminator
    expression = '!(abs(rel_kittel) < 1e-10)'
    fail_mode = HARD
    error_level = ERROR
    message = 'k = 0 eigenvalue is not the Kittel frequency g0 (H0 - Ms)'
    execute_on = 'timestep_end'
  []
[]

[Executioner]
  type = Eigenvalue
  solve_type = KRYLOVSCHUR
  eigen_problem_type = GEN_NON_HERMITIAN
  which_eigen_pairs = TARGET_MAGNITUDE
  n_eigen_pairs = 8
  n_basis_vectors = 24
  eigen_tol = 1e-12
  initial_eigenvalue = 1e12
  petsc_options_iname = '-st_type -eps_target -st_pc_type -st_pc_factor_mat_solver_type'
  petsc_options_value = 'sinvert 0       lu          mumps'
[]

[Outputs]
  file_base = out_sw_eigen_film
  exodus = true
  csv = true
  execute_on = 'timestep_end'
[]
