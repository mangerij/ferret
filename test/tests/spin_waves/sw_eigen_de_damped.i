!include film_inplane_mesh_coarse.i

alpha = 0.05

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z
  potential_H_int = potential_H_int
  Hext_x = Hext_x
  Hext_y = Hext_y
  Hext_z = Hext_z
  mu0 = 1.256637
  g0 = 140.206695
  Hscale = 1.0
[]

[Problem]
  use_hash_table_matrix_assembly = true
[]

[Variables]
  [mag_x]
    block = film
  []
  [mag_y]
    block = film
  []
  [potential_H_int]
    block = 'film vacuum'
  []
  [lm_mag_x]
    block = film_secondary
  []
  [lm_mag_y]
    block = film_secondary
  []
  [lm_potential]
    block = ends_secondary
  []
[]

[AuxVariables]
  [mag_z]
    block = film
    [InitialCondition]
      type = ConstantIC
      value = 1.0
    []
  []
  [Hext_x]
    block = film
  []
  [Hext_y]
    block = film
  []
  [Hext_z]
    block = film
    [InitialCondition]
      type = ConstantIC
      value = 0.1
    []
  []
[]

[Materials]
  [film_properties]
    type = GenericConstantMaterial
    prop_names = 'alpha Ae Ms permittivity'
    prop_values = '${alpha} 1.0 1.261566 1.0'
    block = film
  []
  [vacuum_properties]
    type = GenericConstantMaterial
    prop_names = 'Ms permittivity'
    prop_values = '0.0 1.0'
    block = vacuum
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
  [mag_y_time]
    type = CoefReaction
    variable = mag_y
    coefficient = -1.0
    extra_vector_tags = 'eigen'
    block = film
  []
  [exchange_x]
    type = MasterExchangeCartLLG
    variable = mag_x
    component = 0
    block = film
  []
  [exchange_y]
    type = MasterExchangeCartLLG
    variable = mag_y
    component = 1
    block = film
  []
  [field_x]
    type = MasterInteractionCartLLGHConst
    variable = mag_x
    component = 0
    block = film
  []
  [field_y]
    type = MasterInteractionCartLLGHConst
    variable = mag_y
    component = 1
    block = film
  []
  [magnetostatic_laplace]
    type = Electrostatics
    variable = potential_H_int
    block = 'film vacuum'
  []
  [magnetostatic_source]
    type = MagHStrongCart
    variable = potential_H_int
    block = film
  []
[]

[Constraints]
  [periodic_mag_x]
    type = EqualValueConstraint
    variable = lm_mag_x
    secondary_variable = mag_x
    primary_boundary = film_right
    secondary_boundary = film_left
    primary_subdomain = film_primary
    secondary_subdomain = film_secondary
    periodic = true
  []
  [periodic_mag_y]
    type = EqualValueConstraint
    variable = lm_mag_y
    secondary_variable = mag_y
    primary_boundary = film_right
    secondary_boundary = film_left
    primary_subdomain = film_primary
    secondary_subdomain = film_secondary
    periodic = true
  []
  [periodic_potential]
    type = EqualValueConstraint
    variable = lm_potential
    secondary_variable = potential_H_int
    primary_boundary = right
    secondary_boundary = left
    primary_subdomain = ends_primary
    secondary_subdomain = ends_secondary
    periodic = true
  []
[]

[BCs]
  [phi_pin]
    type = DirichletBC
    variable = potential_H_int
    boundary = phi_pin
    value = 0.0
  []
[]

[VectorPostprocessors]
  [eigenvalues]
    type = Eigenvalues
    execute_on = 'timestep_end'
  []
[]

[Postprocessors]
  [eig_imag_0]
    type = VectorPostprocessorComponent
    vectorpostprocessor = eigenvalues
    vector_name = eigen_values_imag
    index = 0
    execute_on = 'timestep_end'
    outputs = none
  []
  [eig_real_0]
    type = VectorPostprocessorComponent
    vectorpostprocessor = eigenvalues
    vector_name = eigen_values_real
    index = 0
    execute_on = 'timestep_end'
    outputs = none
  []
  [omega_0]
    type = ParsedPostprocessor
    expression = 'abs(eig_imag_0)'
    pp_names = 'eig_imag_0'
    execute_on = 'timestep_end'
  []
  [gamma_0]
    type = ParsedPostprocessor
    expression = 'eig_real_0'
    pp_names = 'eig_real_0'
    execute_on = 'timestep_end'
  []
  [eig_imag_2]
    type = VectorPostprocessorComponent
    vectorpostprocessor = eigenvalues
    vector_name = eigen_values_imag
    index = 2
    execute_on = 'timestep_end'
    outputs = none
  []
  [eig_real_2]
    type = VectorPostprocessorComponent
    vectorpostprocessor = eigenvalues
    vector_name = eigen_values_real
    index = 2
    execute_on = 'timestep_end'
    outputs = none
  []
  [omega_2]
    type = ParsedPostprocessor
    expression = 'abs(eig_imag_2)'
    pp_names = 'eig_imag_2'
    execute_on = 'timestep_end'
  []
  [gamma_2]
    type = ParsedPostprocessor
    expression = 'eig_real_2'
    pp_names = 'eig_real_2'
    execute_on = 'timestep_end'
  []
  [eig_imag_6]
    type = VectorPostprocessorComponent
    vectorpostprocessor = eigenvalues
    vector_name = eigen_values_imag
    index = 6
    execute_on = 'timestep_end'
    outputs = none
  []
  [eig_real_6]
    type = VectorPostprocessorComponent
    vectorpostprocessor = eigenvalues
    vector_name = eigen_values_real
    index = 6
    execute_on = 'timestep_end'
    outputs = none
  []
  [omega_6]
    type = ParsedPostprocessor
    expression = 'abs(eig_imag_6)'
    pp_names = 'eig_imag_6'
    execute_on = 'timestep_end'
  []
  [gamma_6]
    type = ParsedPostprocessor
    expression = 'eig_real_6'
    pp_names = 'eig_real_6'
    execute_on = 'timestep_end'
  []
  [eig_imag_10]
    type = VectorPostprocessorComponent
    vectorpostprocessor = eigenvalues
    vector_name = eigen_values_imag
    index = 10
    execute_on = 'timestep_end'
    outputs = none
  []
  [eig_real_10]
    type = VectorPostprocessorComponent
    vectorpostprocessor = eigenvalues
    vector_name = eigen_values_real
    index = 10
    execute_on = 'timestep_end'
    outputs = none
  []
  [omega_10]
    type = ParsedPostprocessor
    expression = 'abs(eig_imag_10)'
    pp_names = 'eig_imag_10'
    execute_on = 'timestep_end'
  []
  [gamma_10]
    type = ParsedPostprocessor
    expression = 'eig_real_10'
    pp_names = 'eig_real_10'
    execute_on = 'timestep_end'
  []
  [rel_omega_kittel]
    type = ParsedPostprocessor
    expression = 'omega_0 / (g0 / (1 + a^2) * sqrt(H1 * H2 - a^2 * (H1 - H2)^2 / 4)) - 1'
    pp_names = 'omega_0'
    constant_names = 'g0 a H1 H2'
    constant_expressions = '140.206695 ${alpha} 0.1 ${fparse 0.1 + 1.261566}'
    execute_on = 'timestep_end'
  []
  [rel_gamma_kittel]
    type = ParsedPostprocessor
    expression = 'gamma_0 / (g0 * a * (H1 + H2) / (2 * (1 + a^2))) - 1'
    pp_names = 'gamma_0'
    constant_names = 'g0 a H1 H2'
    constant_expressions = '140.206695 ${alpha} 0.1 ${fparse 0.1 + 1.261566}'
    execute_on = 'timestep_end'
  []
[]

[UserObjects]
  [kittel_omega]
    type = Terminator
    expression = '!(abs(rel_omega_kittel) < 1e-10)'
    fail_mode = HARD
    error_level = ERROR
    message = 'damped k = 0 frequency differs from the closed form'
    execute_on = 'timestep_end'
  []
  [kittel_gamma]
    type = Terminator
    expression = '!(abs(rel_gamma_kittel) < 1e-10)'
    fail_mode = HARD
    error_level = ERROR
    message = 'damped k = 0 linewidth differs from the closed form'
    execute_on = 'timestep_end'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Eigenvalue
  solve_type = KRYLOVSCHUR
  eigen_problem_type = GEN_NON_HERMITIAN
  which_eigen_pairs = TARGET_MAGNITUDE
  n_eigen_pairs = 12
  n_basis_vectors = 36
  eigen_tol = 1e-12
  initial_eigenvalue = 1e12
  petsc_options_iname = '-st_type -eps_target -st_pc_type -st_pc_factor_mat_solver_type'
  petsc_options_value = 'sinvert 0       lu          mumps'
[]

[Outputs]
  file_base = out_sw_eigen_de_damped
  csv = true
  execute_on = 'timestep_end'
[]
