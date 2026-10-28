!include film_common.i
!include transient_common.i

eigenvector = out_sw_eigen_film.e
eigenvalues = out_sw_eigen_film_eigenvalues_0001.csv
eigen_index = 2
amplitude = 1e-6
dt = 0.01
num_steps = 200
tolerance = 1e-9

[UserObjects]
  [eigenvector]
    type = SolutionUserObject
    mesh = ${eigenvector}
    system_variables = 'mag_x mag_z potential_H_int'
    timestep = LATEST
  []
  [agreement]
    type = Terminator
    expression = '!(abs(rel_ringdown) < ${tolerance})'
    fail_mode = HARD
    error_level = ERROR
    message = 'ring-down frequency (CN-corrected) and eigenvalue disagree by more than the tolerance'
    execute_on = 'final'
  []
[]

[Functions]
  [ic_mag_x]
    type = SolutionFunction
    solution = eigenvector
    from_variable = mag_x
    scale_factor = ${amplitude}
  []
  [ic_mag_z]
    type = SolutionFunction
    solution = eigenvector
    from_variable = mag_z
    scale_factor = ${amplitude}
  []
  [ic_potential]
    type = SolutionFunction
    solution = eigenvector
    from_variable = potential_H_int
    scale_factor = ${amplitude}
  []
[]

[VectorPostprocessors]
  [eigenvalues]
    type = CSVReaderVectorPostprocessor
    csv_file = ${eigenvalues}
    outputs = none
  []
[]

[Postprocessors]
  [omega_eigen]
    type = VectorPostprocessorComponent
    vectorpostprocessor = eigenvalues
    vector_name = eigen_values_imag
    index = ${eigen_index}
    execute_on = 'initial'
  []
  [d1x]
    type = ChangeOverTimePostprocessor
    postprocessor = gx
    execute_on = 'initial timestep_end'
  []
  [d1z]
    type = ChangeOverTimePostprocessor
    postprocessor = gz
    execute_on = 'initial timestep_end'
  []
  [d2x]
    type = ChangeOverTimePostprocessor
    postprocessor = d1x
    execute_on = 'initial timestep_end'
  []
  [d2z]
    type = ChangeOverTimePostprocessor
    postprocessor = d1z
    execute_on = 'initial timestep_end'
  []
  [num_step]
    type = ParsedPostprocessor
    expression = 'if(t > 1.5 * dt, d2x * (gx - d1x) + d2z * (gz - d1z), 0)'
    pp_names = 'gx gz d1x d1z d2x d2z'
    constant_names = 'dt'
    constant_expressions = '${dt}'
    use_t = true
    execute_on = 'initial timestep_end'
  []
  [den_step]
    type = ParsedPostprocessor
    expression = 'if(t > 1.5 * dt, 2 * ((gx - d1x)^2 + (gz - d1z)^2), 0)'
    pp_names = 'gx gz d1x d1z'
    constant_names = 'dt'
    constant_expressions = '${dt}'
    use_t = true
    execute_on = 'initial timestep_end'
  []
  [num]
    type = CumulativeValuePostprocessor
    postprocessor = num_step
  []
  [den]
    type = CumulativeValuePostprocessor
    postprocessor = den_step
  []
  [omega_raw]
    type = ParsedPostprocessor
    expression = 'acos(1 + num / den) / dt'
    pp_names = 'num den'
    constant_names = 'dt'
    constant_expressions = '${dt}'
  []
  [omega_ringdown]
    type = ParsedPostprocessor
    expression = '2 / dt * tan(omega_raw * dt / 2)'
    pp_names = 'omega_raw'
    constant_names = 'dt'
    constant_expressions = '${dt}'
  []
  [rel_raw]
    type = ParsedPostprocessor
    expression = 'omega_raw / abs(omega_eigen) - 1'
    pp_names = 'omega_raw omega_eigen'
  []
  [rel_ringdown]
    type = ParsedPostprocessor
    expression = 'omega_ringdown / abs(omega_eigen) - 1'
    pp_names = 'omega_ringdown omega_eigen'
  []
[]

[Outputs]
  file_base = out_sw_eigen_vs_ringdown
  csv = true
[]
