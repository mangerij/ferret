!include film_common.i
!include transient_common.i

amplitude = 0.01
dt = 0.01
num_steps = 50

[Functions]
  [ic_mag_x]
    type = ParsedFunction
    expression = '${amplitude} * cos(pi * x / 40.0)'
  []
  [ic_mag_z]
    type = ParsedFunction
    expression = '0.0'
  []
  [ic_potential]
    type = ParsedFunction
    expression = '0.0'
  []
[]

[Postprocessors]
  [gx_quarter]
    type = PointValue
    variable = mag_x
    point = '10.0 0.5 0.0'
    execute_on = 'initial timestep_end'
  []
  [gz_quarter]
    type = PointValue
    variable = mag_z
    point = '10.0 0.5 0.0'
    execute_on = 'initial timestep_end'
  []
[]

[Outputs]
  file_base = out_sw_ringdown_film
  csv = true
[]
