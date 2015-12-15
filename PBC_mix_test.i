#debug file to find how to implement PBC at end caps of cylinder
# translation vector doesnt work for regular quasi-2D pbc but auto_direction does
# for irregular auto_direction complains about being irregular so try translation
# should check auto_direction for quasi-2D for a non-generated mesh (perhaps a cheap superlattice)


[Mesh]
  file = exodus_cylinder_r3_h20.e
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
  [./dot]
    type = TimeDerivative
    variable = u
  [../]
[]

[BCs]
  [./Periodic]
   [./u]
    variable = u
    primary = '2'
    secondary = '3'
    translation = '0 0 10'
   [../]
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.5
  num_steps = 10
  solve_type = NEWTON
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]
