[Mesh]
  file=poissonbox.e
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

  [./forcing]
    type = Magnetizing
    variable = u
  [../]
[]

[BCs]
  [./Periodic]
    #Note: Enable either "auto" or both "manual" conditions for this example
    #active = 'manual_x manual_y'

    # Can use auto_direction with Generated Meshes
    # [./auto]
    #   variable = u
    #   auto_direction = 'x y'
    # [../]

     # Use Translation vectors for everything else
     [./manual_x]
       variable = u
       primary = 4
       secondary = 6
       translation = '1 0 0'
     [../]		       

     [./manual_y]
       variable = u
       primary = 3
       secondary = 5
       translation = '0 1 0'
     [../]
  [../]
  [./bottom]
     type=DirichletBC
     variable=u
     boundary=2
     value=0
  [../]
  [./top]
     type=DirichletBC
     variable=u
     boundary=1
     value=0
  [../]
[]


[Executioner]
  type = Steady
[]

[Output]
  file_base = out_pbc_hex
  interval = 1
  exodus = true
  perf_log = true
[]

