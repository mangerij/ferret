[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 50
  ny = 50
  nz = 50

  xmax = 1
  ymax = 1
  zmax = 1
  elem_type =HEX8
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
    active = 'manual_x manual_y'

    # Can use auto_direction with Generated Meshes
    # [./auto]
    #   variable = u
    #   auto_direction = 'x y'
    # [../]

     # Use Translation vectors for everything else
     [./manual_x]
       variable = u
       primary = 'left'
       secondary = 'right'
       translation = '1 0 0'
     [../]

     [./manual_y]
       variable = u
       primary = 'bottom'
       secondary = 'top'
       translation = '0 1 0'
     [../]
  [../]
[./bottom]
     type=DirichletBC
     variable=u
     boundary='top'
     value=0
  [../]
  [./top]
     type=DirichletBC
     variable=u
     boundary='bottom'
     value=0
  [../]
[]

[Executioner]
   type=Steady
[]

[Output]
  file_base = out_pbc1
  interval = 1
  exodus = true
  perf_log = true
[]

