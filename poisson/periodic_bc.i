[Mesh]
  # type = GeneratedMesh
  # dim = 3
  # nx = 50
  # ny = 50
  # nz = 50

  # xmax = 1
  # ymax = 1
  # zmax = 1
  # elem_type = QUAD4
  # file=brick_hex.e
  #file=poissonbox4.e
   file=poissonstripe.e
  #file=brick.e
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
     [./auto]
       variable = u
       auto_direction = 'x y'
     [../]

     # Use Translation vectors for everything else
     #[./manual_x]
     #  variable = u
     #  primary = 4
     #  secondary = 6
     #  translation = '1 0 0'
     #[../]		       

     #[./manual_y]
     #  variable = u
     #  primary = 3
     #  secondary = 5
     #  translation = '0 1 0'
     #[../]
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

[Materials]
  [./Exterior]
      type=PolarMaterial
      block='exterior'
      P='0 0 0'  
  [../]
  [./Interior]
      type=PolarMaterial
      block='interior'
      P='0 0 1'
  [../]
[]

[Executioner]
  type = Steady
[]

[Output]
  #file_base = out_pbc_hex
  file_base = out_pbc_stripe
  interval = 1
  exodus = true
  perf_log = true
[]

