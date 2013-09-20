# This input file generate vector values on the mesh cell. The result is saved in ExodusII. It can then be read in to another appliction for cell-dependent value.
[Mesh]
  file = poissonstripe_coarse.e
 # uniform_refine=2
[]
[Variables]
  [./dummy]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./dummy]
    type = Diffusion
    variable = dummy
  [../]
[]

[AuxVariables]
  [./polar_angle]
    order = FIRST
    family = LAGRANGE
  [../]
  [./azimuthal_angle]
    order = FIRST
    family = LAGRANGE
  [../]
  [./radial]
    order = FIRST
    family = LAGRANGE
  [../]
  # [./polar_x]
  #    order = FIRST
  #   family = LAGRANGE
  # [../]
  #  [./polar_y]
  #    order = FIRST
  #   family = LAGRANGE
  # [../]
  # [./polar_z]
  #    order = FIRST
  #   family = LAGRANGE
  # [../]
[]

[AuxKernels]
  #active='diff'
  [./polar_angle]
    type=FunctionAux
    variable = polar_angle
    function=polar_angle
  [../]
  [./azimuthal_angle]
    type=FunctionAux
    variable = azimuthal_angle
    function=azimuthal_angle
  [../]
  [./radial]
    type=FunctionAux
    variable=radial
    function=radial
  [../]
  # [./polar_x]
  # type=FunctionAux
  # variable=polar_x
  # function=polar_x
  # [../]
  # [./polar_y]
  # type=FunctionAux
  # variable=polar_y
  # function=polar_y
  # [../]
  # [./polar_z]
  # type=FunctionAux
  # variable=polar_z
  # function=polar_z
  # [../]
[]

[Executioner]
  type = Steady
[]

[Functions]
  # [./polar_angle]
  #  type=SinFunc
  #  amplitude=3.14
  #  wave_length_x=1.0
  #  wave_length_y=1.0
  #  wave_length_z=0.2
  #  phrase_x=1.57
  #  phrase_y=1.57
  #  phrase_z=1.57
  #  vertical_shift=0
  # [../]
  # [./azimuthal_angle]
  #  type=SinFunc
  #  amplitude=1.57
  #  wave_length_x=1.0
  #  wave_length_y=1.0
  #  wave_length_z=0.2
  #  phrase_x=1.57
  #  phrase_y=1.57
  #  phrase_z=1.57
  #  vertical_shift=0
  # [../]
  [./polar_angle]
  type=RandomFunc
  min=0
  max=6.28
  [../]
  [./azimuthal_angle]
  type=RandomFunc
  min=-0.78
  max=0.78
  [../]
  [./radial]
  type=RandomFunc
  min=0.3
  max=0.8
  [../]
  # [./polar_x]
  # type=SphereToCartFunc
  # radial_function=radial
  # polar_function=polar_angle
  # azimuthal_function=azimuthal_angle
  # index=0
  # [../]
  # [./polar_y]
  # type=SphereToCartFunc
  # radial_function=radial
  # polar_function=polar_angle
  # azimuthal_function=azimuthal_angle
  # index=1
  # [../]
  # [./polar_z]
  # type=SphereToCartFunc
  # radial_function=radial
  # polar_function=polar_angle
  # azimuthal_function=azimuthal_angle
  # index=2
  # [../]
[]

[Output]
  file_base = initvalues_random
  exodus = true
  elemental_as_nodal=true
[]
