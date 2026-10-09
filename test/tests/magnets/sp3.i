[Mesh]
  # Scaled-down mesh for a fast regression test. For the actual muMAG SP3 benchmark use at least
  # N = 10 cells per cube edge (h <= 0.1 L, extrapolate in N), a refined vacuum shell >= 0.3 L
  # thick, and a graded vacuum out to a box >= 10 L.
  [box]
    type = CartesianMeshGenerator
    dim = 3
    dx = '16.8 2.1 8.4 2.1 16.8'
    dy = '16.8 2.1 8.4 2.1 16.8'
    dz = '16.8 2.1 8.4 2.1 16.8'
    ix = '2 1 4 1 2'
    iy = '2 1 4 1 2'
    iz = '2 1 4 1 2'
  []
  [center]
    type = TransformGenerator
    input = box
    transform = TRANSLATE
    vector_value = '-23.1 -23.1 -23.1'
  []
  [cube_block]
    type = SubdomainBoundingBoxGenerator
    input = center
    block_id = 1
    bottom_left = '-4.2 -4.2 -4.2'
    top_right = '4.2 4.2 4.2'
  []
  [block_names]
    type = RenameBlockGenerator
    input = cube_block
    old_block = '0 1'
    new_block = 'vacuum cube'
  []
[]

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z
  potential_H_int = potential_H_int
  mu0 = 1.256637
  g0 = 221.1
  Hscale = 1.0
[]

[Functions]
  [flower_mx]
    type = ParsedFunction
    expression = '1/sqrt(1 + (0.45*x*y/17.64)^2 + (0.45*x*z/17.64)^2)'
  []
  [flower_my]
    type = ParsedFunction
    expression = '(0.45*x*y/17.64)/sqrt(1 + (0.45*x*y/17.64)^2 + (0.45*x*z/17.64)^2)'
  []
  [flower_mz]
    type = ParsedFunction
    expression = '(0.45*x*z/17.64)/sqrt(1 + (0.45*x*y/17.64)^2 + (0.45*x*z/17.64)^2)'
  []
[]

[Variables]
  [mag_x]
    block = cube
    [InitialCondition]
      type = FunctionIC
      function = flower_mx
    []
  []
  [mag_y]
    block = cube
    [InitialCondition]
      type = FunctionIC
      function = flower_my
    []
  []
  [mag_z]
    block = cube
    [InitialCondition]
      type = FunctionIC
      function = flower_mz
    []
  []
  [potential_H_int]
    block = 'cube vacuum'
  []
[]

[AuxVariables]
  [mag_s]
    block = cube
  []
[]

[AuxKernels]
  [mag_norm]
    type = VectorMagnitudeAux
    variable = mag_s
    x = mag_x
    y = mag_y
    z = mag_z
    block = cube
    execute_on = 'initial timestep_end'
  []
[]

[Materials]
  [cube_properties]
    type = GenericConstantMaterial
    prop_names = 'alpha Ae Ms permittivity K1 nx ny nz'
    prop_values = '100.0 1.0 1.261566 1.0 -0.1 1 0 0'
    block = cube
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
    type = TimeDerivative
    variable = mag_x
    block = cube
  []
  [mag_y_time]
    type = TimeDerivative
    variable = mag_y
    block = cube
  []
  [mag_z_time]
    type = TimeDerivative
    variable = mag_z
    block = cube
  []
  [exchange_x]
    type = MasterExchangeCartLLG
    variable = mag_x
    component = 0
    block = cube
  []
  [exchange_y]
    type = MasterExchangeCartLLG
    variable = mag_y
    component = 1
    block = cube
  []
  [exchange_z]
    type = MasterExchangeCartLLG
    variable = mag_z
    component = 2
    block = cube
  []
  [anisotropy_x]
    type = MasterAnisotropyCartLLG
    variable = mag_x
    component = 0
    block = cube
  []
  [anisotropy_y]
    type = MasterAnisotropyCartLLG
    variable = mag_y
    component = 1
    block = cube
  []
  [anisotropy_z]
    type = MasterAnisotropyCartLLG
    variable = mag_z
    component = 2
    block = cube
  []
  [demag_x]
    type = MasterInteractionCartLLG
    variable = mag_x
    component = 0
    block = cube
  []
  [demag_y]
    type = MasterInteractionCartLLG
    variable = mag_y
    component = 1
    block = cube
  []
  [demag_z]
    type = MasterInteractionCartLLG
    variable = mag_z
    component = 2
    block = cube
  []
  [magnetostatic_laplace]
    type = Electrostatics
    variable = potential_H_int
    block = 'cube vacuum'
  []
  [magnetostatic_source]
    type = MagHStrongCart
    variable = potential_H_int
    block = cube
  []
[]

[BCs]
  [outer_vacuum]
    type = DirichletBC
    variable = potential_H_int
    boundary = '0 1 2 3 4 5'
    value = 0.0
  []
[]

[Postprocessors]
  [mx]
    type = ElementAverageValue
    variable = mag_x
    block = cube
  []
  [my]
    type = ElementAverageValue
    variable = mag_y
    block = cube
  []
  [mz]
    type = ElementAverageValue
    variable = mag_z
    block = cube
  []
  [mag_norm_min]
    type = NodalExtremeValue
    variable = mag_s
    value_type = min
    block = cube
  []
  [Fexch]
    type = MasterMagneticExchangeEnergy
    energy_scale = 1.0
    block = cube
  []
  [Faniso]
    type = MasterMagneticAnisotropyEnergy
    energy_scale = 1.0
    block = cube
  []
  [Fdemag]
    type = MagnetostaticEnergyCart
    energy_scale = 1.0
    block = cube
  []
[]

[UserObjects]
  [renormalize]
    type = PointwiseRenormalizeVector
    v = 'mag_x mag_y mag_z'
    execute_on = timestep_end
    force_preaux = true
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = implicit-euler
  petsc_options_iname = '-pc_type -snes_atol -snes_rtol'
  petsc_options_value = 'lu       1e-8       1e-8'
  dt = 0.25
  num_steps = 5
[]

[Outputs]
  file_base = out_sp3
  exodus = true
  csv = true
[]
