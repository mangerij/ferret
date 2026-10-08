[Mesh]
  # muMAG standard problem 2 (5d x d x 0.1d film, no anisotropy, field along [111]), FIRST-order demag
  # potential. Reduced units: l_ex = 1, Km = mu0 Ms^2 / 2 = 1; lengths below are for d = 1 l_ex.
  # Scaled-down mesh for a fast regression test. For the actual SP2 benchmark multiply every length by
  # d/l_ex and use 100 x 20 x 3 film cells with a graded vacuum, e.g.
  # dx = '2 1 0.5 0.25 0.25 5 0.25 0.25 0.5 1 2', ix = '1 1 1 1 5 100 5 1 1 1 1',
  # dy = '2 1 0.5 0.25 0.25 1 0.25 0.25 0.5 1 2', iy = '1 1 1 1 5 20 5 1 1 1 1',
  # dz = '1.2 0.6 0.3 0.133333 0.066667 0.033333 0.1 0.033333 0.066667 0.133333 0.3 0.6 1.2',
  # iz = '1 1 1 1 1 1 3 1 1 1 1 1 1' (box half-sizes 6.5d x 4.5d x 2.383333d), then sweep the field.
  [box]
    type = CartesianMeshGenerator
    dim = 3
    dx = '4.0 1.0 5.0 1.0 4.0'
    dy = '4.0 1.0 1.0 1.0 4.0'
    dz = '2.0 0.2 0.1 0.2 2.0'
    ix = '1 1 10 1 1'
    iy = '1 1 2 1 1'
    iz = '1 1 1 1 1'
  []
  [center]
    type = TransformGenerator
    input = box
    transform = TRANSLATE
    vector_value = '-7.5 -5.5 -2.25'
  []
  [film_block]
    type = ParsedSubdomainMeshGenerator
    input = center
    block_id = 1
    combinatorial_geometry = 'abs(x) <= 2.500001 & abs(y) <= 0.500001 & abs(z) <= 0.050001'
  []
  [block_names]
    type = RenameBlockGenerator
    input = film_block
    old_block = '0 1'
    new_block = 'vacuum film'
  []
[]

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z
  potential_H_int = potential_H_int
  Hext_x = Hext
  Hext_y = Hext
  Hext_z = Hext
  mu0 = 1.256637
  g0 = 221.1
  Hscale = 1.0
[]

[Variables]
  # start uniform along [111]
  [mag_x]
    block = film
    [InitialCondition]
      type = ConstantIC
      value = 0.5773502691896258
    []
  []
  [mag_y]
    block = film
    [InitialCondition]
      type = ConstantIC
      value = 0.5773502691896258
    []
  []
  [mag_z]
    block = film
    [InitialCondition]
      type = ConstantIC
      value = 0.5773502691896258
    []
  []
  [potential_H_int]
    block = 'film vacuum'
  []
[]

[AuxVariables]
  [mag_s]
    block = film
  []
  # applied field H/Ms = 0.08 along [111]: each component = 0.08 * Ms / sqrt(3)
  [Hext]
    block = film
    [InitialCondition]
      type = ConstantIC
      value = 0.05827295
    []
  []
[]

[AuxKernels]
  [mag_norm]
    type = VectorMagnitudeAux
    variable = mag_s
    x = mag_x
    y = mag_y
    z = mag_z
    block = film
    execute_on = 'initial timestep_end'
  []
[]

[Materials]
  [film_properties]
    type = GenericConstantMaterial
    prop_names = 'alpha Ae Ms permittivity'
    prop_values = '1.0 1.0 1.261566 1.0'
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
    type = TimeDerivative
    variable = mag_x
    block = film
  []
  [mag_y_time]
    type = TimeDerivative
    variable = mag_y
    block = film
  []
  [mag_z_time]
    type = TimeDerivative
    variable = mag_z
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
  [exchange_z]
    type = MasterExchangeCartLLG
    variable = mag_z
    component = 2
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
  [field_z]
    type = MasterInteractionCartLLGHConst
    variable = mag_z
    component = 2
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
    block = film
  []
  [my]
    type = ElementAverageValue
    variable = mag_y
    block = film
  []
  [mz]
    type = ElementAverageValue
    variable = mag_z
    block = film
  []
  [mag_norm_min]
    type = NodalExtremeValue
    variable = mag_s
    value_type = min
    block = film
  []
  [Fexch]
    type = MasterMagneticExchangeEnergy
    energy_scale = 1.0
    block = film
  []
  [Fdemag]
    type = MagnetostaticEnergyCart
    energy_scale = 1.0
    block = film
  []
  [Fzeeman]
    type = MasterMagneticZeemanEnergyCart
    energy_scale = 1.0
    block = film
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
  dt = 0.01
  num_steps = 2
[]

[Outputs]
  file_base = out_sp2_first
  exodus = true
  csv = true
[]
