[Mesh]
  # Scaled-down mesh for a fast regression test. For the actual muMAG SP1 benchmark use 20 nm
  # in-plane film cells and a graded vacuum, e.g. ix = '8 6 50 6 8', iy = '8 6 100 6 8',
  # iz = '8 6 2 6 8' with dz = '300.0 7.0 3.517533 7.0 300.0'.
  [box]
    type = CartesianMeshGenerator
    dim = 3
    dx = '500.0 35.2 175.876674 35.2 500.0'
    dy = '500.0 35.2 351.753349 35.2 500.0'
    dz = '300.0 7.0 3.517533 7.0 300.0'
    ix = '2 1 5 1 2'
    iy = '2 1 10 1 2'
    iz = '2 1 1 1 2'
  []
  [center]
    type = TransformGenerator
    input = box
    transform = TRANSLATE
    vector_value = '-623.138337 -711.076675 -308.758767'
  []
  [film_block]
    type = ParsedSubdomainMeshGenerator
    input = center
    block_id = 1
    combinatorial_geometry = 'x >= -87.938337 & x <= 87.938337 & y >= -175.876675 & y <= 175.876675 & z >= -1.758767 & z <= 1.758767'
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
  Hext_x = Hext_x
  Hext_y = Hext_y
  Hext_z = Hext_z
  mu0 = 1.256637
  g0 = 140.206695
  Hscale = 1.0
[]

[Variables]
  [mag_x]
    block = film
    [InitialCondition]
      type = ConstantIC
      value = 0.0
    []
  []
  [mag_y]
    block = film
    [InitialCondition]
      type = ConstantIC
      value = 1.0
    []
  []
  [mag_z]
    block = film
    [InitialCondition]
      type = ConstantIC
      value = 0.0
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
  [Hext_x]
    block = film
    [InitialCondition]
      type = ConstantIC
      value = 0.06273559
    []
  []
  [Hext_y]
    block = film
    [InitialCondition]
      type = ConstantIC
      value = 0.00109505
    []
  []
  [Hext_z]
    block = film
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
    prop_names = 'alpha Ae Ms permittivity K1 nx ny nz'
    prop_values = '1.0 1.0 1.261566 1.0 -0.001243 0 1 0'
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
  [anisotropy_x]
    type = MasterAnisotropyCartLLG
    variable = mag_x
    component = 0
    block = film
  []
  [anisotropy_y]
    type = MasterAnisotropyCartLLG
    variable = mag_y
    component = 1
    block = film
  []
  [anisotropy_z]
    type = MasterAnisotropyCartLLG
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
  [Faniso]
    type = MasterMagneticAnisotropyEnergy
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
  dt = 0.05
  num_steps = 5
[]

[Outputs]
  file_base = out_sp1
  exodus = true
  csv = true
[]
