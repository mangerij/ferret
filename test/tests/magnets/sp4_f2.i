[Mesh]
  # muMAG SP4 field 2, mesh scaled down for the test. Real mesh used: 2.5 nm film cells (200 x 50 in-plane,
  # 2 through the 3 nm thickness), refined to 1.25 nm (4 through the thickness) in a 5 nm band along the film
  # edges, and a graded hanging-node vacuum out to a 185 x 152 x 180 l_ex box
  [box]
    type = CartesianMeshGenerator
    dim = 3
    dx = '100.0 20.0 87.938337 20.0 100.0'
    dy = '100.0 20.0 21.984584 20.0 100.0'
    dz = '40.0 4.0 0.527630 4.0 40.0'
    ix = '2 1 16 1 2'
    iy = '2 1 4 1 2'
    iz = '2 1 1 1 2'
  []
  [center]
    type = TransformGenerator
    input = box
    transform = TRANSLATE
    vector_value = '-163.969169 -130.992292 -44.263815'
  []
  [film_block]
    type = ParsedSubdomainMeshGenerator
    input = center
    block_id = 1
    combinatorial_geometry = 'x >= -43.969169 & x <= 43.969169 & y >= -10.992292 & y <= 10.992292 & z >= -0.263815 & z <= 0.263815'
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
      value = 0.991756
    []
  []
  [mag_y]
    block = film
    [InitialCondition]
      type = ConstantIC
      value = 0.128140
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
      value = -0.0445490532
    []
  []
  [Hext_y]
    block = film
    [InitialCondition]
      type = ConstantIC
      value = -0.0079058883
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
    prop_names = 'alpha Ae Ms permittivity'
    prop_values = '0.02 1.0 1.261566 1.0'
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
  [mag_norm_max]
    type = NodalExtremeValue
    variable = mag_s
    value_type = max
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
  [TimeIntegrator]
    type = NewmarkBeta
  []
  petsc_options_iname = '-pc_type -snes_atol -snes_rtol'
  petsc_options_value = 'lu       1e-8       1e-8'
  dt = 0.001
  num_steps = 5
[]

[Outputs]
  file_base = out_sp4_f2
  exodus = true
  csv = true
[]
