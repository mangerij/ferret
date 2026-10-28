[Mesh]
  [box]
    type = CartesianMeshGenerator
    dim = 2
    dx = '40.0'
    dy = '80.0 40.0 20.0 10.0 5.0 2.5 1.25 0.5 1.0 0.5 1.25 2.5 5.0 10.0 20.0 40.0 80.0'
    ix = '20'
    iy = '1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1'
    subdomain_id = '0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0'
  []
  [center]
    type = TransformGenerator
    input = box
    transform = TRANSLATE
    vector_value = '0.0 -159.75 0.0'
  []
  [block_names]
    type = RenameBlockGenerator
    input = center
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

[AuxVariables]
  [mag_y]
    block = film
    [InitialCondition]
      type = ConstantIC
      value = 1.0
    []
  []
  [Hext_x]
    block = film
  []
  [Hext_y]
    block = film
    [InitialCondition]
      type = ConstantIC
      value = 0.238434
    []
  []
  [Hext_z]
    block = film
  []
[]

[Materials]
  [film_properties]
    type = GenericConstantMaterial
    prop_names = 'alpha Ae Ms permittivity'
    prop_values = '0.0 1.0 1.261566 1.0'
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
  [exchange_x]
    type = MasterExchangeCartLLG
    variable = mag_x
    component = 0
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
[]

[BCs]
  [outer_vacuum]
    type = DirichletBC
    variable = potential_H_int
    boundary = 'left right bottom top'
    value = 0.0
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]
