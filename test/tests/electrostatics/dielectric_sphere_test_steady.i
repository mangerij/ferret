[Mesh]
  file = sphere_in_box_nr3.e
[]

[Variables]
  [./potential]
    order=FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./E_Ext_block2]
     type = Electrostatics
     variable = potential
     block = '2'
  [../]

  [./E_Ext_block1]
     type = Electrostatics
     variable = potential
     block = '1'
  [../]
[]

[Materials]
  [./permitivitty_1]
    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '1'
    block = '1'
  [../]
  [./permitivitty_2]
    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '6'
    block = '2'
  [../]
[]

[BCs]

  [./potential_ext_1]
    type = DirichletBC
    variable = potential
    boundary = '1'
    value = 62.0
  [../]
  [./potential_ext_2]
    type = DirichletBC
    variable = potential
    boundary = '2'
    value = -62.0
  [../]
[]

[Preconditioning]
   [./smp]
     type = SMP
     full = true
   [../]
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]

[Outputs]
  file_base = out_lindie_sph_test
  print_linear_residuals = true
  [./out]
    type = Exodus
    elemental_as_nodal = true
  [../]
[]
