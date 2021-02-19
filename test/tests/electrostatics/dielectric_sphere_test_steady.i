# This is your standard linear dielectric (polarizable) sphere in a dielectric medium
# This problem is worked out analytically in Jackson or Griffiths and we've shown that
# the analytical solution matches the numerical solution to within high degree of accuracy
# provided the medium mesh is large enough.

[Mesh]
  file = exodus_sphere3.e
  #block 2  = sphere
  #block 1  = medium
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
  # Please note that in typical textbooks, you apply the field with a Neumann boundary condition
  # but this isn't ideal for the exodiff so we can just as easily apply with a DirichletBC
  # but also note that the field strength depends on the distance between the boundaries!

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
