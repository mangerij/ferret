#This is a polarized (spontaneously) sphere in a dielectric medium.
#A background (optical) dielectric constant of 1 is assigned, while the medium has dielectric constant of 10.
#The nonzero polarization gives rise to nonzero surface charges which are handled naturally.

[Mesh]
  file = exodus_sphere3.e
  #block 2  = sphere
  #block 1  = medium
[]

[GlobalParams]
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
[]


[Variables]
  [./potential]
    order = FIRST
    family = LAGRANGE
  [../]
[]


[AuxVariables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = ConstantIC
      value = 0.0
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = ConstantIC
      value = 0.5
    [../]
  [../]
  [./E_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./E_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./E_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./E_Ext_block2]
     type = Electrostatics
     variable = potential
     block='1'
  [../]

  [./E_Ext_block1]
     type = Electrostatics
     variable = potential
     block = '2'
  [../]

  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential
     block = '1'
  [../]
[]

[AuxKernels]
  [./Ex]
    type = QuasistaticFieldAux
    variable = E_x
    block = '1 2'
    component = 0
    potential_int = potential
  [../]
  [./Ey]
    type = QuasistaticFieldAux
    variable = E_y
    block = '1 2'
    component = 1
    potential_int = potential
  [../]
  [./Ez]
    type = QuasistaticFieldAux
    variable = E_z
    block = '1 2'
    component = 2
    potential_int = potential
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
    prop_values = '10'
    block = '2'
  [../]
[]

[BCs]
  [./potential_ext_1]
    type = DirichletBC
    variable = potential
    boundary = '1 2 3 4 5 6'
    value = 0.0
  [../]
[]

[Preconditioning]
   [./smp]
     type = SMP
     full = true   #to use every off diagonal block
   [../]
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_polar_sph_test
    elemental_as_nodal = true
    interval = 1
  [../]
[]
