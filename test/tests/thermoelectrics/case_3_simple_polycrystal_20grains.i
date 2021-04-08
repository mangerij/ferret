[Mesh]
  [./basic_mesh]
  type = FileMeshGenerator
  file = 20g.e
 []


  [./add_sidesets]
    type = SideSetsFromNormalsGenerator
    input = basic_mesh
    normals = '1  0  0
              -1  0  0
               0  1  0
               0 -1  0
               0  0  1
               0  0 -1'
    fixed_normal = true
    new_boundary = 'right left front back top bottom'
    variance = 0.5
  [../]
[]

[Variables]
  [./potential_E_int]
    order = FIRST
    family = LAGRANGE
    # initial_condition = 1
  [../]
  [./T]
    order = FIRST
    family = LAGRANGE
    initial_condition = 25
  [../]
  []

[Kernels]
  [./residualV_x]
  type = ResidualV
  component = 0
  variable = potential_E_int
  T = 'T'
  potential_E_int = 'potential_E_int'
  ecC = 'ecC'
  sbC = 'sbC'
[../]

[./residualT_x]
 type = ResidualT
 component = 0
 variable = T
 T = 'T'
 potential_E_int = 'potential_E_int'
 ecC = 'ecC'
 sbC = 'sbC'
 thC = 'thC'
[../]

[]

[AuxVariables]
  [./j_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./j_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./j_z]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./q_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./q_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./q_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
  []

[AuxKernels]
  [./Electric_flux_x]
    type = Electric_flux
    variable = j_x
    T = T
    ecC = 'ecC'
    sbC = 'sbC'
    potential_E_int = 'potential_E_int'
    component = 0
  [../]
  [./Electric_flux_y]
    type = Electric_flux
    variable = j_y
    T = T
    ecC = 'ecC'
    sbC = 'sbC'
    potential_E_int = 'potential_E_int'
    component = 1
  [../]
  [./Electric_flux_z]
    type = Electric_flux
    variable = j_z
    T = T
    ecC = 'ecC'
    sbC = 'sbC'
    potential_E_int = 'potential_E_int'
    component = 2
  [../]

  [./Heat_flux_x]
    type = Heat_flux
    variable = q_x
    T = 'T'
    thC = 'thC'
    ecC = 'ecC'
    sbC = 'sbC'
    potential_E_int = 'potential_E_int'
    component = 0
  [../]
  [./heat_flux_y]
    type = Heat_flux
    variable = q_y
    T = 'T'
    thC = 'thC'
    ecC = 'ecC'
    sbC = 'sbC'
    potential_E_int = 'potential_E_int'
    component = 1
  [../]
  [./Heat_flux_z]
    type = Heat_flux
    variable = q_z
    T = 'T'
    thC = 'thC'
    ecC = 'ecC'
    sbC = 'sbC'
    potential_E_int = 'potential_E_int'
    component = 2
  [../]
  []


  [Materials]
    # same sbC
   #  [./ThermoelectricProperties1]
   #   type = GenericConstantMaterial
   #   prop_names = 'ecC sbC thC'
   #   prop_values = '8.422e4 1.941e-4 1.612'
   #   block = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20'
   # [../]
   # different sbC
   [./ThermoelectricProperties_1]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '41038.32138169937 -9.458012562559781e-05 0.7854876996829658'
    block = '1'
  [../]
  [./ThermoelectricProperties_2]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '81333.58011592268 -0.00018744773094871284 1.5567529226652503'
    block = '2'
  [../]
  [./ThermoelectricProperties_3]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '70204.40361474919 -0.00016179856021874634 1.343736625824931'
    block = '3'
  [../]
  [./ThermoelectricProperties_4]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '68213.5029818504 -0.0001572101748845543 1.3056300974441089'
    block = '4'
  [../]
  [./ThermoelectricProperties_5]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '59504.277875304295 -0.00013713821343619762 1.1389325093207139'
    block = '5'
  [../]
  [./ThermoelectricProperties_6]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '44718.36549954493 -0.00010306144316625114 0.8559250200102878'
    block = '6'
  [../]
  [./ThermoelectricProperties_7]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '28307.619070052922 -6.523995323554111e-05 0.5418176435635872'
    block = '7'
  [../]
  [./ThermoelectricProperties_8]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '11621.8224450934 -2.6784561108912714e-05 0.22244571101271146'
    block = '8'
  [../]
  [./ThermoelectricProperties_9]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '37598.99109971035 -8.66535760205863e-05 0.7196577256320718'
    block = '9'
  [../]
  [./ThermoelectricProperties_10]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '12642.07519589958 -2.9135915406365572e-05 0.24197370239598817'
    block = '10'
  [../]
  [./ThermoelectricProperties_11]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '16543.716284239385 -3.812794265935484e-05 0.316652465568676'
    block = '11'
  [../]
  [./ThermoelectricProperties_12]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '59052.90704234521 -0.00013609794890666357 1.1302931150826463'
    block = '12'
  [../]
  [./ThermoelectricProperties_13]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '10142.219112025681 -2.3374551527477853e-05 0.19412559022305154'
    block = '13'
  [../]
  [./ThermoelectricProperties_14]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '5586.057883392847 -1.2874065960182278e-05 0.10691908463582604'
    block = '14'
  [../]
  [./ThermoelectricProperties_15]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '45166.65802963037 -0.00010409461319818636 0.8645054944640722'
    block = '15'
  [../]
  [./ThermoelectricProperties_16]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '64085.0484929901 -0.00014769541572654213 1.2266100471467591'
    block = '16'
  [../]
  [./ThermoelectricProperties_17]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '31084.02733134581 -7.163868089544314e-05 0.5949590602960039'
    block = '17'
  [../]
  [./ThermoelectricProperties_18]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '23020.612595655355 -5.3055104545437004e-05 0.44062250658034235'
    block = '18'
  [../]
  [./ThermoelectricProperties_19]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '38253.31651101918 -8.81615855472432e-05 0.7321817408663372'
    block = '19'
  [../]
  [./ThermoelectricProperties_20]
    type = GenericConstantMaterial
    prop_names = 'ecC sbC thC'
    prop_values = '12804.494846512856 -2.9510240438234925e-05 0.24508247082140497'
    block = '20'
  [../]
  []


[BCs]
  [./sideset_1T]
    type = DirichletBC
    variable = T
    boundary = 'left'
    value = 298
  [../]
  [./sideset_2T]
    type = DirichletBC
    variable = T
    boundary = 'right'
    value = 298
  [../]
  [./sideset_1V]
    type = DirichletBC
    variable = potential_E_int
    boundary = 'right'
    value = 0

  [../]
  [./sideset_2V]
    type = DirichletBC
    variable = potential_E_int
    boundary = 'left'
    value = 0.058
  [../]
  # [./Periodic]
  #   [./All]
  #     auto_direction = 'x y z'
  #   [../]
  # [../]
  []

[Executioner]
  type = Steady
  solve_type = NEWTON
  # line_search = none
[]

[Preconditioning]
[./SMP]
type = SMP
full = true
[../]
[]

[Postprocessors]
  [./potential_x]
    type = PointValue
    point = '0 0 0'
    variable = potential_E_int
  [../]

  [./T]
    type = PointValue
    point = '0 0 0'
    variable = T
  [../]
  []

[Outputs]
  exodus = true
  csv = true
  file_base = polycrystal_20grains
[]

[Debug]
  show_var_residual_norms = true
  []
