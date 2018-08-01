[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 9
  ny = 9
  nz = 9
  xmin = -4
  xmax = 4
  ymin = -4
  ymax = 4
  zmin = -4
  zmax = 4
  elem_type = HEX8
[]

[GlobalParams]
  len_scale = 1.0
  alpha1 = -0.1722883
  alpha11 = -0.07253
  alpha111 = 0.26
  alpha12 = 0.75
  alpha112 = 0.61
  alpha123 = -3.67

  G110 = 0.15
  G11_G110 = 0.6
  G12_G110 = 0
  G44_G110 = 0.3
  G44P_G110 = 0.3

  permittivity = 0.5843763

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_E_int = potential_E_int
[]



[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-6
      max = 0.01e-6
      seed = 5
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-6
      max = 0.01e-6
      seed = 5
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-6
      max = 0.01e-6
      seed = 5
    [../]
  [../]
  [./potential_E_int]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  #Bulk energy density
  [./bed_x]
    type = BulkEnergyDerivativeSixth
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = BulkEnergyDerivativeSixth
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = BulkEnergyDerivativeSixth
    variable = polar_z
    component = 2
  [../]
  ##Wall energy penalty
  [./walled_x]
     type=WallEnergyDerivative
     variable = polar_x
     component = 0
  [../]
  [./walled_y]
     type=WallEnergyDerivative
     variable = polar_y
     component = 1
  [../]
  [./walled_z]
     type=WallEnergyDerivative
     variable = polar_z
     component = 2
  [../]

  ##Electrostatics
  [./polar_x_electric_E]
     type=PolarElectricEStrong
     variable = potential_E_int
     block = '0'
  [../]
  [./FE_E_int]
     type=Electrostatics
     variable = potential_E_int
     block = '0'
  [../]
  [./polar_electric_px]
     type=PolarElectricPStrong
     variable = polar_x
     component = 0
  [../]
  [./polar_electric_py]
     type=PolarElectricPStrong
     variable = polar_y
     component = 1
  [../]
  [./polar_electric_pz]
     type=PolarElectricPStrong
     variable = polar_z
     component = 2
  [../]
  ##Time dependence
  [./polar_x_time]
     type=TimeDerivativeScaled
     variable=polar_x
    time_scale = 1.0
  [../]
  [./polar_y_time]
     type=TimeDerivativeScaled
     variable=polar_y
    time_scale = 1.0
  [../]
  [./polar_z_time]
     type=TimeDerivativeScaled
     variable = polar_z
    time_scale = 1.0
  [../]
[]

[BCs]
  [./potential_cube5]
    type = DirichletBC
    boundary = 'front'
    value = 0.0002
    variable = potential_E_int
  [../]
  [./potential_cube6]
    type = DirichletBC
    boundary = 'back'
    value = 0.0002
    variable = potential_E_int
  [../]

[]


[Postprocessors]
  [./bulk_energy]
   type = BulkEnergy
   execute_on = 'initial timestep_end'
  [../]
  [./wall_energy]
   type = WallEnergy
   execute_on = 'initial timestep_end'
  [../]
  [./electrostatic_energy]
   type = ElectrostaticEnergy
   execute_on = 'initial timestep_end'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121               1e-10       1e-8       1e-6     bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 0.5
  num_steps = 8
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_wall_test
  [../]
[]
