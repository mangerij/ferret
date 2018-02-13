
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 7
  ny = 7
  nz = 5
  xmin = -3
  xmax = 3
  ymin = -3
  ymax = 3
  zmin = -2
  zmax = 2
  elem_type = HEX8
[]

[GlobalParams]
  len_scale = 1.0
  G110 = 0.173
  G11_G110 = 0.6
  G12_G110 = 0
  G44_G110 = 0.3
  G44P_G110 = 0.3

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
[]



[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-4
      max = 0.01e-4
      seed = 5
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-4
      max = 0.01e-4
      seed = 5
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-4
      max = 0.01e-4
      seed = 5
    [../]
  [../]
[]


[Materials]
  [./ferro_rank_two]
    type = ComputeRankTwoLandauTensor
    a_ij = '0 0 0 0 0 0 0 0 0'
  [../]
  [./ferro_rank_four]
    type = ComputeRankFourLandauTensor
    a_ijkl = '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0'
  [../]
[]


[Kernels]
  [./bed_x]
    type = LandauMaterialBulkFourthDerivative
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = LandauMaterialBulkFourthDerivative
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = LandauMaterialBulkFourthDerivative
    variable = polar_z
    component = 2
  [../]

  [./walled_x]
    type = WallEnergyDerivative
    variable = polar_x
    component = 0
 [../]
 [./walled_y]
    type = WallEnergyDerivative
    variable = polar_y
    component = 1
  [../]
  [./walled_z]
     type = WallEnergyDerivative
     variable = polar_z
     component = 2
  [../]

  [./polar_x_time]
     type = TimeDerivativeScaled
     variable=polar_x
     time_scale = 1.0
  [../]
  [./polar_y_time]
     type = TimeDerivativeScaled
     variable=polar_y
     time_scale = 1.0
  [../]
  [./polar_z_time]
     type = TimeDerivativeScaled
     variable = polar_z
     time_scale = 1.0
  [../]
[]


[BCs]

[]



[Postprocessors]

[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121               1e-10      1e-8      1e-6    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 0.1
  num_steps = 10
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = false
  [./out]
    type = Exodus
    file_base = outPTO_LandauMat_test
    elemental_as_nodal = true
    interval = 1
  [../]
[]
