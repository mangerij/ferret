
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 25
  ny = 25
  nz = 20
  xmin = -5
  xmax = 5
  ymin = -5
  ymax = 5
  zmin = -3
  zmax = 3
  elem_type = HEX8
[]



[GlobalParams]
  len_scale = 1.0

  alpha11 = 0.25
  alpha111 = 0.0
  alpha12 = 0.0
  alpha112 = 0.0
  alpha123 = 0.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  # -0.1722883 #room temp PTO
  alpha1 = -0.5

[]


[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.0001
      max = 0.0001
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.0001
      max = 0.0001
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.0001
      max = 0.0001
    [../]
  [../]
[]

[Kernels]
  [./bed_z]
    type = BulkEnergyDerivativeSixth
    variable = polar_z
    component = 2
  [../]

  [./const_field]
    type = ConstField
    variable = polar_z
    field = 0.1
  [../]

  ##Time dependence
  [./polar_x_time]
     type = TimeDerivativeScaled
     variable = polar_x
     time_scale = 1.0
  [../]
  [./polar_y_time]
     type = TimeDerivativeScaled
     variable = polar_y
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
  [./avePz]
    type = ElementAverageValue
    variable = polar_z
    execute_on = 'timestep_end'
  [../]
  [./Fbulk]
    type = BulkEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = Fbulk
  [../]
[]



[UserObjects]
 [./kill]
  type = Terminator
  expression = 'perc_change <= 1.0e-3'
 [../]
[]



[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    120               1e-10      1e-8     1e-4    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmax = 0.7
[]

[Outputs]
  print_linear_residuals = false
  print_perf_graph = true
  [./out]
    type = Exodus
    execute_on = 'timestep_end'
    file_base = out_prob0
    elemental_as_nodal = true
    interval = 1
  [../]
  [./outcsv]
    type = CSV
    file_base = out_prob0
    execute_on = 'timestep_end'
  [../]
[]
