
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 15
  ny = 15
  nz = 10
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

  potential_E_int = potential_int
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

  [./potential_int]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./bed_z]
    type = BulkEnergyDerivativeSixth
    variable = polar_z
    component = 2
  [../]

  ##Electrostatics
  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     permittivity = 1.0
  [../]

  [./polar_electric_pz]
     type = PolarElectricPStrong
     variable = polar_z
     component = 2
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
  [./top_electrode_top]
     type = DirichletBC
     variable = 'potential_int'
     value = 0.6
     boundary = 'front'
  [../]

  [./top_electrode_bottom]
     type = DirichletBC
     variable = 'potential_int'
     value = 0.0
     boundary = 'back'
  [../]
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
  print_perf_log = true
  [./out]
    type = Exodus
    execute_on = 'timestep_end'
    file_base = out_prob2a
    elemental_as_nodal = true
    interval = 1
  [../]
  [./outcsv]
    type = CSV
    file_base = out_prob2a
    execute_on = 'timestep_end'
  [../]
[]
