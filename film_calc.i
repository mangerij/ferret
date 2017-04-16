
[Mesh]
  file = exodus_thinfilm_test_085_20_20_10.e
[]

[GlobalParams]
  len_scale = 1.0
  G110 = 0.173
  G11/G110 = 2.0
  G12/G110 = 0
  G44/G110 = 1.0
  G44P/G110 = 1.0
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int

  #epsilon = -0.001 #negative = tension, positive = compression.
  
  epsilon = 0.01

  #T = 298

  T = 298
  
  lambda = 0.0
[]



[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
  [../]
[]


[AuxVariables]
  [./chern]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./chernMag]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]


[AuxKernels]
  [./cherndens]
    type = ChernSimonsDensity  #Need to code in skyrmion number instead of this
    variable = chern
  [../]
  [./chernMagdens]
    type = ChernSimonsDensityMag  #Need to code in skyrmion number instead of this
    variable = chernMag
  [../]
[]

[Kernels]
  #Bulk energy density
  [./bed_x]
    type = RenormalizedFreeEnergy
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = RenormalizedFreeEnergy
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = RenormalizedFreeEnergy
    variable = polar_z
    component = 2
  [../]
  ##Wall energy penalty
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

  ##Electrostatics
  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     block = '1 2'
     permittivity = 0.08854187
  [../]
  [./polar_electric_px]
     type = PolarElectricPStrong
     variable = polar_x
     component = 0
  [../]
  [./polar_electric_py]
     type = PolarElectricPStrong
     variable = polar_y
     component = 1
  [../]
  [./polar_electric_pz]
     type = PolarElectricPStrong
     variable = polar_z
     component = 2
  [../]


  [./depol_z]
    type = DepolEnergy
    permitivitty = 0.00885
    variable = polar_z
    avePz = avePz
  [../]

  ##Time dependence
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

  [./potential_int_1]
    type = DirichletBC
    variable = potential_int
    boundary = '1'
    value = -0.001
  [../]

  [./bot_potential_int]
    variable = potential_int
    type = DirichletBC
    value = -0.001
    boundary = '7'
  [../]

  [./Periodic]
    [./TB_polar_x_pbc]
      variable = polar_x
      primary = '3'
      secondary = '5'
      translation = '0 20 0'
    [../]
    [./TB_polar_y_pbc]
      variable = polar_y
      primary = '3'
      secondary = '5'
      translation = '0 20 0'
    [../]
    [./TB_polar_z_pbc]
      variable = polar_z
      primary = '3'
      secondary = '5'
      translation = '0 20 0'
    [../]
    [./TB_potential_int_pbc]
      variable = potential_int
      primary = '3'
      secondary = '5'
      translation = '0 20 0'
    [../]
  #

    [./RL_polar_x_pbc]
      variable = polar_x
      primary = '4'
      secondary = '6'
      translation = '20 0 0'
    [../]
    [./RL_polar_y_pbc]
      variable = polar_y
      primary = '4'
      secondary = '6'
      translation = '20 0 0'
    [../]
    [./RL_polar_z_pbc]
      variable = polar_z
      primary = '4'
      secondary = '6'
      translation = '20 0 0'
    [../]
    [./RL_potential_int_pbc]
      variable = potential_int
      primary = '4'
      secondary = '6'
      translation = '20 0 0'
    [../]
  [../]
[]



[Postprocessors]
   [./avePz]
     type = ElementAverageValue
     variable = polar_z
     execute_on = 'initial linear nonlinear timestep_begin timestep_end'
   [../]
   [./avgChern]
     block = '1'
     type = ElementAverageValue
    variable = chern
   [../]
   [./avgChernMag]
     block = '1'
     type = ElementAverageValue
    variable = chernMag
   [../]
   [./Fbulk]
      type = RenormalizedBulkEnergy
      block = '1'
      execute_on = 'timestep_end'
    [../]
    [./Fwall]
      type = WallEnergy
      block = '1'
      execute_on = 'timestep_end'
    [../]
    [./Felec]
      block = '1'
      type = ElectrostaticEnergy
      execute_on = 'timestep_end'
    [../]
    [./Fdepol]
      block = '1'
      type = DepolarizationEnergy
      execute_on = 'timestep_end'
    [../]
    [./Ftotal]
      type = TotalEnergyP
      Fbulk = Fbulk
      Fwall = Fwall
      Fec = Fec
      Fdepol = Fdepol
      execute_on = 'timestep_end'
    [../]
    [./perc_change]
     type = PercentChangePostprocessor
     postprocessor = Ftotal
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
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121               1e-10     1e-8      1e-8    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
    [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    #iteration_window = 3
    optimal_iterations = 6 #should be 5 probably
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.8
[../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 0.7
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_film_p01
    elemental_as_nodal = true
    interval = 1
    execute_on = 'timestep_end'
  [../]
  [./outcsv]
    type = CSV
    file_base = out_film_p01
    execute_on = 'timestep_end'
  [../]
[]
