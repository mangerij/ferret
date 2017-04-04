
[Mesh]
  file = exodus_disk_r8_h1.e
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

  #epsilon = -0.001 #negative = tension, positive = compression. This tunes skyrm->vortex, skyrm->cd as in analytical work
  epsilon = def

  #T = 298
  T = def
  

  #lambda = 0.002
  lambda = def

[]

[Functions]
  #-----------------------------------------------------#
  # This is a skyrmion solution from Svitlana and Igor  #
  # It was fit in Mathematica                           #
  #-----------------------------------------------------#


  [./parsed_function_x_skyrm]
    type = ParsedFunction
    value = '-(0.738217-0.00686984*(x^2+y^2)^(0.5)+0.00644497*(x^2+y^2)-0.0188174*(x^2+y^2)^(1.5)+0.00441745*(x^2+y^2)^2-0.000274842*(x^2+y^2)^(5/2))*sin(-0.028395+0.267482*(x^2+y^2)^(0.5)-0.146762*(x^2+y^2)+0.0632932*(x^2+y^2)^(1.5)-0.00790942*(x^2+y^2)^(2)+0.000294936*(x^2+y^2)^(5/2))*sin(atan(y/x))'
  [../]
  [./parsed_function_y_skyrm]
    type = ParsedFunction
    value = '(0.738217-0.00686984*(x^2+y^2)^(0.5)+0.00644497*(x^2+y^2)-0.0188174*(x^2+y^2)^(1.5)+0.00441745*(x^2+y^2)^2-0.000274842*(x^2+y^2)^(5/2))*sin(-0.028395+0.267482*(x^2+y^2)^(0.5)-0.146762*(x^2+y^2)+0.0632932*(x^2+y^2)^(1.5)-0.00790942*(x^2+y^2)^(2)+0.000294936*(x^2+y^2)^(5/2))*cos(atan(y/x))'
  [../]
  [./parsed_function_z_skyrm]
    type = ParsedFunction
    value = '(0.738217-0.00686984*(x^2+y^2)^(0.5)+0.00644497*(x^2+y^2)-0.0188174*(x^2+y^2)^(1.5)+0.00441745*(x^2+y^2)^2-0.000274842*(x^2+y^2)^(5/2))*cos(-0.028395+0.267482*(x^2+y^2)^(0.5)-0.146762*(x^2+y^2)+0.0632932*(x^2+y^2)^(1.5)-0.00790942*(x^2+y^2)^(2)+0.000294936*(x^2+y^2)^(5/2))'
  [../]

  #-----------------------------------------------------#
  # This is a cylindrical domain solution from Svitlana #
  # and Igor. It was fit in Mathematica.                #
  #-----------------------------------------------------#


  [./parsed_function_z_cd]
    type = ParsedFunction
    value = '0.719527-0.0061793*(x^2+y^2)^(0.5)-0.00641062*(x^2+y^2)+0.00508983*(x^2+y^2)^(1.5)-0.0020986*(x^2+y^2)^2+0.000171088*(x^2+y^2)^(5/2)'
  [../]
[]

[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = parsed_function_x_skyrm
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = parsed_function_y_skyrm
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = parsed_function_z_skyrm
    [../]
  [../]
  [./potential_int]
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
[]


[AuxVariables]
  [./chern]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]


[AuxKernels]
  [./cherndens]
    type = ChernSimonsDensity
    variable = chern
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
     block = '1'
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
  boundary = '2'
  value = -0.001
[../]

[./potential_int_2]
  type = DirichletBC
  variable = potential_int
  boundary = '3'
  value = -0.001
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
    petsc_options_value = '    121               1e-10     1e-8      1e-6    bjacobi'
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
    file_base = out_skyrm
    elemental_as_nodal = true
    interval = 1
    execute_on = 'timestep_end'
  [../]
  [./outcsv]
    type = CSV
    file_base = out_skyrm
    execute_on = 'timestep_end'
  [../]
[]
