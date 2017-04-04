
[Mesh]
  file = exodus_thinfilm_test_085_20_20_10.e
[]

[GlobalParams]
  len_scale = 1.0
  alpha1 = -0.1722883 # (3.766(T-765.1)*10^5) C^{-2} nm^2 (T = 293 K)
  alpha11 = -0.07253
  alpha111 = 0.26
  alpha12 = 0.75
  alpha112 = 0.61
  alpha123 = -3.67
  G110 = 0.173
  G11/G110 = 0.6
  G12/G110 = 0
  G44/G110 = 0.3
  G44P/G110 = 0.3
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int

  epsilon = 0.01 #negative = tension, positive = compression
[]



[Variables]
  [./polar_x]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./polar_y]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./polar_z]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
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
    [../]
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
     permittivity = 0.08854187
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     block = '1'
     permittivity = 0.08854187
  [../]

  [./DIE_E_int]
     type = Electrostatics
     variable = potential_int
     block  = '2'
     permittivity = 2.6562561
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
  [./bot_potential_int]
    variable = potential_int
    type = DirichletBC
    value = 0
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
   [./Fbulk]
      type = BulkEnergy
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
      permittivity = 0.08854187
      execute_on = 'timestep_end'
    [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart  -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121                1e-10      1e-8      1e-8    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
    [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.8
    #iteration_window = 3
    optimal_iterations = 6 #should be 5 probably
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.8
[../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  #dt = 0.5
  dtmin = 1e-13
  dtmax = 0.8
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = outPTO_thinfilm_09_2_2_1_t01_STO
    elemental_as_nodal = true
    interval = 1
  [../]
[]

