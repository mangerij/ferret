

# This file just evolves from the PE state an open circuit BC
# chunk of LNO using coefficients from Chen's appendix 
# (and assuming gradient terms are == PTO).

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 25
  ny = 25
  nz = 20
  xmin = -8
  xmax = 8
  ymin = -8
  ymax = 8
  zmin = -7
  zmax = 7
  elem_type = HEX8
[]

[GlobalParams]
  len_scale = 1.0
  alpha1 = 2.012
  alpha2 = 3.608
  alpha3 = 1.345

  G110 = 0.253 #from PRB, 71, 184110 (2005)
  G11/G110 = 0.6
  G12/G110 = 0
  G44/G110 = 0.3
  G44P/G110 = 0.3

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int

[]



[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-5
      max = 0.1e-5
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-5
      max = 0.1e-5
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-5
      max = 0.1e-5
    [../]
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
  [../]
[]


[Kernels]
  [./bed_x]
    type = LBOBulkEnergyDeriv
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = LBOBulkEnergyDeriv
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = LBOBulkEnergyDeriv
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
  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
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

[Functions]
  [./bc_func_1]
    type = ParsedFunction
    value = '10*sin(omega*t)'
    vars = 'omega'
    vals = '0.1'
  [../]
  [./bc_func_2]
    type = ParsedFunction
    value = '-10*sin(omega*t)'
    vars = 'omega'
    vals = '0.1'
  [../]
[]

[BCs]
  # Boundary Condition System
  [./front_pot]
    type = FunctionDirichletBC
    variable = potential_int
    boundary = 'front'
    function = bc_func_1
  [../]
  [./back_pot]
    type = FunctionDirichletBC
    variable = potential_int
    boundary = 'back'
    function = bc_func_2
  [../]
[]




[Postprocessors]
    [./avePz]
      type = ElementAverageValue
      variable = polar_z
      execute_on = 'initial linear nonlinear timestep_begin timestep_end'
    [../]
    [./avePz]
      type = ElementAverageValue
      variable = polar_z
      execute_on = 'timestep_end'
    [../]
    [./Fbulk]
      type = LBOBulkEnergy
      execute_on = 'timestep_end'
    [../]
    [./Fwall]
      type = WallEnergy
      execute_on = 'timestep_end'
    [../]
    [./Felec]
      type = ElectrostaticEnergy
      execute_on = 'timestep_end'
    [../]
    [./Ftotal]
      type = TotalEnergyFlowNoElast
      Fbulk = Fbulk
      Fwall = Fwall
      Felec = Felec
      execute_on = 'timestep_end'
    [../]
    [./perc_change]
     type = PercentChangePostprocessor
     postprocessor = Ftotal
   [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart  -snes_atol -ksp_rtol -pc_type'
    petsc_options_value = '    121                1e-10      1e-8    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
    [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.08
    optimal_iterations = 6
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.8
[../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  #dt = 0.5
  dtmin = 1e-13
  dtmax = 0.08
  num_steps = 1000
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = outLBO_time_test_E1
    elemental_as_nodal = true
    interval = 1
  [../]
  [./outCSV]
    type = CSV
    file_base = outLBO_time_test_E1
  [../]
[]
