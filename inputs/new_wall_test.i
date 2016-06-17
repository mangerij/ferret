[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 28
  ny = 28
  nz = 16
  xmin = -12
  xmax = 12
  ymin = -12
  ymax = 12
  zmin = -6
  zmax = 6
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
  G110 = 0.173
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
    block = '0'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-6
      max = 0.5e-6
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-6
      max = 0.5e-6
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-6
      max = 0.5e-6
    [../]
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-6
      max = 0.5e-6
    [../]
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
     variable = potential_int
     block = '0'
  [../]
  [./FE_E_int]
     type=Electrostatics
     variable = potential_int
     block = '0'
     permittivity = 0.08854187
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
    variable = potential_int
  [../]
  [./potential_cube6]
    type = DirichletBC
    boundary = 'back'
    value = 0.0002
    variable = potential_int
  [../]
 [./Periodic] #PBC ALONG Y
    [./TB_polar_x_pbc]
      variable = polar_x
      primary = 'bottom'
      secondary = 'top'
      translation = '0 24 0'
    [../]
    [./TB_polar_y_pbc]
      variable = polar_y
      primary = 'bottom'
      secondary = 'top'
      translation = '0 24 0'
    [../]
    [./TB_polar_z_pbc]
      variable = polar_z
      primary = 'bottom'
      secondary = 'top'
      translation = '0 24 0'
    [../]
    [./TB_potential_int_pbc]
      variable = potential_int
      primary = 'bottom'
      secondary = 'top'
      translation = '0 24 0'
    [../]


  [../]
[]


[Postprocessors]
  [./Fbulk]
   type = BulkEnergy
   execute_on = 'timestep_end'
   block = '0'
  [../]
  [./Fwall]
   type = WallEnergy
   execute_on = 'timestep_end'
   block = '0'
  [../]

  [./Ftotal]
    type = TotalEnergyFlow
    Fbulk = Fbulk
    Fwall = Fwall
    Fcoupled = Fcoupled
    Felec = Felec
     execute_on = 'timestep_end'
  [../]
  [./electrostatic_energy]
   type = ElectrostaticEnergy
   execute_on = 'timestep_end'
   block = '0'
  [../]
    [./perc_change]
     type = PercentChangePostprocessor
     postprocessor = Ftotal
   [../]
[]

[UserObjects]
 [./kill]
  type = Terminator
  expression = 'perc_change <= 7.5e-3'
 [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason '
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type    -pc_factor_zeropivot'
    petsc_options_value = '    121            1e-6      1e-8    bjacobi           1e-50     '
  [../]
[]

[Executioner]
  type = Transient
  [./TimeStepper] 
    type = IterationAdaptiveDT
    dt = 1.6
    optimal_iterations = 5 
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.65
  [../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-8
  dtmax = 1.6
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_domains
    elemental_as_nodal = true
  [../]
[]
