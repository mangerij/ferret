[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 6
  ny = 6
  nz = 50
  xmin = -3
  xmax = 3
  ymin = -3
  ymax = 3
  zmin = -12
  zmax = 12
  elem_type = HEX8
[]

[GlobalParams]
  #len_scale = 1.0
  #alpha0 = 0.0003766
  #alpha11 = -0.07253
  #alpha111 = 0.26
  #alpha12 = 0.75
  #alpha112 = 0.61
  #alpha123 = -3.67
  #G110 = 0.173
  #G11/G110 = 2.0
  #G12/G110 = 0
  #G44/G110 = 1.0
  #G44P/G110 = 1.0
  #polar_x = polar_x
  #polar_y = polar_y
  #polar_z = polar_z
  #potential_int = potential_int
  temperature = temperature
  #Tc = 765.1
[]

[Variables]
  #[./polar_x]
  #  order = FIRST
  #  family = LAGRANGE
  #  block = '0'
  #  [./InitialCondition]
  #    type = RandomIC
  #    min = -0.5e-6
  #    max = 0.5e-6
  #  [../]
  #[../]
  #[./polar_y]
  #  order = FIRST
  #  family = LAGRANGE
  #  block = '0'
  #  [./InitialCondition]
  #    type = RandomIC
  #    min = -0.5e-6
  #    max = 0.5e-6
  #  [../]
  #[../]
  #[./polar_z]
  #  order = FIRST
  #  family = LAGRANGE
  #  block = '0'
  #  [./InitialCondition]
  #    type = RandomIC
  #    min = -0.5e-6
  #    max = 0.5e-6
  #  [../]
  #[../]
  #[./potential_int]
  #  order = FIRST
  #  family = LAGRANGE
  #  [./InitialCondition]
  #    type = RandomIC
  #    min = -0.5e-6
  #    max = 0.5e-6
  #  [../]
  #[../]
  [./temperature]
     [./InitialCondition]
      type = RandomIC
      min = 1
      max = 300
    [../]
  [../]
[]

[Kernels]
  ##Bulk energy density
  #[./bed_x]
  #  type = BulkEnergyDerivativeSixthCoupledT
  #  variable = polar_x
  #  component = 0
  #[../]
  #[./bed_y]
  #  type = BulkEnergyDerivativeSixthCoupledT
  #  variable = polar_y
  #  component = 1
  #[../]
  #[./bed_z]
  #  type = BulkEnergyDerivativeSixthCoupledT
  #  variable = polar_z
  #  component = 2
  #[../]
  ###Wall energy penalty
  #[./walled_x]
  #   type=WallEnergyDerivative
  #   variable = polar_x
  #   component = 0
  #[../]
  #[./walled_y]
  #   type=WallEnergyDerivative
  #   variable = polar_y
  #   component = 1
  #[../]
  #[./walled_z]
  #   type=WallEnergyDerivative
  #   variable = polar_z
  #   component = 2
  #[../]
  #
  ###Electrostatics
  #[./polar_x_electric_E]
  #   type=PolarElectricEStrong
  #   variable = potential_int
  #   block = '0'
  #[../]
  #[./FE_E_int]
  #   type=Electrostatics
  #   variable = potential_int
  #   block = '0'
  #   permittivity = 0.08854187
  #[../]
  #
  #[./polar_electric_px]
  #   type=PolarElectricPStrong
  #   variable = polar_x
  #   component = 0
  #[../]
  #[./polar_electric_py]
  #   type=PolarElectricPStrong
  #   variable = polar_y
  #   component = 1
  #[../]
  #[./polar_electric_pz]
  #   type=PolarElectricPStrong
  #   variable = polar_z
  #   component = 2
  #[../]

  ##Thermal operators

#  [./k_op_temp]
#     type = KarmanenkoDriver
#     variable = temperature
#     C1 = 1e3
#     C2 = 1e5
#     dEstep = 20.5
#  [../]

  [./T_diff]
     type = KappaTDiffusion
     c0 = -5.77e10
     c1 = 3.04e9
     c2 = 2.85e8
     c3 = 1.59e7
     c4 = 9.66e6
     #c5 = 1.02
     variable = temperature
     block = '0'
  [../]


  ##Time dependence
  #[./polar_x_time]
  #   type=TimeDerivativeScaled
  #   variable=polar_x
  #  time_scale = 1.0
  #[../]
  #[./polar_y_time]
  #   type=TimeDerivativeScaled
  #   variable=polar_y
  #  time_scale = 1.0
  #[../]
  #[./polar_z_time]
  #   type=TimeDerivativeScaled
  #   variable = polar_z
  #  time_scale = 1.0
  #[../]
  #[./temperature_time]
#    type=TimeDerivativeScaled
#    variable = temperature
#    time_scale = 1.0
#  [../]
[]

[BCs]
  #[./potential_cube5]
  #  type = DirichletBC
  #  boundary = 'front'
  #  value = 0.0002
  #  variable = potential_int
  #[../]
  #[./potential_cube6]
  #  type = DirichletBC
  #  boundary = 'back'
  #  value = 0.0002
  #  variable = potential_int
  #[../]

  [./Tcube5]
    type = PresetBC
    boundary = 'front'
    value = 1
    variable = temperature
  [../]
  [./Tcube6]
    type = PresetBC
    boundary = 'back'
    value = 300
    variable = temperature
  [../]

 #[./Periodic] #PBC ALONG Y
 #   [./TB_polar_x_pbc]
 #     variable = polar_x
 #     primary = 'bottom'
 #     secondary = 'top'
 #     translation = '0 24 0'
 #   [../]
 #   [./TB_polar_y_pbc]
 #     variable = polar_y
 #     primary = 'bottom'
 #     secondary = 'top'
 #     translation = '0 24 0'
 #   [../]
 #   [./TB_polar_z_pbc]
 #     variable = polar_z
 #     primary = 'bottom'
 #     secondary = 'top'
 #     translation = '0 24 0'
 #   [../]
 #   [./TB_potential_int_pbc]
 #     variable = potential_int
 #     primary = 'bottom'
 #     secondary = 'top'
 #     translation = '0 24 0'
 #   [../]
 # [../]
[]


#[Postprocessors]
  #[./Fbulk]
  # type = BulkEnergyCoupledT
  # execute_on = 'timestep_end'
  # block = '0'
  #[../]
  #[./Fwall]
  # type = WallEnergy
  # execute_on = 'timestep_end'
  # block = '0'
  #[../]

  #[./Ftotal]
  #  type = TotalEnergyFlow
  #  Fbulk = Fbulk
  #  Fwall = Fwall
  #  Fcoupled = Fcoupled
  #  Felec = Felec
  #   execute_on = 'timestep_end'
  #[../]
  #[./electrostatic_energy]
  # type = ElectrostaticEnergy
  # execute_on = 'timestep_end'
  # block = '0'
  #[../]
  #  [./perc_change]
  #   type = PercentChangePostprocessor
  #   postprocessor = Ftotal
  # [../]
#[]

#[UserObjects]
# [./kill]
#  type = Terminator
#  expression = 'perc_change <= 7.5e-3'
# [../]
#[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_shift_amount -pc_factor_shift_type -pc_factor_shift_amount'
    petsc_options_value = ' bjacobi    lu              NONZERO                  1e-10 NONZERO               1e-10'
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  #type = Transient
  #[./TimeStepper]
  #  type = IterationAdaptiveDT
  #  dt = 0.008
  #  optimal_iterations = 5
  #  growth_factor = 1.4
  #  linear_iteration_ratio = 1000
  #  cutback_factor =  0.65
  #[../]
  #solve_type = 'PJFNK'       #"PJFNK, JFNK, NEWTON"
  #scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  #dtmin = 1e-12
  #dtmax = 0.008
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_domains_tempKDrive
    elemental_as_nodal = true
  [../]
[]