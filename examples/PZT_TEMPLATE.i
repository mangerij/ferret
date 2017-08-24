
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 10
  ny = 10
  nz = 10
  xmin = -5
  xmax = 5
  ymin = -5
  ymax = 5
  zmin = -5
  zmax = 5
  elem_type = HEX8
[]


[GlobalParams]

# Use the Srwolitz paper N. Ng et al. / Acta Materialia 60 (2012) 3632–3642
# for PZT.

  len_scale = 1.0
  alpha1 = -0.1484 # (3.766(T-765.1)*10^5) C^{-2} nm^2 (T = 300 K)
  alpha11 = -0.0305
  alpha111 = 0.2475
  alpha12 = 0.632
  alpha112 = 0.96839
  alpha123 = -4.901
  G110 = 0.28
  G11_G110 = 1.0
  G12/G110 = 0
  G44/G110 = 0.5
  G44P/G110 = 0.5
  Q_mnkl = '0.08142 -0.02446 -0.02446 0.08142 -0.02446 0.08142 0.06417 0.06417 0.06417'
  permittivity = 0.5843763
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
  #potential_ext = potential_ext
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  prefactor = 0.002
  displacements = 'disp_x disp_y disp_z'
  #use_displaced_mesh = false
  C_ijkl = '176.4 79.37 79.37 176.4 79.37 176.4 111. 111. 111.'
  #initial_from_file_timestep = 120
[]



[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block = '0'
   # initial_from_file_var = polar_x
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block = '0'
  # initial_from_file_var = polar_y
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
    #initial_from_file_var = polar_y
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
    #initial_from_file_var = polar_z
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
    #initial_from_file_var = potential_int
  [../]
  #[./potential_ext]
  #  order=FIRST
  #  family = LAGRANGE
  #  #initial_from_file_var = potential_ext
  #[../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
    block = '0'
    #scaling = 1e6
    #initial_from_file_var = disp_x
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
    block = '0'
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
[]

[Kernels]
  #Elastic problem
  [./TensorMechanics]
     #This is an action block
  [../]
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
     type = PZTWallEnergyDerivative
     variable = polar_x
     component = 0
  [../]
  [./walled_y]
     type = PZTWallEnergyDerivative
     variable = polar_y
     component = 1
  [../]
  [./walled_z]
     type = PZTWallEnergyDerivative
     variable = polar_z
     component = 2
  [../]
  ##Polarization-strain coupling
  [./ferroelectriccouplingu_x]
     type = FerroelectricCouplingX
     variable = disp_x
     component = 0
     block = '0'
  [../]
  [./ferroelectriccouplingu_y]
     type = FerroelectricCouplingX
     variable = disp_y
     component = 1
     block = '0'
  [../]
  [./ferroelectriccouplingu_z]
     type = FerroelectricCouplingX
     variable = disp_z
     component = 2
     block = '0'
  [../]
  [./ferroelectriccouplingp_xx]
     type = FerroelectricCouplingP
     variable = polar_x
     component = 0
  [../]
  [./ferroelectriccouplingp_yy]
     type = FerroelectricCouplingP
     variable = polar_y
     component = 1
  [../]
  [./ferroelectriccouplingp_zz]
     type = FerroelectricCouplingP
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
     type=TimeDerivative
     variable=polar_x
  [../]
  [./polar_y_time]
     type=TimeDerivative
     variable=polar_y
  [../]
  [./polar_z_time]
     type=TimeDerivative
     variable = polar_z
  [../]
[]

[BCs]

 #[./Periodic]
 #  #active='disp_x_x disp_y_x polar_z_x disp_x_y disp_y_y disp_z_y'
 #  #active='polar_x_x polar_y_x polar_z_x polar_x_y polar_y_y polar_z_y'
 #  #active='potential_int_x potential_ext_x potential_int_y potential_ext_y'
 #
 #  [./potential_int_x]
 #     variable = potential_int
 #     auto_direction = 'x y'
 #  [../]
 #
 #  [./polar_x]
 #     variable = polar_x
 #     auto_direction = 'x y'
 #  [../]
 #  [./polar_y]
 #     variable = polar_y
 #     auto_direction = 'x y'
 #  [../]
 #  [./polar_z]
 #     variable = polar_z
 #     auto_direction = 'x y'
 #  [../]
 #
 #  [./disp_x]
 #     variable = disp_x
 #     auto_direction = 'x y'
 #  [../]
 #  [./disp_y]
 #     variable = disp_y
 #     auto_direction = 'x y'
 #
 #  [../]
 #  [./disp_z]
 #     variable = disp_z
 #     auto_direction = 'x y'
 #  [../]
 # [../]

[]

[Materials]
  [./slab_ferroelectric]
    type=LinearFerroelectricMaterial
    block = '0'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = '0'
    fill_method = symmetric9
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = '0'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '0'
  [../]

[]

[Postprocessors]
  [./bulk_energy]
   type = BulkEnergy
   execute_on = 'timestep_end'
  # initial_from_file_var = bulk_energy
   block = '0'
  [../]
  [./wall_energy]
   type = PZTWallEnergy
   execute_on = 'timestep_end'
  # initial_from_file_var = wall_energy
   block = '0'
  [../]
  [./elastic_energy]
   type = ElasticEnergy
   execute_on = 'timestep_end'
  # initial_from_file_var = elastic_energy
   block = '0'
  [../]
  [./coupled_energy]
    type = CoupledEnergy
    execute_on = 'timestep_end'
    block = '0'
  [../]
  [./electrostatic_energy]
   type = ElectrostaticEnergy
   execute_on = 'timestep_end'
   block = '0'
  # initial_from_file_var = electrostatic_energy
  [../]
  [./total_energy_noelastic]
   type = TotalEnergyFlow
   bulk_energy = bulk_energy
   wall_energy = wall_energy
   bulk_energy_fourth = bulk_energy_fourth
   coupled_energy = coupled_energy
   electrostatic_energy = electrostatic_energy
   execute_on = 'timestep_end'
  # initial_from_file_var = total_energy
  [../]
  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = total_energy_noelastic
  [../]
[]

[UserObjects]
  [./kill]
    type = Terminator
    expression = 'perc_change <= 5.0e-5'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -options_left -snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type    -sub_pc_factor_zeropivot -pc_factor_zeropivot -pc_side '
    petsc_options_value = '    121            1e-8      1e-8    gamg        1e-50    1e-50      left        '
  [../]
[]

#Limits exist on -snes_rtol =< 1e-10.

[Executioner]
  type = Transient
  [./TimeStepper] #iterative DT halfs the time it takes to find a solution? oh well, our time is fake in this simulation anyway...
    type = IterationAdaptiveDT
    dt = 0.25 #max seems to be about 1.0 but could depend on refinement...
    #there is also a cutback on this for 0.2*optimal and yes i think it does count the 0th one.
    #iteration_window = 10
    optimal_iterations = 4 #i think this is 3 or more then cut? less than 3 grow, does it count the 0th iteration? no the cutting has to do with the iteration ratio
    growth_factor = 1.4
    linear_iteration_ratio = 100
    #linear_iteration_ratio = 1000
    cutback_factor =  0.55
  [../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  #dt = 2.0
  dtmin = 1e-13
  dtmax = 0.25
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_PbZrTiO3_test
    output_initial = true
    elemental_as_nodal = true
    interval = 2
  [../]
[]
