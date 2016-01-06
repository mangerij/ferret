[Mesh]
  file = exodus_cylinder_r3_h15.e
[]

[GlobalParams]

# Use the Srwolitz paper N. Ng et al. / Acta Materialia 60 (2012) 3632–3642 -PZT
  len_scale = 1.0
  alpha1 = -0.1484 # T = 300 K
  alpha11 = -0.0305
  alpha111 = 0.2475
  alpha12 = 0.632
  alpha112 = 0.96839
  alpha123 = -4.901
  G110 = 0.1483
  Q_mnkl = '0.08142 -0.02446 -0.02446 0.08142 -0.02446 0.08142 0.06417 0.06417 0.06417'
  permittivity = 0.5843763
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int #can add a potential_ext
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  displacements = 'disp_x disp_y disp_z'
  C_ijkl = '176.4 79.37 79.37 176.4 79.37 176.4 111. 111. 111.'
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
    [../]
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
    block = '1'
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
    block = '1'
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    block = '1'
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
     block = '1'
  [../]
  [./ferroelectriccouplingu_y]
     type = FerroelectricCouplingX
     variable = disp_y
     component = 1
     block = '1'
  [../]
  [./ferroelectriccouplingu_z]
     type = FerroelectricCouplingX
     variable = disp_z
     component = 2
     block = '1'
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
     block = '1'
  [../]
  [./FE_E_int]
     type=Electrostatics
     variable = potential_int
     block = '1'
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
  [./Periodic]
    [./polar_x_pbc]
      variable = polar_x #could leave empty and default is ALL variables in system...
      primary = '2'
      secondary = '3'
      translation = '0 0 15'
    [../]
    [./polar_y_pbc]
      variable = polar_y 
      primary = '2'
      secondary = '3'
      translation = '0 0 15'
    [../]
    [./polar_z_pbc]
      variable = polar_z 
      primary = '2'
      secondary = '3'
      translation = '0 0 15'
    [../]
    [./potential_int_pbc]
      variable = potential_int 
      primary = '2'
      secondary = '3'
      translation = '0 0 15'
    [../]
    [./disp_x_pbc]
      variable = disp_x
      primary = '2'
      secondary = '3'
      translation = '0 0 15'
    [../]
    [./disp_y_pbc]
      variable = disp_y
      primary = '2'
      secondary = '3'
      translation = '0 0 15'
    [../]
    [./disp_z_pbc]
      variable = disp_z
      primary = '2'
      secondary = '3'
      translation = '0 0 15'
    [../]
  [../]
[]

[Materials]
  [./slab_ferroelectric]
    type = LinearFerroelectricMaterial
    block = '1'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = '1'
    fill_method = symmetric9
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = '1'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '1'
  [../]
[]

[Postprocessors]
  [./bulk_energy]
   type = BulkEnergy
   execute_on = 'timestep_end'
   block = '1'
  [../]
  [./wall_energy]
   type = PZTWallEnergy
   execute_on = 'timestep_end'
   block = '1'
  [../]
  [./elastic_energy]
   type = ElasticEnergy
   execute_on = 'timestep_end'
   block = '1'
  [../]
  [./coupled_energy]
    type = CoupledEnergy
    execute_on = 'timestep_end'
    block = '1'
  [../]
  [./electrostatic_energy]
   type = ElectrostaticEnergy
   execute_on = 'timestep_end'
   block = '1'
  [../]
  [./total_energy_noelastic]
   type = TotalEnergyFlow
   bulk_energy = bulk_energy
   wall_energy = wall_energy
   bulk_energy_fourth = bulk_energy_fourth
   coupled_energy = coupled_energy
   electrostatic_energy = electrostatic_energy
   execute_on = 'timestep_end'
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
    petsc_options = '-snes_converged_reason -snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type    -sub_pc_factor_zeropivot -pc_factor_zeropivot'
    petsc_options_value = '    121            1e-8      1e-8    gamg        1e-50    1e-50'
  [../]
[]

[Executioner]
  type = Transient
  [./TimeStepper] #iterative DT halfs the time it takes to find a solution? oh well, our time is fake in this simulation anyway...
    type = IterationAdaptiveDT
    dt = 0.85 #max seems to be about 1.0 but could depend on refinement...
    #there is also a cutback on this for 0.2*optimal and yes i think it does count the 0th one.
    #iteration_window = 10
    optimal_iterations = 5 #i think this is 3 or more then cut? less than 3 grow, does it count the 0th iteration? no the cutting has to do with the iteration ratio
    growth_factor = 1.4
    linear_iteration_ratio = 100
    cutback_factor =  0.35
  [../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 1.285
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_PZT_cylinder_free_r3_h15
    interval = 3
  [../]
[]
