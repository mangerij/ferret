[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 5
  ny = 5
  nz = 2
  xmin = -2.5
  xmax = 2.5
  ymin = -2.5
  ymax = 2.5
  zmin = -1
  zmax = 1
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
  G110 = 0.15
  G11/G110 = 0.6
  G12/G110 = 0
  G44/G110 = 0.3
  G44P/G110 = 0.3
  Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.03375 0.03375 0.03375'
  permittivity = 0.5843763
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  displacements = 'disp_x disp_y disp_z'
  C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
[]



[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = ConstantIC
      value = 0.05
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = ConstantIC
      value = 0.05
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = ConstantIC
      value = 0.05
    [../]
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.0005
    [../]
  [../]

  [./disp_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.005
    [../]
    block = '0'
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.005
    [../]
    block = '0'
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = ConstantIC
      value = 0.005
    [../]
  [../]
[]

[AuxVariables]
  [./elastic_strain_yy]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_yy
  [../]
[]

[AuxKernels]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    block = '0'
    variable = elastic_strain_yy
    execute_on = 'timestep_end'
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

  [./disp_x_slab_f]
    type = DirichletBC
    variable = disp_x
    boundary = 'front'
    value = 0.0
  [../]
  [./disp_y_slab_f]
    type = DirichletBC
    variable = disp_y
    boundary = 'front'
    value = 0.0
  [../]
  [./disp_z_slab_f]
    type = DirichletBC
    variable = disp_z
    boundary = 'front'
    value = 0.0
  [../]

  [./disp_x_slab_b]
    type = DirichletBC
    variable = disp_x
    boundary = 'back'
    value = 0.0
  [../]
  [./disp_y_slab_b]
    type = DirichletBC
    variable = disp_y
    boundary = 'back'
    value = 0.0
  [../]
  [./disp_z_slab_b]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0.0
  [../]

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

[]

[Materials]
  [./slab_ferroelectric]
    type=LinearFerroelectricMaterial
    block = '0'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    #in GPA. from N. Pandech et al. Ceramic. Internat., times 1e9 to convert to N/m^2
    # C11 C12 C13 C22 C23 C33 C44 C55 C66
    #C_ijkl = '3.80e-7 1.5e-7 1.50e-7 3.80e-7 1.50e-7 3.80e-7 1.1e-7 1.1e-7 1.1e-7'
    #in m^4/C^2. from http://arxiv.org/pdf/1205.5640.pdf
    # Q11 Q12 Q13 Q22 Q23 Q33 Q44 Q55 Q66
    #euler_angle_1 = 0.0 #currently will only rotate C_ijkl
    #euler_angle_2 = 0.0
    #euler_angle_3 = 0.0
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
   block = '0'
  [../]
  [./wall_energy]
   type = WallEnergy
   execute_on = 'timestep_end'
   block = '0'
  [../]
  [./elastic_energy]
   type = ElasticEnergy
   execute_on = 'timestep_end'
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
    expression = 'perc_change <= 5.0e-4'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason '
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type    -pc_factor_zeropivot'
    petsc_options_value = '    121            1e-8      1e-8    gamg           1e-50     '
  [../]
[]

[Executioner]
  type = Transient
  [./TimeStepper] #iterative DT halfs the time it takes to find a solution? oh well, our time is fake in this simulation anyway...
    type = IterationAdaptiveDT
    dt = 0.35 #max seems to be about 1.0 but could depend on refinement...
    #there is also a cutback on this for 0.2*optimal and yes i think it does count the 0th one.
    #iteration_window = 10
    optimal_iterations = 4 #i think this is 3 or more then cut? less than 3 grow, does it count the 0th iteration? no the cutting has to do with the iteration ratio
    growth_factor = 1.4
    linear_iteration_ratio = 100
    #linear_iteration_ratio = 1000
    cutback_factor =  0.55
    num_steps = 100
  [../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 0.25
  num_steps = 100
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_coupled_test
    output_initial = true
    elemental_as_nodal = true
    interval = 1
  [../]
[]
