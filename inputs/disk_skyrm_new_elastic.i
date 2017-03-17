
[Mesh]
  file = exodus_disk_r8_h5.e
[]


[GlobalParams]
  len_scale = 1.0

  alpha1 = -0.1722883 #room temp PTO
  alpha11 = -0.07253
  alpha111 = 0.26
  alpha12 = 0.75
  alpha112 = 0.61
  alpha123 = -3.67

  G110 = 0.173
  G11/G110 = 1.6
  G12/G110 = 0
  G44/G110 = 0.8
  G44P/G110 = 0.8

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  potential_int = potential_int

  # Elastic interactions are turned off for now.
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z

  displacements = 'disp_x disp_y disp_z'
  prefactor = 0.00 #negative = tension, positive = compression
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
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]


[Materials]
  [./eigen_strain_zz] #Use for stress-free strain (ie epitaxial)
    type = ComputeEigenstrain
    block = '1'
    # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
    eigen_base = '1 0 0 0 1 0 0 0 0'
    eigenstrain_name = eigenstrain
  [../]
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #from MaterialsProject
    C_ijkl = '281 115.74 115.74 281 115.74 281 97.18 97.18 97.18'
    block = '1'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    block = '1'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
    block = '1'
  [../]

  [./slab_ferroelectric]
    block = '1'
    type = ComputeElectrostrictiveTensor
    Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.03375 0.03375 0.03375'
    C_ijkl = '281 115.74 115.74 281 115.74 281 97.18 97.18 97.18'
  [../]

  #[./elasticity_tensor_2]
  #  type = ComputeElasticityTensor
  #  #averaged from BulkMod
  #  C_ijkl = '281 115.74 115.74 281 115.74 281 97.18 97.18 97.18'
  #  fill_method = symmetric9
  #  block = '2'
  #[../]
  #[./strain_2]
  #  type = ComputeSmallStrain
  #  block = '2'
  #  eigenstrain_names = eigenstrain
  #[../]
  #[./stress_2]
  #  type = ComputeLinearElasticStress
  #  block = '2'
  #[../]
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
  ##Polarization-strain coupling

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


    [./ferroelectriccouplingX_xx]
      type = FerroelectricCouplingX
      block = '1'
      variable = disp_x
      component = 0
    [../]
    [./ferroelectriccouplingX_yy]
      type = FerroelectricCouplingX
      block = '1'
      variable = disp_y
      component = 1
    [../]
    [./ferroelectriccouplingX_zz]
      type = FerroelectricCouplingX
      block = '1'
      variable = disp_z
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

    #[./METAL_E_int]
    #   type = Electrostatics
    #   variable = potential_int
    #   block  = '2'
    #   permittivity = 18.0
    #[../]
    #[./METAL_TF]
    #   type = ThomasFermiPotential
    #   variable = potential_int
    #   block  = '2'
    #   TFconstant = -0.00001
    #[../]

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

   [./top_electrode_top_bottom]
      type = DirichletBC
      variable = 'potential_int'
      value = 0.001
     boundary = '2 3'
   [../]

   [./side_neumann_x]
    variable = polar_x
    type = NeumannBC
    value = 0
    boundary = '1'
  [../]
  [./side_neumann_y]
    variable = polar_y
    type = NeumannBC
    value = 0
    boundary = '1'
  [../]
  [./side_neumann_z]
    variable = polar_z
    type = NeumannBC
    value = 0
    boundary = '1'
  [../]
[]



[Postprocessors]
   [./avePz]
     type = ElementAverageValue
     block = '1'
     variable = polar_z
     execute_on = 'initial linear nonlinear timestep_begin timestep_end'
   [../]

   [./Felastic]
     type = ElasticEnergy
     block = '1'
     execute_on = 'timestep_end'
   [../]

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
    [./Fcoupled]
      type = CoupledEnergy
      block = '1'
      execute_on = 'timestep_end'
    [../]
    [./Felec]
      block = '1'
      type = ElectrostaticEnergy
      permittivity = 0.08854187
      execute_on = 'timestep_end'
    [../]
    [./Ftotal]
      type = TotalEnergyFlow
      #type = TotalEnergyFlowNoElastNoElec
      Fbulk = Fbulk
      Fwall = Fwall
      Fcoupled = Fcoupled
      Felec = Felec
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
    petsc_options_value = '    250               1e-10      1e-6      1e-4    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
    [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.2
    #iteration_window = 3
    optimal_iterations = 6 #should be 5 probably
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.65
[../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  #dt = 0.5
  dtmin = 1e-13
  dtmax = 1.5
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    execute_on = 'timestep_end'
    file_base = out_nanodisk_elastic_SHORT2
    elemental_as_nodal = true
    interval = 1
  [../]
[]
