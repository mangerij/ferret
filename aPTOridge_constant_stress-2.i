
[Mesh]
  file = exodus_antiparallel.e
  #block = '1 2' PTO
  #block = '4'   substrate
  #block = '3'   vacuum
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
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./polar_x]
    block = '1 2'
    order = FIRST
    family = LAGRANGE
  [../]
  [./polar_y]
    block = '1 2'
    order = FIRST
    family = LAGRANGE
  [../]
  [./polar_z]
    block = '1 2'
    order = FIRST
    family = LAGRANGE
  [../]
  [./potential_int]
    block = '1 2 3 4'
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_x]
    block = '1 2 4'
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    block = '1 2 4'
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    block = '1 2 4'
    order = FIRST
    family = LAGRANGE
  [../]
[]


[ICs]
  [./polarzblock1]
    type = RandomIC
    variable = polar_z
    block = '1'
    min = -0.6e-3 #-p
    max = 0.6e-3
  [../]
  [./polarzblock2]
    variable = polar_z
    type = RandomIC
    block = '2'
    min = -0.6e-3 #+P
    max = 0.6e-3
  [../]

  [./polarxblock1]
    type = RandomIC
    variable = polar_x
    block = '1'
    min = -0.6e-3
    max = 0.6e-3
  [../]
  [./polarxblock2]
    variable = polar_x
    type = RandomIC
    block = '2'
    min = -0.6e-3
    max = 0.6e-3
  [../]

  [./polaryblock1]
    type = RandomIC
    variable = polar_y
    block = '1'
    min = -0.6e-3
    max = 0.6e-3
  [../]
  [./polaryblock2]
    variable = polar_y
    type = RandomIC
    block = '2'
    min = -0.6e-3
    max = 0.6e-3
  [../]
[]

[AuxVariables]
  [./stress_xx_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xx_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./matl_e11]
    type = RankTwoAux
    block = '1 2 4'
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    variable = strain_xx_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e12]
    type = RankTwoAux
    block = '1 2 4'
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    variable = strain_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e13]
    type = RankTwoAux
    block = '1 2 4'
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 2
    variable = strain_xz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e22]
    type = RankTwoAux
    block = '1 2 4'
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    variable = strain_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e23]
    type = RankTwoAux
    block = '1 2 4'
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 2
    variable = strain_yz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e33]
    type = RankTwoAux
    block = '1 2 4'
    rank_two_tensor = elastic_strain
    index_i = 2
    index_j = 2
    variable = strain_zz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s11]
    type = RankTwoAux
    block = '1 2 4'
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = stress_xx_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s12]
    type = RankTwoAux
    block = '1 2 4'
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = stress_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s13]
    type = RankTwoAux
    block = '1 2 4'
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
    variable = stress_xz_elastic
    execute_on = 'timestep_end'
  [../]
 [./matl_s22]
    type = RankTwoAux
    block = '1 2 4'
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = stress_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s23]
    type = RankTwoAux
    block = '1 2 4'
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    variable = stress_yz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s33]
    type = RankTwoAux
    block = '1 2 4'
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    variable = stress_zz_elastic
    execute_on = 'timestep_end'
  [../]
[]

[Materials]
  [./vacmat]
    type = GenericConstantMaterial
    block = '3'
  [../]
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    block = '1 2'
    fill_method = symmetric9
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
  [../]
  [./strain_1]
    block = '1 2'
    type = ComputeSmallStrain
  [../]
  [./stress_1]
    block = '1 2'
    type = ComputeLinearElasticStress
  [../]
  [./elasticity_tensor_2]
    type = ComputeElasticityTensor
    block = '4'
    fill_method = symmetric9
    C_ijkl = '319 99.6 99.6 319 99.6 319 109.53 109.53 109.53'
  [../]
  [./strain_2]
    block = '4'
    type = ComputeSmallStrain
  [../]
  [./stress_2]
    block = '4'
    type = ComputeLinearElasticStress
  [../]
  [./slab_ferroelectric]
    type = ComputeElectrostrictiveTensor
    block = '1 2'
    Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.03375 0.03375 0.03375'
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
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
    block = '1 2'
    variable = disp_x
    component = 0
  [../]
  [./ferroelectriccouplingX_yy]
    type = FerroelectricCouplingX
    block = '1 2'
    variable = disp_y
    component = 1
  [../]
  [./ferroelectriccouplingX_zz]
    type = FerroelectricCouplingX
    block = '1 2'
    variable = disp_z
    component = 2
  [../]

  ##Electrostatics

  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
     block = '1 2'
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     block = '1 2'
     permittivity = 0.08854187
  [../]
  [./VAC_E_int]
     type = Electrostatics
     variable = potential_int
     block = '3'
     permittivity = 0.008854187
  [../]
  [./SUB_E_int]
     type = Electrostatics
     variable = potential_int
     block = '4'
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

  [./disp_x_sub]
    type = DirichletBC
    variable = disp_x
    value = 0.0
    boundary = '1'
  [../]
  [./disp_y_sub]
    type = DirichletBC
    variable = disp_y
    value = 0.0
    boundary = '1'
  [../]
  [./disp_z_sub]
    type = DirichletBC
    variable = disp_z
    value = 0.0
    boundary = '1'
  [../]

  [./potential_sub]
    type = DirichletBC
    variable = potential_int
    value = 0.0
    boundary = '8'
  [../]

  [./stress_dispx_left]
    type = StressBC
    variable = disp_x
    component = 0
    boundary_stress = '-2.0 0.0 0 0 0 0'
    boundary = '9'
  [../]

  [./stress_dispx_right]
    type = StressBC
    variable = disp_x
    component = 0
    boundary_stress = '-2.0 0.0 0 0 0 0'
    boundary = '10'
  [../]


 [./Periodic]
    #FE section
    [./FE_disp_x_pbc]
      variable = disp_x
      primary = '2'
      secondary = '3'
      translation = '0 40 0'
    [../]
    [./FE_disp_y_pbc]
      variable = disp_y
      primary = '2'
      secondary = '3'
      translation = '0 40 0'
    [../]
    [./FE_disp_z_pbc]
      variable = disp_z
      primary = '2'
      secondary = '3'
      translation = '0 40 0'
    [../]
    [./FE_polar_x_pbc]
      variable = polar_x
      primary = '2'
      secondary = '3'
      translation = '0 40 0'
    [../]
    [./FE_polar_y_pbc]
      variable = polar_y
      primary = '2'
      secondary = '3'
      translation = '0 40 0'
    [../]
    [./FE_polar_z_pbc]
      variable = polar_z
      primary = '2'
      secondary = '3'
      translation = '0 40 0'
    [../]
    [./FE_potential_int_pbc]
      variable = potential_int
      primary = '2'
      secondary = '3'
      translation = '0 40 0'
    [../]

    #Substrate section
    [./Sub_disp_x_pbc]
      variable = disp_x
      primary = '6'
      secondary = '7'
      translation = '0 -40 0'
    [../]
    [./Sub_disp_y_pbc]
      variable = disp_y
      primary = '6'
      secondary = '7'
      translation = '0 -40 0'
    [../]
    [./Sub_disp_z_pbc]
      variable = disp_z
      primary = '6'
      secondary = '7'
      translation = '0 -40 0'
   [../]

    #Vacuum section
    [./Vac_potential_int_pbc_1]
      variable = potential_int
      primary = '4'
      secondary = '5'
      translation = '0 -40 0'
    [../]
 [../]
[]

[Postprocessors]
   [./Fbulk]
      type = BulkEnergy
      block = '1 2'
      execute_on = 'timestep_end'
    [../]
    [./Fwall]
      type = WallEnergy
      block = '1 2'
      execute_on = 'timestep_end'
    [../]
    [./Felastic]
      type = ElasticEnergy
      block = '1 2'
      execute_on = 'timestep_end'
    [../]
    [./Fcoupled]
      type = CoupledEnergy
      block = '1 2'
      execute_on = 'timestep_end'
    [../]
    [./Felec]
      type = ElectrostaticEnergy
      block = '1 2'
      execute_on = 'timestep_end'
    [../]
    [./Ftotal]
      type = TotalEnergyFlow
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
  expression = 'perc_change <= 1.0e-3' #this seems to be the max percent change before numerical noise
 [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121                1e-6      1e-8    bjacobi'
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
  #dt = 0.5
  dtmin = 1e-13
  dtmax = 0.1
  #num_steps = 2
[]

[Debug]
  show_var_residual_norms = true
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_aPTOridge_test_const_stress1p_lowPFull
    elemental_as_nodal = true
    interval = 1
  [../]
  [./outCSV]
    type = CSV
    file_base = out_aPTOridge_test_const_stress1p_lowPFull
    interval = 1
  [../]
[]
