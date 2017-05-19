
[Mesh]
  file = exodus_trilayer_8_8_4.e
[]

[GlobalParams]
  #NOTE: We use a nanometer, nanonewton, attocoulomb unit system

  len_scale = 1.0
  alpha1 = -0.1722883 # T = 298 K
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

  #miscellaneous GlobalParams options
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  displacements = 'disp_x disp_y disp_z'

  #applied misfit strain, negative compression
  prefactor = -0.008

  q = 0.16 # charge of the electron in attocoulombs
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

  [./n]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]

  [./p]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
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

  #semiconducting charge carriers (store their values)
  [./nm]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pp]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./NAm]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rho]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]


[AuxKernels]
  [./matl_e11]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    variable = strain_xx_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    variable = strain_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e13]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 2
    variable = strain_xz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    variable = strain_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e23]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 2
    variable = strain_yz_elastic
  #this is a comment
    execute_on = 'timestep_end'
  [../]
  [./matl_e33]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 2
    index_j = 2
    variable = strain_zz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s11]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = stress_xx_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s12]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = stress_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s13]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
    variable = stress_xz_elastic
    execute_on = 'timestep_end'
  [../]
 [./matl_s22]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = stress_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s23]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    variable = stress_yz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s33]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    variable = stress_zz_elastic
    execute_on = 'timestep_end'
  [../]

  #calculate the semiconducting charge carriers
  [./nm_calc]
    type = SemiconductingChargeCarriersAux
    charge_type = 0
    variable = nm
    kT = 0.0041124
    q = -0.16
    NA = 1e-7 #1e20 / m^3
    NC = 0.000666
    NV = 0.015
    EA = -0.00801
    EC = -0.66483
    EV = -0.71289 #-4.45 eV
    EF = -0.71289 #~-4.45 eV
    block  = '2'
    execute_on = 'timestep_end'
  [../]
  [./pp_calc]
    type = SemiconductingChargeCarriersAux
    charge_type = 1
    variable = pp
    kT = 0.0041124
    q = -0.16
    NA = 1e-7 #1e20 / m^3
    NC = 0.000666
    NV = 0.015
    EA = -0.00801
    EC = -0.66483
    EV = -0.71289 #-4.45 eV
    EF = -0.71289 #~-4.45 eV
    block  = '2'
    execute_on = 'timestep_end'
  [../]
  [./NAm_calc]
    type = SemiconductingChargeCarriersAux
    charge_type = 2
    variable = NAm
    kT = 0.0041124
    q = -0.16
    NA = 1e-7 #1e20 / m^3
    NC = 0.000666
    NV = 0.015
    EA = -0.00801
    EC = -0.66483
    EV = -0.71289 #-4.45 eV
    EF = -0.71289 #~-4.45 eV
    block  = '2'
    execute_on = 'timestep_end'
  [../]
  [./rho_calc]
    type = SemiconductingChargeCarriersAux
    charge_type = 3
    variable = rho
    kT = 0.0041124
    q = -0.16
    NA = 1e-7 #1e20 / m^3
    NC = 0.000666
    NV = 0.015
    EA = -0.00801
    EC = -0.66483
    EV = -0.71289 #-4.45 eV
    EF = -0.71289 #~-4.45 eV
    block  = '2'
    execute_on = 'timestep_end'
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
    C_ijkl = '281 115.74 115.74 281 115.74 281 97.18 97.18 97.18' #PTO from MaterialsProject
  [../]
  [./elasticity_tensor_2]
    type = ComputeElasticityTensor
    C_ijkl = '319 99.6 99.6 319 99.6 319 109.53 109.53 109.53' #STO from MaterialsProject
    fill_method = symmetric9
    block = '2'
  [../]
  [./strain_2]
    type = ComputeSmallStrain
    block = '2'
  [../]
  [./stress_2]
    type = ComputeLinearElasticStress
    block = '2'
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
     #permittivity = 0.0575522155 #er = 6.5
     block = '1'
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     block = '1'
     permittivity = 0.0575522155 #er = 6.5
  [../]

  [./semi_E_int]
     type = Electrostatics
     variable = potential_int
     block  = '2'
     permittivity = 0.44270935 #er = 50
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
  [./bot_disp_x]
    variable = disp_x
    type = DirichletBC
    value = 0
    boundary = '1'
  [../]
  [./bot_disp_y]
    variable = disp_y
    type = DirichletBC
    value = 0
    boundary = '1'
  [../]
  [./bot_disp_z]
    variable = disp_z
    type = DirichletBC
    value = 0
    boundary = '1'
  [../]

  [./bot_potential_int]
    variable = potential_int
    type = DirichletBC
    value = 0.00016
    boundary = '1'
  [../]

  [./top_potential_int]
    variable = potential_int
    type = DirichletBC
    value = 0.00016
    boundary = '2'
  [../]


  [./Periodic]
    [./TB_FE_disp_x_pbc]
      variable = disp_x
      primary = '3'
      secondary = '4'
      translation = '0 8 0'
    [../]
    [./TB_FE_disp_y_pbc]
      variable = disp_y
      primary = '3'
      secondary = '4'
      translation = '0 8 0'
    [../]
    [./TB_FE_disp_z_pbc]
      variable = disp_z
      primary = '3'
      secondary = '4'
      translation = '0 8 0'
    [../]

    [./TB_FE_polar_x_pbc]
      variable = polar_x
      primary = '3'
      secondary = '4'
      translation = '0 8 0'
    [../]
    [./TB_FE_polar_y_pbc]
      variable = polar_y
      primary = '3'
      secondary = '4'
      translation = '0 8 0'
    [../]
    [./TB_FE_polar_z_pbc]
      variable = polar_z
      primary = '3'
      secondary = '4'
      translation = '0 8 0'
    [../]
    [./TB_FE_potential_int_pbc]
      variable = potential_int
      primary = '3'
      secondary = '4'
      translation = '0 8 0'
    [../]
  #
    [./TBsemi_disp_x_pbc]
      variable = disp_x
      primary = '5'
      secondary = '6'
      translation = '8 0 0'
    [../]
    [./TBsemi_disp_y_pbc]
      variable = disp_y
      primary = '5'
      secondary = '6'
      translation = '8 0 0'
    [../]
    [./TBsemi_disp_z_pbc]
      variable = disp_z
      primary = '5'
      secondary = '6'
      translation = '8 0 0'
    [../]
    [./TBsemi_potential_int_pbc]
      variable = potential_int
      primary = '5'
      secondary = '6'
      translation = '8 0 0'
    [../]

    [./RL_FE_disp_x_pbc]
      variable = disp_x
      primary = '5'
      secondary = '6'
      translation = '8 0 0'
    [../]
    [./RL_FE_disp_y_pbc]
      variable = disp_y
      primary = '5'
      secondary = '6'
      translation = '8 0 0'
    [../]
    [./RL_FE_disp_z_pbc]
      variable = disp_z
      primary = '5'
      secondary = '6'
      translation = '8 0 0'
    [../]

    [./RL_FE_polar_x_pbc]
      variable = polar_x
      primary = '5'
      secondary = '6'
      translation = '8 0 0'
    [../]
    [./RL_FE_polar_y_pbc]
      variable = polar_y
      primary = '5'
      secondary = '6'
      translation = '8 0 0'
    [../]
    [./RL_FE_polar_z_pbc]
      variable = polar_z
      primary = '5'
      secondary = '6'
      translation = '8 0 0'
    [../]
    [./RL_FE_potential_int_pbc]
      variable = potential_int
      primary = '5'
      secondary = '6'
      translation = '8 0 0'
    [../]

    [./RLsemi_disp_x_pbc]
      variable = disp_x
      primary = '11'
      secondary = '12'
      translation = '8 0 0'
    [../]
    [./RLsemi_disp_y_pbc]
      variable = disp_y
      primary = '11'
      secondary = '12'
      translation = '8 0 0'
    [../]
    [./RLsub_disp_z_pbc]
      variable = disp_z
      primary = '11'
      secondary = '12'
      translation = '8 0 0'
    [../]
    [./RLsub_potential_int_pbc]
      variable = potential_int
      primary = '11'
      secondary = '12'
      translation = '8 0 0'
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
    [./Felastic]
      type = ElasticEnergy
      block = '1 2'
      execute_on = 'timestep_end'
    [../]
    [./Fcoupled]
      block = '1'
      type = CoupledEnergy
      execute_on = 'timestep_end'
    [../]
    [./Felec]
      block = '1'
      type = ElectrostaticEnergy
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


#[UserObjects]
# [./kill]
#  type = Terminator
#  expression = 'perc_change <= 1.0e-3'
# [../]
#[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart  -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '        200            1e-10     1e-8      1e-4      bjacobi'
  [../]
[]

[Executioner]

  type = Transient

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.25
    optimal_iterations = 6 #should be 5 probably
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.85
  [../]

  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  #dt = 0.5
  dtmin = 1e-13
  dtmax = 0.85
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = outPTO_trilayer_085_8_8_4_c01
    elemental_as_nodal = true
    interval = 1
  [../]
[]
