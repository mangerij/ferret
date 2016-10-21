
[Mesh]
  file = exodus_thinfilm_test_085_2_2_1.e
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
  G11/G110 = 1.0
  G12/G110 = 0
  G44/G110 = 0.5
  G44P/G110 = 0.5
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  use_displaced_mesh = false
  displacements = 'disp_x disp_y disp_z'
  prefactor = -0.01 #negative = tension, positive = compression
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
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 1e-4
      max = 1e-4
    [../]
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 1e-4
      max = 1e-4
    [../]
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 1e-4
      max = 1e-4
    [../]
  [../]
[]


[AuxVariables]
  [./strain_xx_elastic_old]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy_elastic_old]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xy_elastic_old]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xz_elastic_old]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zz_elastic_old]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yz_elastic_old]
    order = CONSTANT
    family = MONOMIAL
  [../]


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
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    variable = strain_xx_elastic
  [../]
  [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    variable = strain_xy_elastic
  [../]
  [./matl_e13]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 2
    variable = strain_xz_elastic
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    variable = strain_yy_elastic
  [../]
  [./matl_e23]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 2
    variable = strain_yz_elastic
  [../]
  [./matl_e33]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 2
    index_j = 2
    variable = strain_zz_elastic
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


  [./matl_e11_old]
    type = OldVar
    variable = strain_xx_elastic_old
  [../]
  [./matl_e12_old]
    type = OldVar
    variable = strain_xy_elastic_old
  [../]
  [./matl_e13_old]
    type = OldVar
    variable = strain_xz_elastic_old
  [../]
  [./matl_e22_old]
    type = OldVar
    variable = strain_yy_elastic_old
  [../]
  [./matl_e23_old]
    type = OldVar
    variable = strain_yz_elastic_old
  [../]
  [./matl_e33_old]
    type = OldVar
    variable = strain_zz_elastic_old
  [../]
[]

[Materials]
  [./eigen_strain_zz] #Use for stress-free strain (ie epitaxial)
    type = ComputeEigenstrain
    block = '1'
   # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
    eigen_base = '1 0 0 0 1 0 0 0 0'
  [../]

  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
    block = '1'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    block = '1'
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

  [./elasticity_tensor_2]
    type = ComputeElasticityTensor
    C_ijkl = '319 99.6 99.6 319 99.6 319 109.53 109.53 109.53'
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


 #------------------------#
  #-Ferroelectric problem -#
  #------------------------#

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

  [./potential_int_2]
    type = DirichletBC
    variable = potential_int
    boundary = '7'
    value = 0.0
  [../]

  [./bot_disp_x]
    variable = disp_x
    type = DirichletBC
    value = 0
    boundary = '7'
  [../]
  [./bot_disp_y]
    variable = disp_y
    type = DirichletBC
    value = 0
    boundary = '7'
  [../]
  [./bot_disp_z]
    variable = disp_z
    type = DirichletBC
    value = 0
    boundary = '7'
  [../]

###The many BCs -- so many BCs####

  [./TB_disp_x_3_match_PBC_x]
    variable = disp_x
    type = MatchedGradValueBC
    component = 0
    boundary = '3'
    variable = disp_x
    u = disp_x
    v = strain_xx_elastic_old
  [../]

  [./TB_disp_x_5_match_PBC_x]
    variable = disp_x
    type = MatchedGradValueBC
    component = 0
    boundary = '5'
    variable = disp_x
    u = disp_x
    v = strain_xx_elastic_old
  [../]

  [./TB_disp_y_3_match_PBC_x]
    variable = disp_y
    type = MatchedGradValueBC
    component = 0
    boundary = '3'
    variable = disp_y
    u = disp_y
    v = strain_xy_elastic_old
  [../]

  [./TB_disp_y_5_match_PBC_x]
    variable = disp_y
    type = MatchedGradValueBC
    component = 0
    boundary = '5'
    variable = disp_y
    u = disp_y
    v = strain_xy_elastic_old
  [../]

  [./TB_disp_z_3_match_PBC_x]
    variable = disp_z
    type = MatchedGradValueBC
    component = 0
    boundary = '3'
    variable = disp_z
    u = disp_z
    v = strain_xz_elastic_old
  [../]

  [./TB_disp_z_5_match_PBC_x]
    variable = disp_z
    type = MatchedGradValueBC
    component = 0
    boundary = '5'
    variable = disp_z
    u = disp_z
    v = strain_xz_elastic_old
  [../]
#
  [./TB_disp_x_3_match_PBC_y]
    variable = disp_x
    type = MatchedGradValueBC
    component = 1
    boundary = '3'
    variable = disp_x
    u = disp_x
    v = strain_xy_elastic_old
  [../]

  [./TB_disp_x_5_match_PBC_y]
    variable = disp_x
    type = MatchedGradValueBC
    component = 1
    boundary = '5'
    variable = disp_x
    u = disp_x
    v = strain_xy_elastic_old
  [../]

  [./TB_disp_y_3_match_PBC_y]
    variable = disp_y
    type = MatchedGradValueBC
    component = 1
    boundary = '3'
    variable = disp_y
    u = disp_y
    v = strain_yy_elastic_old
  [../]

  [./TB_disp_y_5_match_PBC_y]
    variable = disp_y
    type = MatchedGradValueBC
    component = 1
    boundary = '5'
    variable = disp_y
    u = disp_y
    v = strain_yy_elastic_old
  [../]

  [./TB_disp_z_3_match_PBC_y]
    variable = disp_z
    type = MatchedGradValueBC
    component = 1
    boundary = '3'
    variable = disp_z
    u = disp_z
    v = strain_yz_elastic_old
  [../]

  [./TB_disp_z_5_match_PBC_y]
    variable = disp_z
    type = MatchedGradValueBC
    component = 1
    boundary = '5'
    variable = disp_z
    u = disp_z
    v = strain_yz_elastic_old
  [../]
#
  [./TB_disp_x_3_match_PBC_z]
    variable = disp_x
    type = MatchedGradValueBC
    component = 2
    boundary = '3'
    variable = disp_x
    u = disp_x
    v = strain_xz_elastic_old
  [../]

  [./TB_disp_x_5_match_PBC_z]
    variable = disp_x
    type = MatchedGradValueBC
    component = 2
    boundary = '5'
    variable = disp_x
    u = disp_x
    v = strain_xz_elastic_old
  [../]

  [./TB_disp_y_3_match_PBC_z]
    variable = disp_y
    type = MatchedGradValueBC
    component = 2
    boundary = '3'
    variable = disp_y
    u = disp_y
    v = strain_yz_elastic_old
  [../]

  [./TB_disp_y_5_match_PBC_z]
    variable = disp_y
    type = MatchedGradValueBC
    component = 2
    boundary = '5'
    variable = disp_y
    u = disp_y
    v = strain_yz_elastic_old
  [../]

  [./TB_disp_z_3_match_PBC_z]
    variable = disp_z
    type = MatchedGradValueBC
    component = 2
    boundary = '3'
    variable = disp_z
    u = disp_z
    v = strain_zz_elastic_old
  [../]

  [./TB_disp_z_5_match_PBC_z]
    variable = disp_z
    type = MatchedGradValueBC
    component = 2
    boundary = '5'
    variable = disp_z
    u = disp_z
    v = strain_zz_elastic_old
  [../]

#Next surfaces


  [./RL_disp_x_4_match_PBC_x]
    variable = disp_x
    type = MatchedGradValueBC
    component = 0
    boundary = '4'
    variable = disp_x
    u = disp_x
    v = strain_xx_elastic_old
  [../]

  [./RL_disp_x_6_match_PBC_x]
    variable = disp_x
    type = MatchedGradValueBC
    component = 0
    boundary = '6'
    variable = disp_x
    u = disp_x
    v = strain_xx_elastic_old
  [../]

  [./RL_disp_y_4_match_PBC_x]
    variable = disp_y
    type = MatchedGradValueBC
    component = 0
    boundary = '4'
    variable = disp_y
    u = disp_y
    v = strain_xy_elastic_old
  [../]

  [./RL_disp_y_6_match_PBC_x]
    variable = disp_y
    type = MatchedGradValueBC
    component = 0
    boundary = '6'
    variable = disp_y
    u = disp_y
    v = strain_xy_elastic_old
  [../]

  [./RL_disp_z_4_match_PBC_x]
    variable = disp_z
    type = MatchedGradValueBC
    component = 0
    boundary = '4'
    variable = disp_z
    u = disp_z
    v = strain_xz_elastic_old
  [../]

  [./RL_disp_z_6_match_PBC_x]
    variable = disp_z
    type = MatchedGradValueBC
    component = 0
    boundary = '6'
    variable = disp_z
    u = disp_z
    v = strain_xz_elastic_old
  [../]
#
  [./RL_disp_x_4_match_PBC_y]
    variable = disp_x
    type = MatchedGradValueBC
    component = 1
    boundary = '4'
    variable = disp_x
    u = disp_x
    v = strain_xy_elastic_old
  [../]

  [./RL_disp_x_6_match_PBC_y]
    variable = disp_x
    type = MatchedGradValueBC
    component = 1
    boundary = '6'
    variable = disp_x
    u = disp_x
    v = strain_xy_elastic_old
  [../]

  [./RL_disp_y_4_match_PBC_y]
    variable = disp_y
    type = MatchedGradValueBC
    component = 1
    boundary = '4'
    variable = disp_y
    u = disp_y
    v = strain_yy_elastic_old
  [../]

  [./RL_disp_y_6_match_PBC_y]
    variable = disp_y
    type = MatchedGradValueBC
    component = 1
    boundary = '6'
    variable = disp_y
    u = disp_y
    v = strain_yy_elastic_old
  [../]

  [./RL_disp_z_4_match_PBC_y]
    variable = disp_z
    type = MatchedGradValueBC
    component = 1
    boundary = '4'
    variable = disp_z
    u = disp_z
    v = strain_yz_elastic_old
  [../]

  [./RL_disp_z_6_match_PBC_y]
    variable = disp_z
    type = MatchedGradValueBC
    component = 1
    boundary = '6'
    variable = disp_z
    u = disp_z
    v = strain_yz_elastic_old
  [../]
#
  [./RL_disp_x_4_match_PBC_z]
    variable = disp_x
    type = MatchedGradValueBC
    component = 2
    boundary = '4'
    variable = disp_x
    u = disp_x
    v = strain_xz_elastic_old
  [../]

  [./RL_disp_x_6_match_PBC_z]
    variable = disp_x
    type = MatchedGradValueBC
    component = 2
    boundary = '6'
    variable = disp_x
    u = disp_x
    v = strain_xz_elastic_old
  [../]

  [./RL_disp_y_4_match_PBC_z]
    variable = disp_y
    type = MatchedGradValueBC
    component = 2
    boundary = '4'
    variable = disp_y
    u = disp_y
    v = strain_yz_elastic_old
  [../]

  [./RL_disp_y_6_match_PBC_z]
    variable = disp_y
    type = MatchedGradValueBC
    component = 2
    boundary = '6'
    variable = disp_y
    u = disp_y
    v = strain_yz_elastic_old
  [../]

  [./RL_disp_z_4_match_PBC_z]
    variable = disp_z
    type = MatchedGradValueBC
    component = 2
    boundary = '4'
    variable = disp_z
    u = disp_z
    v = strain_zz_elastic_old
  [../]

  [./RL_disp_z_6_match_PBC_z]
    variable = disp_z
    type = MatchedGradValueBC
    component = 2
    boundary = '6'
    variable = disp_z
    u = disp_z
    v = strain_zz_elastic_old
  [../]

###############substrate

  [./TBsub_disp_x_8_match_PBC_x]
    variable = disp_x
    type = MatchedGradValueBC
    component = 0
    boundary = '8'
    variable = disp_x
    u = disp_x
    v = strain_xx_elastic_old
  [../]

  [./TBsub_disp_x_10_match_PBC_x]
    variable = disp_x
    type = MatchedGradValueBC
    component = 0
    boundary = '10'
    variable = disp_x
    u = disp_x
    v = strain_xx_elastic_old
  [../]

  [./TBsub_disp_y_8_match_PBC_x]
    variable = disp_y
    type = MatchedGradValueBC
    component = 0
    boundary = '8'
    variable = disp_y
    u = disp_y
    v = strain_xy_elastic_old
  [../]

  [./TBsub_disp_y_10_match_PBC_x]
    variable = disp_y
    type = MatchedGradValueBC
    component = 0
    boundary = '10'
    variable = disp_y
    u = disp_y
    v = strain_xy_elastic_old
  [../]

  [./TBsub_disp_z_8_match_PBC_x]
    variable = disp_z
    type = MatchedGradValueBC
    component = 0
    boundary = '8'
    variable = disp_z
    u = disp_z
    v = strain_xz_elastic_old
  [../]

  [./TBsub_disp_z_10_match_PBC_x]
    variable = disp_z
    type = MatchedGradValueBC
    component = 0
    boundary = '10'
    variable = disp_z
    u = disp_z
    v = strain_xz_elastic_old
  [../]
#
  [./TBsub_disp_x_8_match_PBC_y]
    variable = disp_x
    type = MatchedGradValueBC
    component = 1
    boundary = '8'
    variable = disp_x
    u = disp_x
    v = strain_xy_elastic_old
  [../]

  [./TBsub_disp_x_10_match_PBC_y]
    variable = disp_x
    type = MatchedGradValueBC
    component = 1
    boundary = '10'
    variable = disp_x
    u = disp_x
    v = strain_xy_elastic_old
  [../]

  [./TBsub_disp_y_8_match_PBC_y]
    variable = disp_y
    type = MatchedGradValueBC
    component = 1
    boundary = '8'
    variable = disp_y
    u = disp_y
    v = strain_yy_elastic_old
  [../]

  [./TBsub_disp_y_10_match_PBC_y]
    variable = disp_y
    type = MatchedGradValueBC
    component = 1
    boundary = '10'
    variable = disp_y
    u = disp_y
    v = strain_yy_elastic_old
  [../]

  [./TBsub_disp_z_8_match_PBC_y]
    variable = disp_z
    type = MatchedGradValueBC
    component = 1
    boundary = '8'
    variable = disp_z
    u = disp_z
    v = strain_yz_elastic_old
  [../]

  [./TBsub_disp_z_10_match_PBC_y]
    variable = disp_z
    type = MatchedGradValueBC
    component = 1
    boundary = '10'
    variable = disp_z
    u = disp_z
    v = strain_yz_elastic_old
  [../]
#
  [./TBsub_disp_x_8_match_PBC_z]
    variable = disp_x
    type = MatchedGradValueBC
    component = 2
    boundary = '8'
    variable = disp_x
    u = disp_x
    v = strain_xz_elastic_old
  [../]

  [./TBsub_disp_x_10_match_PBC_z]
    variable = disp_x
    type = MatchedGradValueBC
    component = 2
    boundary = '10'
    variable = disp_x
    u = disp_x
    v = strain_xz_elastic_old
  [../]

  [./TBsub_disp_y_8_match_PBC_z]
    variable = disp_y
    type = MatchedGradValueBC
    component = 2
    boundary = '8'
    variable = disp_y
    u = disp_y
    v = strain_yz_elastic_old
  [../]

  [./TBsub_disp_y_10_match_PBC_z]
    variable = disp_y
    type = MatchedGradValueBC
    component = 2
    boundary = '10'
    variable = disp_y
    u = disp_y
    v = strain_yz_elastic_old
  [../]

  [./TBsub_disp_z_8_match_PBC_z]
    variable = disp_z
    type = MatchedGradValueBC
    component = 2
    boundary = '8'
    variable = disp_z
    u = disp_z
    v = strain_zz_elastic_old
  [../]

  [./TBsub_disp_z_10_match_PBC_z]
    variable = disp_z
    type = MatchedGradValueBC
    component = 2
    boundary = '10'
    variable = disp_z
    u = disp_z
    v = strain_zz_elastic_old
  [../]

#Next surfaces


  [./RLsub_disp_x_9_match_PBC_x]
    variable = disp_x
    type = MatchedGradValueBC
    component = 0
    boundary = '9'
    variable = disp_x
    u = disp_x
    v = strain_xx_elastic_old
  [../]

  [./RLsub_disp_x_11_match_PBC_x]
    variable = disp_x
    type = MatchedGradValueBC
    component = 0
    boundary = '11'
    variable = disp_x
    u = disp_x
    v = strain_xx_elastic_old
  [../]

  [./RLsub_disp_y_9_match_PBC_x]
    variable = disp_y
    type = MatchedGradValueBC
    component = 0
    boundary = '9'
    variable = disp_y
    u = disp_y
    v = strain_xy_elastic_old
  [../]

  [./RLsub_disp_y_11_match_PBC_x]
    variable = disp_y
    type = MatchedGradValueBC
    component = 0
    boundary = '11'
    variable = disp_y
    u = disp_y
    v = strain_xy_elastic_old
  [../]

  [./RLsub_disp_z_9_match_PBC_x]
    variable = disp_z
    type = MatchedGradValueBC
    component = 0
    boundary = '9'
    variable = disp_z
    u = disp_z
    v = strain_xz_elastic_old
  [../]

  [./RLsub_disp_z_11_match_PBC_x]
    variable = disp_z
    type = MatchedGradValueBC
    component = 0
    boundary = '11'
    variable = disp_z
    u = disp_z
    v = strain_xz_elastic_old
  [../]
#
  [./RLsub_disp_x_9_match_PBC_y]
    variable = disp_x
    type = MatchedGradValueBC
    component = 1
    boundary = '9'
    variable = disp_x
    u = disp_x
    v = strain_xy_elastic_old
  [../]

  [./RLsub_disp_x_11_match_PBC_y]
    variable = disp_x
    type = MatchedGradValueBC
    component = 1
    boundary = '11'
    variable = disp_x
    u = disp_x
    v = strain_xy_elastic_old
  [../]

  [./RLsub_disp_y_9_match_PBC_y]
    variable = disp_y
    type = MatchedGradValueBC
    component = 1
    boundary = '9'
    variable = disp_y
    u = disp_y
    v = strain_yy_elastic_old
  [../]

  [./RLsub_disp_y_11_match_PBC_y]
    variable = disp_y
    type = MatchedGradValueBC
    component = 1
    boundary = '11'
    variable = disp_y
    u = disp_y
    v = strain_yy_elastic_old
  [../]

  [./RLsub_disp_z_9_match_PBC_y]
    variable = disp_z
    type = MatchedGradValueBC
    component = 1
    boundary = '9'
    variable = disp_z
    u = disp_z
    v = strain_yz_elastic_old
  [../]

  [./RLsub_disp_z_11_match_PBC_y]
    variable = disp_z
    type = MatchedGradValueBC
    component = 1
    boundary = '11'
    variable = disp_z
    u = disp_z
    v = strain_yz_elastic_old
  [../]
#
  [./RLsub_disp_x_9_match_PBC_z]
    variable = disp_x
    type = MatchedGradValueBC
    component = 2
    boundary = '9'
    variable = disp_x
    u = disp_x
    v = strain_xz_elastic_old
  [../]

  [./RLsub_disp_x_11_match_PBC_z]
    variable = disp_x
    type = MatchedGradValueBC
    component = 2
    boundary = '11'
    variable = disp_x
    u = disp_x
    v = strain_xz_elastic_old
  [../]

  [./RLsub_disp_y_9_match_PBC_z]
    variable = disp_y
    type = MatchedGradValueBC
    component = 2
    boundary = '9'
    variable = disp_y
    u = disp_y
    v = strain_yz_elastic_old
  [../]

  [./RLsub_disp_y_11_match_PBC_z]
    variable = disp_y
    type = MatchedGradValueBC
    component = 2
    boundary = '11'
    variable = disp_y
    u = disp_y
    v = strain_yz_elastic_old
  [../]

  [./RLsub_disp_z_9_match_PBC_z]
    variable = disp_z
    type = MatchedGradValueBC
    component = 2
    boundary = '9'
    variable = disp_z
    u = disp_z
    v = strain_zz_elastic_old
  [../]

  [./RLsub_disp_z_11_match_PBC_z]
    variable = disp_z
    type = MatchedGradValueBC
    component = 2
    boundary = '11'
    variable = disp_z
    u = disp_z
    v = strain_zz_elastic_old
  [../]

  [./Periodic]
    [./TB_polar_x_pbc]
      variable = polar_x
      primary = '3'
      secondary = '5'
      translation = '0 2 0'
    [../]
    [./TB_polar_y_pbc]
      variable = polar_y
      primary = '3'
      secondary = '5'
      translation = '0 2 0'
    [../]
    [./TB_polar_z_pbc]
      variable = polar_z
      primary = '3'
      secondary = '5'
      translation = '0 2 0'
    [../]
    [./TB_potential_int_pbc]
      variable = potential_int
      primary = '3'
      secondary = '5'
      translation = '0 2 0'
    [../]

    [./RL_polar_x_pbc]
      variable = polar_x
      primary = '4'
      secondary = '6'
      translation = '2 0 0'
    [../]
    [./RL_polar_y_pbc]
      variable = polar_y
      primary = '4'
      secondary = '6'
      translation = '2 0 0'
    [../]
    [./RL_polar_z_pbc]
      variable = polar_z
      primary = '4'
      secondary = '6'
      translation = '2 0 0'
    [../]
    [./RL_potential_int_pbc]
      variable = potential_int
      primary = '4'
      secondary = '6'
      translation = '2 0 0'
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
     permittivity = 0.08854187
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
    [./perc_change]
     type = PercentChangePostprocessor
     postprocessor = Felastic
   [../]
   [./num_NLin]
    type = NumNonlinearIterations
   [../]
   [./num_Lin]
    type = NumLinearIterations
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
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121                1e-6      1e-8    bjacobi'
    #petsc_options_value = '    121                1e-8      1e-5    bjacobi'
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
  num_steps = 80
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_FEthinfilm_09_20_20_10_c01_withPBC
    elemental_as_nodal = true
    interval = 1
  [../]
[]
