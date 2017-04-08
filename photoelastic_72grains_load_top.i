
[Mesh]
  file = 3D_HCP_32_r100.e
[]

[MeshModifiers]
  [./add_side_sets]
    type = SideSetsFromNormals
    normals = '1  0  0
               0  0  1
               0  0  -1'
    new_boundary = 'left side1 side2'
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
  n_a = 2.437
  n_b = 2.437
  n_g = 2.365
[]


[Variables]
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

  [./dn_1]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./dn_2]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./dn_3]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./pfdn_1]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./pfdn_2]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./pfdn_3]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./dn_bire_12]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./dn_bire_23]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./dn_bire_13]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./n_1]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./n_2]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./n_3]
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

  [./dn_s1]
    type = ChangeInRefractiveIndex
    index_i = 0
    index_j = 0
    index_k = 0
    index_l = 0
    variable = dn_1
    execute_on = 'timestep_end'
  [../]

  [./dn_s2]
    type = ChangeInRefractiveIndex
    index_i = 1
    index_j = 1
    index_k = 1
    index_l = 1
    variable = dn_2
    execute_on = 'timestep_end'
  [../]

  [./dn_s3]
    type = ChangeInRefractiveIndex
    index_i = 2
    index_j = 2
    index_k = 2
    index_l = 2
    variable = dn_3
    execute_on = 'timestep_end'
  [../]

  [./pfdn_s1]
    type = ChangeInRefractiveIndex
    index_i = 0
    index_j = 0
    index_k = 0
    index_l = 0
    variable = pfdn_1
    execute_on = 'timestep_end'
  [../]

  [./pfdn_s2]
    type = ChangeInRefractiveIndex
    index_i = 1
    index_j = 1
    index_k = 1
    index_l = 1
    variable = pfdn_2
    execute_on = 'timestep_end'
  [../]

  [./pfdn_s3]
    type = ChangeInRefractiveIndex
    index_i = 2
    index_j = 2
    index_k = 2
    index_l = 2
    variable = pfdn_3
    execute_on = 'timestep_end'
  [../]

  [./dn_bire_s12]
    type = Birefringence
    variable = dn_bire_12
    per1 = n_1
    per2 = n_2
    execute_on = 'timestep_end'
  [../]

  [./dn_bire_s23]
    type = Birefringence
    variable = dn_bire_23
    per1 = n_2
    per2 = n_3
    execute_on = 'timestep_end'
  [../]

  [./dn_bire_s13]
    type = Birefringence
    variable = dn_bire_13
    per1 = n_1
    per2 = n_3
    execute_on = 'timestep_end'
  [../]


  [./n_1_c]
    type = RefractiveIndex
    variable = n_1
    component = 0
    var1 = dn_1
    execute_on = 'timestep_end'
  [../]

  [./n_2_c]
    type = RefractiveIndex
    variable = n_2
    component = 1
    var1 = dn_2
    execute_on = 'timestep_end'
  [../]

  [./n_3_c]
    type = RefractiveIndex
    variable = n_3
    component = 2
    var1 = dn_3
    execute_on = 'timestep_end'
  [../]

[]

################################################
# Block list:                                  #
#                                              #
# No 99?                                       #
#                                              #
################################################


[Materials]
  [./eigen_strain_zz] #Use for stress-free strain (ie epitaxial)
    type = ComputeEigenstrain
    block = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32'
    # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
    eigen_base = '1 0 0 0 1 0 0 0 0'
    eigenstrain_name = eigenstrain
    prefactor = 0.0
  [../]

  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 0.0
    euler_angle_2 = 15.0
    euler_angle_3 = 65.0
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
  [./photoelastic_tensor_1]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 0.0
    euler_angle_2 = 15.0
    euler_angle_3 = 65.0
    block = '1'
  [../]
  [./beta_tensor_1]
    type = ComputeBetaTensor
    block = '1'
    euler_angle_1 = 0.0
    euler_angle_2 = 15.0
    euler_angle_3 = 65.0
  [../]
  [./delta_beta_tensor_1]
    type = ComputeDeltaBetaTensor
    block = '1'
  [../]


  [./elasticity_tensor_2]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 12.0
    euler_angle_2 = -15.0
    euler_angle_3 = 65.0
    block = '2'
  [../]
  [./strain_2]
    type = ComputeSmallStrain
    block = '2'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_2]
    type = ComputeLinearElasticStress
    block = '2'
  [../]
  [./photoelastic_tensor_2]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 12.0
    euler_angle_2 = -15.0
    euler_angle_3 = 65.0
    block = '2'
  [../]
  [./beta_tensor_2]
    type = ComputeBetaTensor
    block = '2'
    euler_angle_1 = 12.0
    euler_angle_2 = -15.0
    euler_angle_3 = 65.0
  [../]
  [./delta_beta_tensor_2]
    type = ComputeDeltaBetaTensor
    block = '2'
  [../]



  [./elasticity_tensor_3]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 1.0
    euler_angle_2 = 50.0
    euler_angle_3 = -132.0
    block = '3'
  [../]
  [./strain_3]
    type = ComputeSmallStrain
    block = '3'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_3]
    type = ComputeLinearElasticStress
    block = '3'
  [../]
  [./photoelastic_tensor_3]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 1.0
    euler_angle_2 = 50.0
    euler_angle_3 = -132.0
    block = '3'
  [../]
  [./beta_tensor_3]
    type = ComputeBetaTensor
    block = '3'
    euler_angle_1 = 1.0
    euler_angle_2 = 50.0
    euler_angle_3 = -132.0
  [../]
  [./delta_beta_tensor_3]
    type = ComputeDeltaBetaTensor
    block = '3'
  [../]


  [./elasticity_tensor_4]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 100.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
    block = '4'
  [../]
  [./strain_4]
    type = ComputeSmallStrain
    block = '4'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_4]
    type = ComputeLinearElasticStress
    block = '4'
  [../]
  [./photoelastic_tensor_4]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 100.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
    block = '4'
  [../]
  [./beta_tensor_4]
    type = ComputeBetaTensor
    block = '4'
    euler_angle_1 = 100.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]
  [./delta_beta_tensor_4]
    type = ComputeDeltaBetaTensor
    block = '4'
  [../]


  [./elasticity_tensor_5]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 200.0
    euler_angle_2 = 33.0
    euler_angle_3 = 5.0
    block = '5'
  [../]
  [./strain_5]
    type = ComputeSmallStrain
    block = '5'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_5]
    type = ComputeLinearElasticStress
    block = '5'
  [../]
  [./photoelastic_tensor_5]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 200.0
    euler_angle_2 = 33.0
    euler_angle_3 = 5.0
    block = '5'
  [../]
  [./beta_tensor_5]
    type = ComputeBetaTensor
    block = '5'
    euler_angle_1 = 200.0
    euler_angle_2 = 33.0
    euler_angle_3 = 5.0
  [../]
  [./delta_beta_tensor_5]
    type = ComputeDeltaBetaTensor
    block = '5'
  [../]


  [./elasticity_tensor_6]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 110.0
    euler_angle_2 = 0.0
    euler_angle_3 = 100.0
    block = '6'
  [../]
  [./strain_6]
    type = ComputeSmallStrain
    block = '6'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_6]
    type = ComputeLinearElasticStress
    block = '6'
  [../]
  [./photoelastic_tensor_6]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 110.0
    euler_angle_2 = 0.0
    euler_angle_3 = 100.0
    block = '6'
  [../]
  [./beta_tensor_6]
    type = ComputeBetaTensor
    block = '6'
    euler_angle_1 = 110.0
    euler_angle_2 = 0.0
    euler_angle_3 = 100.0
  [../]
  [./delta_beta_tensor_6]
    type = ComputeDeltaBetaTensor
    block = '6'
  [../]


  [./elasticity_tensor_7]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 33.0
    euler_angle_2 = -37.0
    euler_angle_3 = 25.0
    block = '7'
  [../]
  [./strain_7]
    type = ComputeSmallStrain
    block = '7'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_7]
    type = ComputeLinearElasticStress
    block = '7'
  [../]
  [./photoelastic_tensor_7]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 33.0
    euler_angle_2 = -37.0
    euler_angle_3 = 25.0
    block = '7'
  [../]
  [./beta_tensor_7]
    type = ComputeBetaTensor
    block = '7'
    euler_angle_1 = 33.0
    euler_angle_2 = -37.0
    euler_angle_3 = 25.0
  [../]
  [./delta_beta_tensor_7]
    type = ComputeDeltaBetaTensor
    block = '7'
  [../]


  [./elasticity_tensor_8]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 10.0
    euler_angle_2 = 137.0
    euler_angle_3 = -228.0
    block = '8'
  [../]
  [./strain_8]
    type = ComputeSmallStrain
    block = '8'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_8]
    type = ComputeLinearElasticStress
    block = '8'
  [../]
  [./photoelastic_tensor_8]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 10.0
    euler_angle_2 = 137.0
    euler_angle_3 = -228.0
    block = '8'
  [../]
  [./beta_tensor_8]
    type = ComputeBetaTensor
    block = '8'
    euler_angle_1 = 10.0
    euler_angle_2 = 137.0
    euler_angle_3 = -228.0
  [../]
  [./delta_beta_tensor_8]
    type = ComputeDeltaBetaTensor
    block = '8'
  [../]

  [./elasticity_tensor_9]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -250.0
    euler_angle_2 = 33.0
    euler_angle_3 = -5.0
    block = '9'
  [../]
  [./strain_9]
    type = ComputeSmallStrain
    block = '9'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_9]
    type = ComputeLinearElasticStress
    block = '9'
  [../]
  [./photoelastic_tensor_9]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = -250.0
    euler_angle_2 = 33.0
    euler_angle_3 = -5.0
    block = '9'
  [../]
  [./beta_tensor_9]
    type = ComputeBetaTensor
    block = '9'
    euler_angle_1 = -250.0
    euler_angle_2 = 33.0
    euler_angle_3 = -5.0
  [../]
  [./delta_beta_tensor_9]
    type = ComputeDeltaBetaTensor
    block = '9'
  [../]



  [./elasticity_tensor_10]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -25.0
    euler_angle_2 = 63.0
    euler_angle_3 = 125.0
    block = '10'
  [../]
  [./strain_10]
    type = ComputeSmallStrain
    block = '10'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_10]
    type = ComputeLinearElasticStress
    block = '10'
  [../]
  [./photoelastic_tensor_10]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = -25.0
    euler_angle_2 = 63.0
    euler_angle_3 = 125.0
    block = '10'
  [../]
  [./beta_tensor_10]
    type = ComputeBetaTensor
    block = '10'
    euler_angle_1 = -25.0
    euler_angle_2 = 63.0
    euler_angle_3 = 125.0
  [../]
  [./delta_beta_tensor_10]
    type = ComputeDeltaBetaTensor
    block = '10'
  [../]


  [./elasticity_tensor_11]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -25.0
    euler_angle_2 = 63.0
    euler_angle_3 = 125.0
    block = '11'
  [../]
  [./strain_11]
    type = ComputeSmallStrain
    block = '11'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_11]
    type = ComputeLinearElasticStress
    block = '11'
  [../]
  [./photoelastic_tensor_11]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = -25.0
    euler_angle_2 = 63.0
    euler_angle_3 = 125.0
    block = '11'
  [../]
  [./beta_tensor_11]
    type = ComputeBetaTensor
    block = '11'
    euler_angle_1 = -25.0
    euler_angle_2 = 63.0
    euler_angle_3 = 125.0
  [../]
  [./delta_beta_tensor_11]
    type = ComputeDeltaBetaTensor
    block = '11'
  [../]


  [./elasticity_tensor_12]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -125.0
    euler_angle_2 = 23.0
    euler_angle_3 = 5.7
    block = '12'
  [../]
  [./strain_12]
    type = ComputeSmallStrain
    block = '12'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_12]
    type = ComputeLinearElasticStress
    block = '12'
  [../]
  [./photoelastic_tensor_12]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = -125.0
    euler_angle_2 = 23.0
    euler_angle_3 = 5.7
    block = '12'
  [../]
  [./beta_tensor_12]
    type = ComputeBetaTensor
    block = '12'
    euler_angle_1 = -125.0
    euler_angle_2 = 23.0
    euler_angle_3 = 5.7
  [../]
  [./delta_beta_tensor_12]
    type = ComputeDeltaBetaTensor
    block = '12'
  [../]


  [./elasticity_tensor_13]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 25.0
    euler_angle_2 = -310.0
    euler_angle_3 = 57.9
    block = '13'
  [../]
  [./strain_13]
    type = ComputeSmallStrain
    block = '13'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_13]
    type = ComputeLinearElasticStress
    block = '13'
  [../]
  [./photoelastic_tensor_13]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 25.0
    euler_angle_2 = -310.0
    euler_angle_3 = 57.9
    block = '13'
  [../]
  [./beta_tensor_13]
    type = ComputeBetaTensor
    block = '13'
    euler_angle_1 = 25.0
    euler_angle_2 = -310.0
    euler_angle_3 = 57.9
  [../]
  [./delta_beta_tensor_13]
    type = ComputeDeltaBetaTensor
    block = '13'
  [../]


  [./elasticity_tensor_14]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 77.0
    euler_angle_2 = 5.0
    euler_angle_3 = -32.0
    block = '14'
  [../]
  [./strain_14]
    type = ComputeSmallStrain
    block = '14'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_14]
    type = ComputeLinearElasticStress
    block = '14'
  [../]
  [./photoelastic_tensor_14]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 77.0
    euler_angle_2 = 5.0
    euler_angle_3 = -32.0
    block = '14'
  [../]
  [./beta_tensor_14]
    type = ComputeBetaTensor
    block = '14'
    euler_angle_1 = 77.0
    euler_angle_2 = 5.0
    euler_angle_3 = -32.0
  [../]
  [./delta_beta_tensor_14]
    type = ComputeDeltaBetaTensor
    block = '14'
  [../]


  [./elasticity_tensor_15]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -250.0
    euler_angle_2 = 15.0
    euler_angle_3 = -77.9
    block = '15'
  [../]
  [./strain_15]
    type = ComputeSmallStrain
    block = '15'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_15]
    type = ComputeLinearElasticStress
    block = '15'
  [../]
  [./photoelastic_tensor_15]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = -250.0
    euler_angle_2 = 15.0
    euler_angle_3 = -77.9
    block = '15'
  [../]
  [./beta_tensor_15]
    type = ComputeBetaTensor
    block = '15'
    euler_angle_1 = -250.0
    euler_angle_2 = 15.0
    euler_angle_3 = -77.9
  [../]
  [./delta_beta_tensor_15]
    type = ComputeDeltaBetaTensor
    block = '15'
  [../]


  [./elasticity_tensor_16]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 31.2
    euler_angle_2 = -31.0
    euler_angle_3 = -31.9
    block = '16'
  [../]
  [./strain_16]
    type = ComputeSmallStrain
    block = '16'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_16]
    type = ComputeLinearElasticStress
    block = '16'
  [../]
  [./photoelastic_tensor_16]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 31.2
    euler_angle_2 = -31.0
    euler_angle_3 = -31.9
    block = '16'
  [../]
  [./beta_tensor_16]
    type = ComputeBetaTensor
    block = '16'
    euler_angle_1 = 31.2
    euler_angle_2 = -31.0
    euler_angle_3 = -31.9
  [../]
  [./delta_beta_tensor_16]
    type = ComputeDeltaBetaTensor
    block = '16'
  [../]

  [./elasticity_tensor_17]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -185.0
    euler_angle_2 = -110.0
    euler_angle_3 = 7.9
    block = '17'
  [../]
  [./strain_17]
    type = ComputeSmallStrain
    block = '17'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_17]
    type = ComputeLinearElasticStress
    block = '17'
  [../]
  [./photoelastic_tensor_17]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = -185.0
    euler_angle_2 = -110.0
    euler_angle_3 = 7.9
    block = '17'
  [../]
  [./beta_tensor_17]
    type = ComputeBetaTensor
    block = '17'
    euler_angle_1 = -185.0
    euler_angle_2 = -110.0
    euler_angle_3 = 7.9
  [../]
  [./delta_beta_tensor_17]
    type = ComputeDeltaBetaTensor
    block = '17'
  [../]

  [./elasticity_tensor_18]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -55.0
    euler_angle_2 = -51.0
    euler_angle_3 = -51.9
    block = '18'
  [../]
  [./strain_18]
    type = ComputeSmallStrain
    block = '18'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_18]
    type = ComputeLinearElasticStress
    block = '18'
  [../]
  [./photoelastic_tensor_18]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = -55.0
    euler_angle_2 = -51.0
    euler_angle_3 = -51.9
    block = '18'
  [../]
  [./beta_tensor_18]
    type = ComputeBetaTensor
    block = '18'
    euler_angle_1 = -55.0
    euler_angle_2 = -51.0
    euler_angle_3 = -51.9
  [../]
  [./delta_beta_tensor_18]
    type = ComputeDeltaBetaTensor
    block = '18'
  [../]

  [./elasticity_tensor_19]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 0.0
    euler_angle_2 = -37.0
    euler_angle_3 = 128.0
    block = '19'
  [../]
  [./strain_19]
    type = ComputeSmallStrain
    block = '19'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_19]
    type = ComputeLinearElasticStress
    block = '19'
  [../]
  [./photoelastic_tensor_19]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 0.0
    euler_angle_2 = -37.0
    euler_angle_3 = 128.0
    block = '19'
  [../]
  [./beta_tensor_19]
    type = ComputeBetaTensor
    block = '19'
    euler_angle_1 = 0.0
    euler_angle_2 = -37.0
    euler_angle_3 = 128.0
  [../]
  [./delta_beta_tensor_19]
    type = ComputeDeltaBetaTensor
    block = '19'
  [../]


  [./elasticity_tensor_20]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -55.0
    euler_angle_2 = -11.0
    euler_angle_3 = 11.9
    block = '20'
  [../]
  [./strain_20]
    type = ComputeSmallStrain
    block = '20'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_20]
    type = ComputeLinearElasticStress
    block = '20'
  [../]
  [./photoelastic_tensor_20]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = -55.0
    euler_angle_2 = -11.0
    euler_angle_3 = 11.9
    block = '20'
  [../]
  [./beta_tensor_20]
    type = ComputeBetaTensor
    block = '20'
    euler_angle_1 = -55.0
    euler_angle_2 = -11.0
    euler_angle_3 = 11.9
  [../]
  [./delta_beta_tensor_20]
    type = ComputeDeltaBetaTensor
    block = '20'
  [../]


  [./elasticity_tensor_21]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -260.0
    euler_angle_2 = -105.0
    euler_angle_3 = 76.9
    block = '21'
  [../]
  [./strain_21]
    type = ComputeSmallStrain
    block = '21'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_21]
    type = ComputeLinearElasticStress
    block = '21'
  [../]
  [./photoelastic_tensor_21]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = -260.0
    euler_angle_2 = -105.0
    euler_angle_3 = 76.9
    block = '21'
  [../]
  [./beta_tensor_21]
    type = ComputeBetaTensor
    block = '21'
    euler_angle_1 = -260.0
    euler_angle_2 = -105.0
    euler_angle_3 = 76.9
  [../]
  [./delta_beta_tensor_21]
    type = ComputeDeltaBetaTensor
    block = '21'
  [../]


  [./elasticity_tensor_22]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -81.0
    euler_angle_2 = -81.0
    euler_angle_3 = 81.9
    block = '22'
  [../]
  [./strain_22]
    type = ComputeSmallStrain
    block = '22'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_22]
    type = ComputeLinearElasticStress
    block = '22'
  [../]
  [./photoelastic_tensor_22]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = -81.0
    euler_angle_2 = -81.0
    euler_angle_3 = 81.9
    block = '22'
  [../]
  [./beta_tensor_22]
    type = ComputeBetaTensor
    block = '22'
    euler_angle_1 = -81.0
    euler_angle_2 = -81.0
    euler_angle_3 = 81.9
  [../]
  [./delta_beta_tensor_22]
    type = ComputeDeltaBetaTensor
    block = '22'
  [../]


  [./elasticity_tensor_23]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -4.0
    euler_angle_2 = -41.0
    euler_angle_3 = 101.9
    block = '23'
  [../]
  [./strain_23]
    type = ComputeSmallStrain
    block = '23'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_23]
    type = ComputeLinearElasticStress
    block = '23'
  [../]
  [./photoelastic_tensor_23]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = -4.0
    euler_angle_2 = -41.0
    euler_angle_3 = 101.9
    block = '23'
  [../]
  [./beta_tensor_23]
    type = ComputeBetaTensor
    block = '23'
    euler_angle_1 = -4.0
    euler_angle_2 = -41.0
    euler_angle_3 = 101.9
  [../]
  [./delta_beta_tensor_23]
    type = ComputeDeltaBetaTensor
    block = '23'
  [../]



  [./elasticity_tensor_24]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -125.0
    euler_angle_2 = -101.0
    euler_angle_3 = 31.9
    block = '24'
  [../]
  [./strain_24]
    type = ComputeSmallStrain
    block = '24'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_24]
    type = ComputeLinearElasticStress
    block = '24'
  [../]
  [./photoelastic_tensor_24]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = -125.0
    euler_angle_2 = -101.0
    euler_angle_3 = 31.9
    block = '24'
  [../]
  [./beta_tensor_24]
    type = ComputeBetaTensor
    block = '24'
    euler_angle_1 = -125.0
    euler_angle_2 = -101.0
    euler_angle_3 = 31.9
  [../]
  [./delta_beta_tensor_24]
    type = ComputeDeltaBetaTensor
    block = '24'
  [../]


  [./elasticity_tensor_25]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -205.0
    euler_angle_2 = -201.0
    euler_angle_3 = 81.5
    block = '25'
  [../]
  [./strain_25]
    type = ComputeSmallStrain
    block = '25'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_25]
    type = ComputeLinearElasticStress
    block = '25'
  [../]
  [./photoelastic_tensor_25]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = -205.0
    euler_angle_2 = -201.0
    euler_angle_3 = 81.5
    block = '25'
  [../]
  [./beta_tensor_25]
    type = ComputeBetaTensor
    block = '25'
    euler_angle_1 = -205.0
    euler_angle_2 = -201.0
    euler_angle_3 = 81.5
  [../]
  [./delta_beta_tensor_25]
    type = ComputeDeltaBetaTensor
    block = '25'
  [../]


  [./elasticity_tensor_26]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -7.0
    euler_angle_2 = -17.0
    euler_angle_3 = 77.9
    block = '26'
  [../]
  [./strain_26]
    type = ComputeSmallStrain
    block = '26'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_26]
    type = ComputeLinearElasticStress
    block = '26'
  [../]
  [./photoelastic_tensor_26]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = -7.0
    euler_angle_2 = -17.0
    euler_angle_3 = 77.9
    block = '26'
  [../]
  [./beta_tensor_26]
    type = ComputeBetaTensor
    block = '26'
    euler_angle_1 = -7.0
    euler_angle_2 = -17.0
    euler_angle_3 = 77.9
  [../]
  [./delta_beta_tensor_26]
    type = ComputeDeltaBetaTensor
    block = '26'
  [../]


  [./elasticity_tensor_27]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 12.0
    euler_angle_2 = 12.0
    euler_angle_3 = 44.9
    block = '27'
  [../]
  [./strain_27]
    type = ComputeSmallStrain
    block = '27'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_27]
    type = ComputeLinearElasticStress
    block = '27'
  [../]
  [./photoelastic_tensor_27]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 12.0
    euler_angle_2 = 12.0
    euler_angle_3 = 44.9
    block = '27'
  [../]
  [./beta_tensor_27]
    type = ComputeBetaTensor
    block = '27'
    euler_angle_1 = 12.0
    euler_angle_2 = 12.0
    euler_angle_3 = 44.9
  [../]
  [./delta_beta_tensor_27]
    type = ComputeDeltaBetaTensor
    block = '27'
  [../]


  [./elasticity_tensor_28]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 12.0
    euler_angle_2 = -22.0
    euler_angle_3 = 77.9
    block = '28'
  [../]
  [./strain_28]
    type = ComputeSmallStrain
    block = '28'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_28]
    type = ComputeLinearElasticStress
    block = '28'
  [../]
  [./photoelastic_tensor_28]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 12.0
    euler_angle_2 = -22.0
    euler_angle_3 = 77.9
    block = '28'
  [../]
  [./beta_tensor_28]
    type = ComputeBetaTensor
    block = '28'
    euler_angle_1 = 12.0
    euler_angle_2 = -22.0
    euler_angle_3 = 77.9
  [../]
  [./delta_beta_tensor_28]
    type = ComputeDeltaBetaTensor
    block = '28'
  [../]


  [./elasticity_tensor_29]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 100.0
    euler_angle_2 = -33.0
    euler_angle_3 = 5.0
    block = '29'
  [../]
  [./strain_29]
    type = ComputeSmallStrain
    block = '29'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_29]
    type = ComputeLinearElasticStress
    block = '29'
  [../]
  [./photoelastic_tensor_29]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 100.0
    euler_angle_2 = -33.0
    euler_angle_3 = 5.0
    block = '29'
  [../]
  [./beta_tensor_29]
    type = ComputeBetaTensor
    block = '29'
    euler_angle_1 = 100.0
    euler_angle_2 = -33.0
    euler_angle_3 = 5.0
  [../]
  [./delta_beta_tensor_29]
    type = ComputeDeltaBetaTensor
    block = '29'
  [../]


  [./elasticity_tensor_30]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 0.0
    euler_angle_2 = 200.0
    euler_angle_3 = 30.0
    block = '30'
  [../]
  [./strain_30]
    type = ComputeSmallStrain
    block = '30'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_30]
    type = ComputeLinearElasticStress
    block = '30'
  [../]
  [./photoelastic_tensor_30]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 0.0
    euler_angle_2 = 200.0
    euler_angle_3 = 30.0
    block = '30'
  [../]
  [./beta_tensor_30]
    type = ComputeBetaTensor
    block = '30'
    euler_angle_1 = 0.0
    euler_angle_2 = 200.0
    euler_angle_3 = 30.0
  [../]
  [./delta_beta_tensor_30]
    type = ComputeDeltaBetaTensor
    block = '30'
  [../]


  [./elasticity_tensor_31]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 63.0
    euler_angle_2 = 37.0
    euler_angle_3 = -15.5
    block = '31'
  [../]
  [./strain_31]
    type = ComputeSmallStrain
    block = '31'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_31]
    type = ComputeLinearElasticStress
    block = '31'
  [../]
  [./photoelastic_tensor_31]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 63.0
    euler_angle_2 = 37.0
    euler_angle_3 = -15.5
    block = '31'
  [../]
  [./beta_tensor_31]
    type = ComputeBetaTensor
    block = '31'
    euler_angle_1 = 63.0
    euler_angle_2 = 37.0
    euler_angle_3 = -15.5
  [../]
  [./delta_beta_tensor_31]
    type = ComputeDeltaBetaTensor
    block = '31'
  [../]



  [./elasticity_tensor_32]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -18.5
    euler_angle_2 = -44.0
    euler_angle_3 = 176.5
    block = '32'
  [../]
  [./strain_32]
    type = ComputeSmallStrain
    block = '32'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_32]
    type = ComputeLinearElasticStress
    block = '32'
  [../]
  [./photoelastic_tensor_32]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = -18.5
    euler_angle_2 = -44.0
    euler_angle_3 = 176.5
    block = '32'
  [../]
  [./beta_tensor_32]
    type = ComputeBetaTensor
    block = '32'
    euler_angle_1 = -18.5
    euler_angle_2 = -44.0
    euler_angle_3 = 176.5
  [../]
  [./delta_beta_tensor_32]
    type = ComputeDeltaBetaTensor
    block = '32'
  [../]
[]

[Kernels]
  #Elastic problem
  [./TensorMechanics]
  #This is an action block
  [../]
[]


[BCs]
  [./center_disp_z_top]
    type = DirichletBC
    variable = 'disp_z'
    value = -0.5
    boundary = 'side1'
  [../]

  [./center_disp_z_bottom]
    type = DirichletBC
    variable = 'disp_z'
    value = 0.0
    boundary = 'side2'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type  -pc_hypre_type'
    petsc_options_value = '    250              1e-10      1e-8      1e-6      hypre       boomeramg '
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    execute_on = 'timestep_end'
    file_base = out_32grain_r100_struct_load_50top
    elemental_as_nodal = true
  [../]
[]
