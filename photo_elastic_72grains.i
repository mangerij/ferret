
[Mesh]
  file = 3D_HCP_72.e
[]

[MeshModifiers]
  [./add_side_sets]
    type = SideSetsFromNormals
    normals = '1  0  0
               -1  0  0
               0  0  -1'
    new_boundary = 'side1 side2 back1'
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
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
    block = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15  16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70'
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

  [./elasticity_tensor_5]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -56.0
    euler_angle_2 = -56.0
    euler_angle_3 = 10.0
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

  [./elasticity_tensor_9]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -100.0
    euler_angle_2 = -195.0
    euler_angle_3 = -18.0
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


  [./elasticity_tensor_33]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 105.0
    euler_angle_2 = -63.0
    euler_angle_3 = 115.5
    block = '33'
  [../]
  [./strain_33]
    type = ComputeSmallStrain
    block = '33'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_33]
    type = ComputeLinearElasticStress
    block = '33'
  [../]


  [./elasticity_tensor_34]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
    block = '34'
  [../]
  [./strain_34]
    type = ComputeSmallStrain
    block = '34'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_34]
    type = ComputeLinearElasticStress
    block = '34'
  [../]

  [./elasticity_tensor_35]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 70.0
    euler_angle_2 = 70.0
    euler_angle_3 = 70.0
    block = '35'
  [../]
  [./strain_35]
    type = ComputeSmallStrain
    block = '35'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_35]
    type = ComputeLinearElasticStress
    block = '35'
  [../]

  [./elasticity_tensor_36]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 14.0
    euler_angle_2 = -37.0
    euler_angle_3 = 8.0
    block = '36'
  [../]
  [./strain_36]
    type = ComputeSmallStrain
    block = '36'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_36]
    type = ComputeLinearElasticStress
    block = '36'
  [../]

  [./elasticity_tensor_37]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 44.0
    euler_angle_2 = 13.0
    euler_angle_3 = 38.0
    block = '37'
  [../]
  [./strain_37]
    type = ComputeSmallStrain
    block = '37'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_37]
    type = ComputeLinearElasticStress
    block = '37'
  [../]

  [./elasticity_tensor_38]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 124.0
    euler_angle_2 = 0.0
    euler_angle_3 = 86.0
    block = '38'
  [../]
  [./strain_38]
    type = ComputeSmallStrain
    block = '38'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_38]
    type = ComputeLinearElasticStress
    block = '38'
  [../]

  [./elasticity_tensor_39]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 77.0
    euler_angle_2 = -17.0
    euler_angle_3 = -16.0
    block = '39'
  [../]
  [./strain_39]
    type = ComputeSmallStrain
    block = '39'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_39]
    type = ComputeLinearElasticStress
    block = '39'
  [../]

  [./elasticity_tensor_40]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 24.0
    euler_angle_2 = 24.0
    euler_angle_3 = 24.0
    block = '40'
  [../]
  [./strain_40]
    type = ComputeSmallStrain
    block = '40'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_40]
    type = ComputeLinearElasticStress
    block = '40'
  [../]

  [./elasticity_tensor_41]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -240.0
    euler_angle_2 = -240.0
    euler_angle_3 = 240.0
    block = '41'
  [../]
  [./strain_41]
    type = ComputeSmallStrain
    block = '41'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_41]
    type = ComputeLinearElasticStress
    block = '41'
  [../]

  [./elasticity_tensor_42]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 1.2
    euler_angle_2 = 8.6
    euler_angle_3 = -140.0
    block = '42'
  [../]
  [./strain_42]
    type = ComputeSmallStrain
    block = '42'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_42]
    type = ComputeLinearElasticStress
    block = '42'
  [../]

  [./elasticity_tensor_43]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 12.0
    euler_angle_2 = -86.7
    euler_angle_3 = 5.0
    block = '43'
  [../]
  [./strain_43]
    type = ComputeSmallStrain
    block = '43'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_43]
    type = ComputeLinearElasticStress
    block = '43'
  [../]

  [./elasticity_tensor_44]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 12.2
    euler_angle_2 = 81.6
    euler_angle_3 = 70.0
    block = '44'
  [../]
  [./strain_44]
    type = ComputeSmallStrain
    block = '44'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_44]
    type = ComputeLinearElasticStress
    block = '44'
  [../]

  [./elasticity_tensor_45]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -120.2
    euler_angle_2 = 5.0
    euler_angle_3 = 5.0
    block = '45'
  [../]
  [./strain_45]
    type = ComputeSmallStrain
    block = '45'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_45]
    type = ComputeLinearElasticStress
    block = '45'
  [../]

  [./elasticity_tensor_46]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 10.8
    euler_angle_2 = 55.0
    euler_angle_3 = -75.2
    block = '46'
  [../]
  [./strain_46]
    type = ComputeSmallStrain
    block = '46'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_46]
    type = ComputeLinearElasticStress
    block = '46'
  [../]

  [./elasticity_tensor_48]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 18.0
    euler_angle_2 = 185.0
    euler_angle_3 = 66.0
    block = '48'
  [../]
  [./strain_48]
    type = ComputeSmallStrain
    block = '48'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_48]
    type = ComputeLinearElasticStress
    block = '48'
  [../]

  [./elasticity_tensor_49]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 66.0
    euler_angle_2 = 66.0
    euler_angle_3 = 66.0
    block = '49'
  [../]
  [./strain_49]
    type = ComputeSmallStrain
    block = '49'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_49]
    type = ComputeLinearElasticStress
    block = '49'
  [../]

  [./elasticity_tensor_50]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 100.0
    euler_angle_2 = -45.0
    euler_angle_3 = 0.0
    block = '50'
  [../]
  [./strain_50]
    type = ComputeSmallStrain
    block = '50'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_50]
    type = ComputeLinearElasticStress
    block = '50'
  [../]

  [./elasticity_tensor_51]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -25.0
    euler_angle_2 = 147.0
    euler_angle_3 = 147.0
    block = '51'
  [../]
  [./strain_51]
    type = ComputeSmallStrain
    block = '51'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_51]
    type = ComputeLinearElasticStress
    block = '51'
  [../]

  [./elasticity_tensor_52]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -17.0
    euler_angle_2 = 260.0
    euler_angle_3 = 7.5
    block = '52'
  [../]
  [./strain_52]
    type = ComputeSmallStrain
    block = '52'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_52]
    type = ComputeLinearElasticStress
    block = '52'
  [../]

  [./elasticity_tensor_53]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -17.0
    euler_angle_2 = -45.0
    euler_angle_3 = 42.0
    block = '53'
  [../]
  [./strain_53]
    type = ComputeSmallStrain
    block = '53'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_53]
    type = ComputeLinearElasticStress
    block = '53'
  [../]

  [./elasticity_tensor_54]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 15.0
    euler_angle_2 = -99.0
    euler_angle_3 = 99.0
    block = '54'
  [../]
  [./strain_54]
    type = ComputeSmallStrain
    block = '54'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_54]
    type = ComputeLinearElasticStress
    block = '54'
  [../]

  [./elasticity_tensor_55]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -99.0
    euler_angle_2 = 37.0
    euler_angle_3 = 9.8
    block = '55'
  [../]
  [./strain_55]
    type = ComputeSmallStrain
    block = '55'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_55]
    type = ComputeLinearElasticStress
    block = '55'
  [../]

  [./elasticity_tensor_56]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 10.0
    euler_angle_2 = 10.0
    euler_angle_3 = 198.73
    block = '56'
  [../]
  [./strain_56]
    type = ComputeSmallStrain
    block = '56'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_56]
    type = ComputeLinearElasticStress
    block = '56'
  [../]

  [./elasticity_tensor_57]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -100.0
    euler_angle_2 = 32.0
    euler_angle_3 = 32.73
    block = '57'
  [../]
  [./strain_57]
    type = ComputeSmallStrain
    block = '57'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_57]
    type = ComputeLinearElasticStress
    block = '57'
  [../]

  [./elasticity_tensor_58]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 5.0
    euler_angle_2 = 35.0
    euler_angle_3 = 185.45
    block = '58'
  [../]
  [./strain_58]
    type = ComputeSmallStrain
    block = '58'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_58]
    type = ComputeLinearElasticStress
    block = '58'
  [../]

  [./elasticity_tensor_59]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 67.0
    euler_angle_2 = -23.0
    euler_angle_3 = -68.0
    block = '59'
  [../]
  [./strain_59]
    type = ComputeSmallStrain
    block = '59'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_59]
    type = ComputeLinearElasticStress
    block = '59'
  [../]

  [./elasticity_tensor_60]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 18.0
    euler_angle_2 = 8.0
    euler_angle_3 = -105.0
    block = '60'
  [../]
  [./strain_60]
    type = ComputeSmallStrain
    block = '60'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_60]
    type = ComputeLinearElasticStress
    block = '60'
  [../]

  [./elasticity_tensor_61]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 0.0
    euler_angle_2 = 45.0
    euler_angle_3 = -16.0
    block = '61'
  [../]
  [./strain_61]
    type = ComputeSmallStrain
    block = '61'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_61]
    type = ComputeLinearElasticStress
    block = '61'
  [../]

  [./elasticity_tensor_62]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 195.0
    euler_angle_2 = 195.0
    euler_angle_3 = -16.0
    block = '62'
  [../]
  [./strain_62]
    type = ComputeSmallStrain
    block = '62'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_62]
    type = ComputeLinearElasticStress
    block = '62'
  [../]

  [./elasticity_tensor_63]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 45.0
    euler_angle_2 = 45.0
    euler_angle_3 = -32.0
    block = '63'
  [../]
  [./strain_63]
    type = ComputeSmallStrain
    block = '63'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_63]
    type = ComputeLinearElasticStress
    block = '63'
  [../]

  [./elasticity_tensor_64]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 116.0
    euler_angle_2 = 116.0
    euler_angle_3 = -2.5
    block = '64'
  [../]
  [./strain_64]
    type = ComputeSmallStrain
    block = '64'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_64]
    type = ComputeLinearElasticStress
    block = '64'
  [../]

  [./elasticity_tensor_65]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 20.0
    euler_angle_2 = 26.0
    euler_angle_3 = -26.5
    block = '65'
  [../]
  [./strain_65]
    type = ComputeSmallStrain
    block = '65'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_65]
    type = ComputeLinearElasticStress
    block = '65'
  [../]

  [./elasticity_tensor_66]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 16.0
    euler_angle_2 = -37.0
    euler_angle_3 = 138.5
    block = '66'
  [../]
  [./strain_66]
    type = ComputeSmallStrain
    block = '66'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_66]
    type = ComputeLinearElasticStress
    block = '66'
  [../]

  [./elasticity_tensor_67]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 165.0
    euler_angle_2 = 100.0
    euler_angle_3 = 177.78
    block = '67'
  [../]
  [./strain_67]
    type = ComputeSmallStrain
    block = '67'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_67]
    type = ComputeLinearElasticStress
    block = '67'
  [../]

  [./elasticity_tensor_68]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 0.0
    euler_angle_2 = 22.0
    euler_angle_3 = 77.3
    block = '68'
  [../]
  [./strain_68]
    type = ComputeSmallStrain
    block = '68'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_68]
    type = ComputeLinearElasticStress
    block = '68'
  [../]

  [./elasticity_tensor_69]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -33.0
    euler_angle_2 = 28.0
    euler_angle_3 = 39.3
    block = '69'
  [../]
  [./strain_69]
    type = ComputeSmallStrain
    block = '69'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_69]
    type = ComputeLinearElasticStress
    block = '69'
  [../]

  [./elasticity_tensor_70]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -7.0
    euler_angle_2 = 6.0
    euler_angle_3 = 144.0
    block = '70'
  [../]
  [./strain_70]
    type = ComputeSmallStrain
    block = '70'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_70]
    type = ComputeLinearElasticStress
    block = '70'
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
    variable = 'disp_x'
    value = -0.05
    boundary = 'side1'
  [../]

  [./center_disp_z_bottom]
    type = DirichletBC
    variable = 'disp_x'
    value = 0.05
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
    file_base = out_144grain_structure
    elemental_as_nodal = true
  [../]
[]
