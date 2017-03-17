[Mesh]
  file = 3D_HCP_grainy_nospheres.e
[]

[MeshModifiers]
  [./add_right]
    type = SideSetsFromNormals
    normals = '1  0  0'
    new_boundary = 'right'
  [../]
  [./add_left]
    type = SideSetsFromNormals
    normals = '-1  0  0'
    new_boundary = 'left'
  [../]
[]

[GlobalParams]
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  displacements = 'disp_x disp_y disp_z'
  use_displaced_mesh = true
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
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 0
    variable = strain_xx_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 1
    variable = strain_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e13]
    type = RankTwoAux
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 2
    variable = strain_xz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 1
    variable = strain_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e23]
    type = RankTwoAux
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 2
    variable = strain_yz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e33]
    type = RankTwoAux
    rank_two_tensor = total_strain
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

[Materials]
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
    block = '1'
    euler_angle_1 = 35.0
    euler_angle_2 = 44.0
    euler_angle_3 = 97.0
  [../]
  [./elasticity_tensor_2]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
    block = '2'
    euler_angle_1 = 95.0
    euler_angle_2 = 123.0
    euler_angle_3 = 10.0
  [../]
  [./elasticity_tensor_85]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
    block = '85'
    euler_angle_1 = 1.0
    euler_angle_2 = 33.0
    euler_angle_3 = 35.0
  [../]
  [./elasticity_tensor_all]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
    block = '3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    block = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113'
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
    block = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113'
  [../]
[]

[Kernels]
  [./TensorMechanics]
  [../]
[]



[BCs]
  [./disp_x_right]
    type = DirichletBC
    variable = disp_x
    boundary = 'right'
    value = -0.01
  [../]
  [./disp_x_left]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0.01
  [../]
[]


[Postprocessors]
    [./Felastic]
      type = ElasticEnergy
      block = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113'
      execute_on = 'timestep_end'
    [../]
    [./grain_vol_1]
      type = VolumePostprocessor
      block = '1'
    [../]
    [./grain_vol_2]
      type = VolumePostprocessor
      block = '2'
    [../]
    [./grain_vol_3]
      type = VolumePostprocessor
      block = '3'
    [../]
    [./grain_vol_4]
      type = VolumePostprocessor
      block = '4'
    [../]
    [./grain_vol_5]
      type = VolumePostprocessor
      block = '5'
    [../]
    [./grain_vol_6]
      type = VolumePostprocessor
      block = '6'
    [../]
    [./grain_vol_7]
      type = VolumePostprocessor
      block = '7'
    [../]
    [./grain_vol_8]
      type = VolumePostprocessor
      block = '8'
    [../]
    [./grain_vol_9]
      type = VolumePostprocessor
      block = '9'
    [../]
    [./grain_vol_10]
      type = VolumePostprocessor
      block = '10'
    [../]
    [./grain_vol_11]
      type = VolumePostprocessor
      block = '11'
    [../]
    [./grain_vol_12]
      type = VolumePostprocessor
      block = '12'
    [../]
    [./grain_vol_13]
      type = VolumePostprocessor
      block = '13'
    [../]
    [./grain_vol_14]
      type = VolumePostprocessor
      block = '14'
    [../]
    [./grain_vol_15]
      type = VolumePostprocessor
      block = '15'
    [../]
    [./grain_vol_16]
      type = VolumePostprocessor
      block = '16'
    [../]
    [./grain_vol_17]
      type = VolumePostprocessor
      block = '17'
    [../]
    [./grain_vol_18]
      type = VolumePostprocessor
      block = '18'
    [../]
    [./grain_vol_19]
      type = VolumePostprocessor
      block = '19'
    [../]
    [./grain_vol_20]
      type = VolumePostprocessor
      block = '20'
    [../]
    [./grain_vol_21]
      type = VolumePostprocessor
      block = '21'
    [../]
    [./grain_vol_22]
      type = VolumePostprocessor
      block = '22'
    [../]
    [./grain_vol_23]
      type = VolumePostprocessor
      block = '23'
    [../]
    [./grain_vol_24]
      type = VolumePostprocessor
      block = '24'
    [../]
    [./grain_vol_25]
      type = VolumePostprocessor
      block = '25'
    [../]
    [./grain_vol_26]
      type = VolumePostprocessor
      block = '26'
    [../]
    [./grain_vol_27]
      type = VolumePostprocessor
      block = '27'
    [../]
    [./grain_vol_28]
      type = VolumePostprocessor
      block = '28'
    [../]
    [./grain_vol_29]
      type = VolumePostprocessor
      block = '29'
    [../]
    [./grain_vol_30]
      type = VolumePostprocessor
      block = '30'
    [../]
    [./grain_vol_31]
      type = VolumePostprocessor
      block = '31'
    [../]
    [./grain_vol_32]
      type = VolumePostprocessor
      block = '32'
    [../]
    [./grain_vol_33]
      type = VolumePostprocessor
      block = '33'
    [../]
    [./grain_vol_34]
      type = VolumePostprocessor
      block = '34'
    [../]
    [./grain_vol_35]
      type = VolumePostprocessor
      block = '35'
    [../]
    [./grain_vol_36]
      type = VolumePostprocessor
      block = '36'
    [../]
    [./grain_vol_37]
      type = VolumePostprocessor
      block = '37'
    [../]
    [./grain_vol_38]
      type = VolumePostprocessor
      block = '38'
    [../]
    [./grain_vol_39]
      type = VolumePostprocessor
      block = '39'
    [../]
    [./grain_vol_40]
      type = VolumePostprocessor
      block = '40'
    [../]
    [./grain_vol_41]
      type = VolumePostprocessor
      block = '41'
    [../]
    [./grain_vol_42]
      type = VolumePostprocessor
      block = '42'
    [../]
    [./grain_vol_43]
      type = VolumePostprocessor
      block = '43'
    [../]
    [./grain_vol_44]
      type = VolumePostprocessor
      block = '44'
    [../]
    [./grain_vol_45]
      type = VolumePostprocessor
      block = '45'
    [../]
    [./grain_vol_46]
      type = VolumePostprocessor
      block = '46'
    [../]
    [./grain_vol_47]
      type = VolumePostprocessor
      block = '47'
    [../]
    [./grain_vol_48]
      type = VolumePostprocessor
      block = '48'
    [../]
    [./grain_vol_49]
      type = VolumePostprocessor
      block = '49'
    [../]
    [./grain_vol_50]
      type = VolumePostprocessor
      block = '50'
    [../]
    [./grain_vol_51]
      type = VolumePostprocessor
      block = '51'
    [../]
    [./grain_vol_52]
      type = VolumePostprocessor
      block = '52'
    [../]
    [./grain_vol_53]
      type = VolumePostprocessor
      block = '53'
    [../]
    [./grain_vol_54]
      type = VolumePostprocessor
      block = '54'
    [../]
    [./grain_vol_55]
      type = VolumePostprocessor
      block = '55'
    [../]
    [./grain_vol_56]
      type = VolumePostprocessor
      block = '56'
    [../]
    [./grain_vol_57]
      type = VolumePostprocessor
      block = '57'
    [../]
    [./grain_vol_58]
      type = VolumePostprocessor
      block = '58'
    [../]
    [./grain_vol_59]
      type = VolumePostprocessor
      block = '59'
    [../]
    [./grain_vol_60]
      type = VolumePostprocessor
      block = '60'
    [../]
    [./grain_vol_61]
      type = VolumePostprocessor
      block = '61'
    [../]
    [./grain_vol_62]
      type = VolumePostprocessor
      block = '62'
    [../]
    [./grain_vol_63]
      type = VolumePostprocessor
      block = '63'
    [../]
    [./grain_vol_64]
      type = VolumePostprocessor
      block = '64'
    [../]
    [./grain_vol_65]
      type = VolumePostprocessor
      block = '65'
    [../]
    [./grain_vol_66]
      type = VolumePostprocessor
      block = '66'
    [../]
    [./grain_vol_67]
      type = VolumePostprocessor
      block = '67'
    [../]
    [./grain_vol_68]
      type = VolumePostprocessor
      block = '68'
    [../]
    [./grain_vol_69]
      type = VolumePostprocessor
      block = '69'
    [../]
    [./grain_vol_70]
      type = VolumePostprocessor
      block = '70'
    [../]
    [./grain_vol_71]
      type = VolumePostprocessor
      block = '71'
    [../]
    [./grain_vol_72]
      type = VolumePostprocessor
      block = '72'
    [../]
    [./grain_vol_73]
      type = VolumePostprocessor
      block = '73'
    [../]
    [./grain_vol_74]
      type = VolumePostprocessor
      block = '74'
    [../]
    [./grain_vol_75]
      type = VolumePostprocessor
      block = '75'
    [../]
    [./grain_vol_76]
      type = VolumePostprocessor
      block = '76'
    [../]
    [./grain_vol_77]
      type = VolumePostprocessor
      block = '77'
    [../]
    [./grain_vol_78]
      type = VolumePostprocessor
      block = '78'
    [../]
    [./grain_vol_79]
      type = VolumePostprocessor
      block = '79'
    [../]
    [./grain_vol_80]
      type = VolumePostprocessor
      block = '80'
    [../]
    [./grain_vol_81]
      type = VolumePostprocessor
      block = '81'
    [../]
    [./grain_vol_82]
      type = VolumePostprocessor
      block = '82'
    [../]
    [./grain_vol_83]
      type = VolumePostprocessor
      block = '83'
    [../]
    [./grain_vol_84]
      type = VolumePostprocessor
      block = '84'
    [../]
    [./grain_vol_85]
      type = VolumePostprocessor
      block = '85'
    [../]
    [./grain_vol_86]
      type = VolumePostprocessor
      block = '86'
    [../]
    [./grain_vol_87]
      type = VolumePostprocessor
      block = '87'
    [../]
    [./grain_vol_88]
      type = VolumePostprocessor
      block = '88'
    [../]
    [./grain_vol_89]
      type = VolumePostprocessor
      block = '89'
    [../]
    [./grain_vol_90]
      type = VolumePostprocessor
      block = '90'
    [../]
    [./grain_vol_91]
      type = VolumePostprocessor
      block = '91'
    [../]
    [./grain_vol_92]
      type = VolumePostprocessor
      block = '92'
    [../]
    [./grain_vol_93]
      type = VolumePostprocessor
      block = '93'
    [../]
    [./grain_vol_94]
      type = VolumePostprocessor
      block = '94'
    [../]
    [./grain_vol_95]
      type = VolumePostprocessor
      block = '95'
    [../]
    [./grain_vol_96]
      type = VolumePostprocessor
      block = '96'
    [../]
    [./grain_vol_97]
      type = VolumePostprocessor
      block = '97'
    [../]
    [./grain_vol_98]
      type = VolumePostprocessor
      block = '98'
    [../]
    [./grain_vol_99]
      type = VolumePostprocessor
      block = '99'
    [../]
    [./grain_vol_100]
      type = VolumePostprocessor
      block = '100'
    [../]
    [./grain_vol_101]
      type = VolumePostprocessor
      block = '101'
    [../]
    [./grain_vol_102]
      type = VolumePostprocessor
      block = '102'
    [../]
    [./grain_vol_103]
      type = VolumePostprocessor
      block = '103'
    [../]
    [./grain_vol_104]
      type = VolumePostprocessor
      block = '104'
    [../]
    [./grain_vol_105]
      type = VolumePostprocessor
      block = '105'
    [../]
    [./grain_vol_106]
      type = VolumePostprocessor
      block = '106'
    [../]
    [./grain_vol_107]
      type = VolumePostprocessor
      block = '107'
    [../]
    [./grain_vol_108]
      type = VolumePostprocessor
      block = '108'
    [../]
    [./grain_vol_109]
      type = VolumePostprocessor
      block = '109'
    [../]
    [./grain_vol_110]
      type = VolumePostprocessor
      block = '110'
    [../]
    [./grain_vol_111]
      type = VolumePostprocessor
      block = '111'
    [../]
    [./grain_vol_112]
      type = VolumePostprocessor
      block = '112'
    [../]
    [./grain_vol_113]
      type = VolumePostprocessor
      block = '113'
    [../]

    #---------------------#
    #------- BREAK -------#
    #---------------------#

    #[./grain_size_1]
    #  type = GrainSize
    #  block = '1'
    #  vol = grain_vol_1
    #[../]
    #[./grain_size_2]
    #  type = GrainSize
    #  block = '2'
    #  vol = grain_vol_2
    #[../]
    #[./grain_size_3]
    #  type = GrainSize
    #  block = '3'
    #  vol = grain_vol_3
    #[../]
    #[./grain_size_4]
    #  type = GrainSize
    #  block = '4'
    #  vol = grain_vol_4
    #[../]
    #[./grain_size_5]
    #  type = GrainSize
    #  block = '5'
    #  vol = grain_vol_5
    #[../]
    #[./grain_size_6]
    #  type = GrainSize
    #  block = '6'
    #  vol = grain_vol_6
    #[../]
    #[./grain_size_7]
    #  type = GrainSize
    #  block = '7'
    #  vol = grain_vol_7
    #[../]
    #[./grain_size_8]
    #  type = GrainSize
    #  block = '8'
    #  vol = grain_vol_8
    #[../]
    #[./grain_size_9]
    #  type = GrainSize
    #  block = '9'
    #  vol = grain_vol_9
    #[../]
    #[./grain_size_10]
    #  type = GrainSize
    #  block = '10'
    #  vol = grain_vol_10
    #[../]
    #[./grain_size_11]
    #  type = GrainSize
    #  block = '11'
    #  vol = grain_vol_11
    #[../]
    #[./grain_size_12]
    #  type = GrainSize
    #  block = '12'
    #  vol = grain_vol_12
    #[../]
    #[./grain_size_13]
    #  type = GrainSize
    #  block = '13'
    #  vol = grain_vol_13
    #[../]
    #[./grain_size_14]
    #  type = GrainSize
    #  block = '14'
    #  vol = grain_vol_14
    #[../]
    #[./grain_size_15]
    #  type = GrainSize
    #  block = '15'
    #  vol = grain_vol_15
    #[../]
    #[./grain_size_16]
    #  type = GrainSize
    #  block = '16'
    #  vol = grain_vol_16
    #[../]
    #[./grain_size_17]
    #  type = GrainSize
    #  block = '17'
    #  vol = grain_vol_17
    #[../]
    #[./grain_size_18]
    #  type = GrainSize
    #  block = '18'
    #  vol = grain_vol_18
    #[../]
    #[./grain_size_19]
    #  type = GrainSize
    #  block = '19'
    #  vol = grain_vol_19
    #[../]
    #[./grain_size_20]
    #  type = GrainSize
    #  block = '20'
    #  vol = grain_vol_20
    #[../]
    #[./grain_size_21]
    #  type = GrainSize
    #  block = '21'
    #  vol = grain_vol_21
    #[../]
    #[./grain_size_22]
    #  type = GrainSize
    #  block = '22'
    #  vol = grain_vol_22
    #[../]
    #[./grain_size_23]
    #  type = GrainSize
    #  block = '23'
    #  vol = grain_vol_23
    #[../]
    #[./grain_size_24]
    #  type = GrainSize
    #  block = '24'
    #  vol = grain_vol_24
    #[../]
    #[./grain_size_25]
    #  type = GrainSize
    #  block = '25'
    #  vol = grain_vol_25
    #[../]
    #[./grain_size_26]
    #  type = GrainSize
    #  block = '26'
    #  vol = grain_vol_26
    #[../]
    #[./grain_size_27]
    #  type = GrainSize
    #  block = '27'
    #  vol = grain_vol_27
    #[../]
    #[./grain_size_28]
    #  type = GrainSize
    #  block = '28'
    #  vol = grain_vol_28
    #[../]
    #[./grain_size_29]
    #  type = GrainSize
    #  block = '29'
    #  vol = grain_vol_29
    #[../]
    #[./grain_size_30]
    #  type = GrainSize
    #  block = '30'
    #  vol = grain_vol_30
    #[../]
    #[./grain_size_31]
    #  type = GrainSize
    #  block = '31'
    #  vol = grain_vol_31
    #[../]
    #[./grain_size_32]
    #  type = GrainSize
    #  block = '32'
    #  vol = grain_vol_32
    #[../]
    #[./grain_size_33]
    #  type = GrainSize
    #  block = '33'
    #  vol = grain_vol_33
    #[../]
    #[./grain_size_34]
    #  type = GrainSize
    #  block = '34'
    #  vol = grain_vol_34
    #[../]
    #[./grain_size_35]
    #  type = GrainSize
    #  block = '35'
    #  vol = grain_vol_35
    #[../]
    #[./grain_size_36]
    #  type = GrainSize
    #  block = '36'
    #  vol = grain_vol_36
    #[../]
    #[./grain_size_37]
    #  type = GrainSize
    #  block = '37'
    #  vol = grain_vol_37
    #[../]
    #[./grain_size_38]
    #  type = GrainSize
    #  block = '38'
    #  vol = grain_vol_38
    #[../]
    #[./grain_size_39]
    #  type = GrainSize
    #  block = '39'
    #  vol = grain_vol_39
    #[../]
    #[./grain_size_40]
    #  type = GrainSize
    #  block = '40'
    #  vol = grain_vol_40
    #[../]
    #[./grain_size_41]
    #  type = GrainSize
    #  block = '41'
    #  vol = grain_vol_41
    #[../]
    #[./grain_size_42]
    #  type = GrainSize
    #  block = '42'
    #  vol = grain_vol_42
    #[../]
    #[./grain_size_43]
    #  type = GrainSize
    #  block = '43'
    #  vol = grain_vol_43
    #[../]
    #[./grain_size_44]
    #  type = GrainSize
    #  block = '44'
    #  vol = grain_vol_44
    #[../]
    #[./grain_size_45]
    #  type = GrainSize
    #  block = '45'
    #  vol = grain_vol_45
    #[../]
    #[./grain_size_46]
    #  type = GrainSize
    #  block = '46'
    #  vol = grain_vol_46
    #[../]
    #[./grain_size_47]
    #  type = GrainSize
    #  block = '47'
    #  vol = grain_vol_47
    #[../]
    #[./grain_size_48]
    #  type = GrainSize
    #  block = '48'
    #  vol = grain_vol_48
    #[../]
    #[./grain_size_49]
    #  type = GrainSize
    #  block = '49'
    #  vol = grain_vol_49
    #[../]
    #[./grain_size_50]
    #  type = GrainSize
    #  block = '50'
    #  vol = grain_vol_50
    #[../]
    #[./grain_size_51]
    #  type = GrainSize
    #  block = '51'
    #  vol = grain_vol_51
    #[../]
    #[./grain_size_52]
    #  type = GrainSize
    #  block = '52'
    #  vol = grain_vol_52
    #[../]
    #[./grain_size_53]
    #  type = GrainSize
    #  block = '53'
    #  vol = grain_vol_53
    #[../]
    #[./grain_size_54]
    #  type = GrainSize
    #  block = '54'
    #  vol = grain_vol_54
    #[../]
    #[./grain_size_55]
    #  type = GrainSize
    #  block = '55'
    #  vol = grain_vol_55
    #[../]
    #[./grain_size_56]
    #  type = GrainSize
    #  block = '56'
    #  vol = grain_vol_56
    #[../]
    #[./grain_size_57]
    #  type = GrainSize
    #  block = '57'
    #  vol = grain_vol_57
    #[../]
    #[./grain_size_58]
    #  type = GrainSize
    #  block = '58'
    #  vol = grain_vol_58
    #[../]
    #[./grain_size_59]
    #  type = GrainSize
    #  block = '59'
    #  vol = grain_vol_59
    #[../]
    #[./grain_size_60]
    #  type = GrainSize
    #  block = '60'
    #  vol = grain_vol_60
    #[../]
    #[./grain_size_61]
    #  type = GrainSize
    #  block = '61'
    #  vol = grain_vol_61
    #[../]
    #[./grain_size_62]
    #  type = GrainSize
    #  block = '62'
    #  vol = grain_vol_62
    #[../]
    #[./grain_size_63]
    #  type = GrainSize
    #  block = '63'
    #  vol = grain_vol_63
    #[../]
    #[./grain_size_64]
    #  type = GrainSize
    #  block = '64'
    #  vol = grain_vol_64
    #[../]
    #[./grain_size_65]
    #  type = GrainSize
    #  block = '65'
    #  vol = grain_vol_65
    #[../]
    #[./grain_size_66]
    #  type = GrainSize
    #  block = '66'
    #  vol = grain_vol_66
    #[../]
    #[./grain_size_67]
    #  type = GrainSize
    #  block = '67'
    #  vol = grain_vol_67
    #[../]
    #[./grain_size_68]
    #  type = GrainSize
    #  block = '68'
    #  vol = grain_vol_68
    #[../]
    #[./grain_size_69]
    #  type = GrainSize
    #  block = '69'
    #  vol = grain_vol_69
    #[../]
    #[./grain_size_70]
    #  type = GrainSize
    #  block = '70'
    #  vol = grain_vol_70
    #[../]
    #[./grain_size_71]
    #  type = GrainSize
    #  block = '71'
    #  vol = grain_vol_71
    #[../]
    #[./grain_size_72]
    #  type = GrainSize
    #  block = '72'
    #  vol = grain_vol_72
    #[../]
    #[./grain_size_73]
    #  type = GrainSize
    #  block = '73'
    #  vol = grain_vol_73
    #[../]
    #[./grain_size_74]
    #  type = GrainSize
    #  block = '74'
    #  vol = grain_vol_74
    #[../]
    #[./grain_size_75]
    #  type = GrainSize
    #  block = '75'
    #  vol = grain_vol_75
    #[../]
    #[./grain_size_76]
    #  type = GrainSize
    #  block = '76'
    #  vol = grain_vol_76
    #[../]
    #[./grain_size_77]
    #  type = GrainSize
    #  block = '77'
    #  vol = grain_vol_77
    #[../]
    #[./grain_size_78]
    #  type = GrainSize
    #  block = '78'
    #  vol = grain_vol_78
    #[../]
    #[./grain_size_79]
    #  type = GrainSize
    #  block = '79'
    #  vol = grain_vol_79
    #[../]
    #[./grain_size_80]
    #  type = GrainSize
    #  block = '80'
    #  vol = grain_vol_80
    #[../]
    #[./grain_size_81]
    #  type = GrainSize
    #  block = '81'
    #  vol = grain_vol_81
    #[../]
    #[./grain_size_82]
    #  type = GrainSize
    #  block = '82'
    #  vol = grain_vol_82
    #[../]
    #[./grain_size_83]
    #  type = GrainSize
    #  block = '83'
    #  vol = grain_vol_83
    #[../]
    #[./grain_size_84]
    #  type = GrainSize
    #  block = '84'
    #  vol = grain_vol_84
    #[../]
    #[./grain_size_85]
    #  type = GrainSize
    #  block = '85'
    #  vol = grain_vol_85
    #[../]
    #[./grain_size_86]
    #  type = GrainSize
    #  block = '86'
    #  vol = grain_vol_86
    #[../]
    #[./grain_size_87]
    #  type = GrainSize
    #  block = '87'
    #  vol = grain_vol_87
    #[../]
    #[./grain_size_88]
    #  type = GrainSize
    #  block = '88'
    #  vol = grain_vol_88
    #[../]
    #[./grain_size_89]
    #  type = GrainSize
    #  block = '89'
    #  vol = grain_vol_89
    #[../]
    #[./grain_size_90]
    #  type = GrainSize
    #  block = '90'
    #  vol = grain_vol_90
    #[../]
    #[./grain_size_91]
    #  type = GrainSize
    #  block = '91'
    #  vol = grain_vol_91
    #[../]
    #[./grain_size_92]
    #  type = GrainSize
    #  block = '92'
    #  vol = grain_vol_92
    #[../]
    #[./grain_size_93]
    #  type = GrainSize
    #  block = '93'
    #  vol = grain_vol_93
    #[../]
    #[./grain_size_94]
    #  type = GrainSize
    #  block = '94'
    #  vol = grain_vol_94
    #[../]
    #[./grain_size_95]
    #  type = GrainSize
    #  block = '95'
    #  vol = grain_vol_95
    #[../]
    #[./grain_size_96]
    #  type = GrainSize
    #  block = '96'
    #  vol = grain_vol_96
    #[../]
    #[./grain_size_97]
    #  type = GrainSize
    #  block = '97'
    #  vol = grain_vol_97
    #[../]
    #[./grain_size_98]
    #  type = GrainSize
    #  block = '98'
    #  vol = grain_vol_98
    #[../]
    #[./grain_size_99]
    #  type = GrainSize
    #  block = '99'
    #  vol = grain_vol_99
    #[../]
    #[./grain_size_100]
    #  type = GrainSize
    #  block = '100'
    #  vol = grain_vol_100
    #[../]
    #[./grain_size_101]
    #  type = GrainSize
    #  block = '101'
    #  vol = grain_vol_101
    #[../]
    #[./grain_size_102]
    #  type = GrainSize
    #  block = '102'
    #  vol = grain_vol_102
    #[../]
    #[./grain_size_103]
    #  type = GrainSize
    #  block = '103'
    #  vol = grain_vol_103
    #[../]
    #[./grain_size_104]
    #  type = GrainSize
    #  block = '104'
    #  vol = grain_vol_104
    #[../]
    #[./grain_size_105]
    #  type = GrainSize
    #  block = '105'
    #  vol = grain_vol_105
    #[../]
    #[./grain_size_106]
    #  type = GrainSize
    #  block = '106'
    #  vol = grain_vol_106
    #[../]
    #[./grain_size_107]
    #  type = GrainSize
    #  block = '107'
    #  vol = grain_vol_107
    #[../]
    #[./grain_size_108]
    #  type = GrainSize
    #  block = '108'
    #  vol = grain_vol_108
    #[../]
    #[./grain_size_109]
    #  type = GrainSize
    #  block = '109'
    #  vol = grain_vol_109
    #[../]
    #[./grain_size_110]
    #  type = GrainSize
    #  block = '110'
    #  vol = grain_vol_110
    #[../]
    #[./grain_size_111]
    #  type = GrainSize
    #  block = '111'
    #  vol = grain_vol_111
    #[../]
    #[./grain_size_112]
    #  type = GrainSize
    #  block = '112'
    #  vol = grain_vol_112
    #[../]
    #[./grain_size_113]
    #  type = GrainSize
    #  block = '113'
    #  vol = grain_vol_113
    #[../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -ksp_snes_ew'
    petsc_options_iname = '-ksp_gmres_restart  -snes_atol -snes_rtol -ksp_rtol -pc_type -pc_hypre_type'
    petsc_options_value = '    121                1e-10      1e-8      1e-8    hypre    boomeramg'
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK' #"PJFNK, JFNK, NEWTON"
[]

[Outputs]
  # This file takes 231 seconds to solve on 8 processors.
  # Don't need to print linear residuals. We don't have convergence issues
  print_linear_residuals = false
  print_perf_log = true
  [./console]
    type = Console
    execute_postprocessors_on = none
  [../]
  [./out]
    type = Exodus
    file_base = out_HCP_grains_test
    elemental_as_nodal = true
  [../]
  [./out_csv]
    type = CSV
    file_base = out_HCP_grains_test_csv
  [../]
[]
