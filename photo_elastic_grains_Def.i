
[Mesh]
  file = 3D_HCP_256.e
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
    block = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15  16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 100 101 102 103 104 105 106 107 108 109 110 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243'
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

  [./elasticity_tensor_47]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 0.77
    euler_angle_2 = 0.77
    euler_angle_3 = -236.0
    block = '47'
  [../]
  [./strain_47]
    type = ComputeSmallStrain
    block = '47'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_47]
    type = ComputeLinearElasticStress
    block = '47'
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

  [./elasticity_tensor_71]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -99.0
    euler_angle_2 = -18.0
    euler_angle_3 = 63.0
    block = '71'
  [../]
  [./strain_71]
    type = ComputeSmallStrain
    block = '71'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_71]
    type = ComputeLinearElasticStress
    block = '71'
  [../]

  [./elasticity_tensor_72]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -63.0
    euler_angle_2 = 25.0
    euler_angle_3 = 78.0
    block = '72'
  [../]
  [./strain_72]
    type = ComputeSmallStrain
    block = '72'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_72]
    type = ComputeLinearElasticStress
    block = '72'
  [../]

  [./elasticity_tensor_73]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 256.0
    euler_angle_2 = -3.0
    euler_angle_3 = 33.0
    block = '73'
  [../]
  [./strain_73]
    type = ComputeSmallStrain
    block = '73'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_73]
    type = ComputeLinearElasticStress
    block = '73'
  [../]

  [./elasticity_tensor_74]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 282.0
    euler_angle_2 = 182.0
    euler_angle_3 = 6.3
    block = '74'
  [../]
  [./strain_74]
    type = ComputeSmallStrain
    block = '74'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_74]
    type = ComputeLinearElasticStress
    block = '74'
  [../]

  [./elasticity_tensor_75]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 81.0
    euler_angle_2 = 85.0
    euler_angle_3 = 89.0
    block = '75'
  [../]
  [./strain_75]
    type = ComputeSmallStrain
    block = '75'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_75]
    type = ComputeLinearElasticStress
    block = '75'
  [../]

  [./elasticity_tensor_76]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 8.2
    euler_angle_2 = -8.2
    euler_angle_3 = 63.0
    block = '76'
  [../]
  [./strain_76]
    type = ComputeSmallStrain
    block = '76'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_76]
    type = ComputeLinearElasticStress
    block = '76'
  [../]

  [./elasticity_tensor_77]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 76.0
    euler_angle_2 = 48.0
    euler_angle_3 = -29.0
    block = '77'
  [../]
  [./strain_77]
    type = ComputeSmallStrain
    block = '77'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_77]
    type = ComputeLinearElasticStress
    block = '77'
  [../]

  [./elasticity_tensor_78]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -9.2
    euler_angle_2 = 6.8
    euler_angle_3 = 110.0
    block = '78'
  [../]
  [./strain_78]
    type = ComputeSmallStrain
    block = '78'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_78]
    type = ComputeLinearElasticStress
    block = '78'
  [../]

  [./elasticity_tensor_79]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 12.0
    euler_angle_2 = -12.0
    euler_angle_3 = 146.0
    block = '79'
  [../]
  [./strain_79]
    type = ComputeSmallStrain
    block = '79'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_79]
    type = ComputeLinearElasticStress
    block = '79'
  [../]

  [./elasticity_tensor_80]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -129.0
    euler_angle_2 = 0.0
    euler_angle_3 = 13.0
    block = '80'
  [../]
  [./strain_80]
    type = ComputeSmallStrain
    block = '80'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_80]
    type = ComputeLinearElasticStress
    block = '80'
  [../]


  [./elasticity_tensor_81]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 24.0
    euler_angle_2 = 24.0
    euler_angle_3 = -69.0
    block = '81'
  [../]
  [./strain_81]
    type = ComputeSmallStrain
    block = '81'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_81]
    type = ComputeLinearElasticStress
    block = '81'
  [../]

  [./elasticity_tensor_82]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 84.0
    euler_angle_2 = 18.0
    euler_angle_3 = 125.0
    block = '82'
  [../]
  [./strain_82]
    type = ComputeSmallStrain
    block = '82'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_82]
    type = ComputeLinearElasticStress
    block = '82'
  [../]

  [./elasticity_tensor_83]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 105.0
    euler_angle_2 = 150.0
    euler_angle_3 = -26
    block = '83'
  [../]
  [./strain_83]
    type = ComputeSmallStrain
    block = '83'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_83]
    type = ComputeLinearElasticStress
    block = '83'
  [../]

  [./elasticity_tensor_84]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -105.0
    euler_angle_2 = -75.0
    euler_angle_3 = 10.0
    block = '84'
  [../]
  [./strain_84]
    type = ComputeSmallStrain
    block = '84'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_84]
    type = ComputeLinearElasticStress
    block = '84'
  [../]

  [./elasticity_tensor_85]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -46.0
    euler_angle_2 = -48.0
    euler_angle_3 = 49.0
    block = '85'
  [../]
  [./strain_85]
    type = ComputeSmallStrain
    block = '85'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_85]
    type = ComputeLinearElasticStress
    block = '85'
  [../]

  [./elasticity_tensor_86]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 87.0
    euler_angle_2 = 86.0
    euler_angle_3 = 38.0
    block = '86'
  [../]
  [./strain_86]
    type = ComputeSmallStrain
    block = '86'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_86]
    type = ComputeLinearElasticStress
    block = '86'
  [../]

  [./elasticity_tensor_87]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 108.0
    euler_angle_2 = 10.6
    euler_angle_3 = -2.4
    block = '87'
  [../]
  [./strain_87]
    type = ComputeSmallStrain
    block = '87'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_87]
    type = ComputeLinearElasticStress
    block = '87'
  [../]

  [./elasticity_tensor_88]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 56.0
    euler_angle_2 = 56.0
    euler_angle_3 = 98.0
    block = '88'
  [../]
  [./strain_88]
    type = ComputeSmallStrain
    block = '88'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_88]
    type = ComputeLinearElasticStress
    block = '88'
  [../]

  [./elasticity_tensor_89]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 23.0
    euler_angle_2 = -36.0
    euler_angle_3 = 39.0
    block = '89'
  [../]
  [./strain_89]
    type = ComputeSmallStrain
    block = '89'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_89]
    type = ComputeLinearElasticStress
    block = '89'
  [../]

  [./elasticity_tensor_90]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 5.0
    euler_angle_2 = 7.0
    euler_angle_3 = 198.0
    block = '90'
  [../]
  [./strain_90]
    type = ComputeSmallStrain
    block = '90'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_90]
    type = ComputeLinearElasticStress
    block = '90'
  [../]

  [./elasticity_tensor_91]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -9.0
    euler_angle_2 = 15.5
    euler_angle_3 = 126.0
    block = '91'
  [../]
  [./strain_91]
    type = ComputeSmallStrain
    block = '91'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_91]
    type = ComputeLinearElasticStress
    block = '91'
  [../]

  [./elasticity_tensor_92]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 37.0
    euler_angle_2 = 0.4
    euler_angle_3 = 116.0
    block = '92'
  [../]
  [./strain_92]
    type = ComputeSmallStrain
    block = '92'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_92]
    type = ComputeLinearElasticStress
    block = '92'
  [../]

  [./elasticity_tensor_93]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '93'
  [../]
  [./strain_93]
    type = ComputeSmallStrain
    block = '93'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_93]
    type = ComputeLinearElasticStress
    block = '93'
  [../]

  [./elasticity_tensor_94]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 16.0
    euler_angle_2 = 78.0
    euler_angle_3 = 137.0
    block = '94'
  [../]
  [./strain_94]
    type = ComputeSmallStrain
    block = '94'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_94]
    type = ComputeLinearElasticStress
    block = '94'
  [../]

  [./elasticity_tensor_95]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 145.0
    euler_angle_2 = 136.0
    euler_angle_3 = 157.0
    block = '95'
  [../]
  [./strain_95]
    type = ComputeSmallStrain
    block = '95'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_95]
    type = ComputeLinearElasticStress
    block = '95'
  [../]

  [./elasticity_tensor_96]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -15.0
    euler_angle_2 = -16.0
    euler_angle_3 = -17.0
    block = '96'
  [../]
  [./strain_96]
    type = ComputeSmallStrain
    block = '96'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_96]
    type = ComputeLinearElasticStress
    block = '96'
  [../]


  [./elasticity_tensor_97]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 0.5
    euler_angle_2 = 0.5
    euler_angle_3 = 3.5
    block = '97'
  [../]
  [./strain_97]
    type = ComputeSmallStrain
    block = '97'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_97]
    type = ComputeLinearElasticStress
    block = '97'
  [../]

  [./elasticity_tensor_98]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 65.5
    euler_angle_2 = 85.5
    euler_angle_3 = 13.5
    block = '98'
  [../]
  [./strain_98]
    type = ComputeSmallStrain
    block = '98'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_98]
    type = ComputeLinearElasticStress
    block = '98'
  [../]


  [./elasticity_tensor_100]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 25.0
    euler_angle_2 = 45.0
    euler_angle_3 = -69.0
    block = '100'
  [../]
  [./strain_100]
    type = ComputeSmallStrain
    block = '100'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_100]
    type = ComputeLinearElasticStress
    block = '100'
  [../]


  [./elasticity_tensor_101]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 50.0
    euler_angle_2 = 50.0
    euler_angle_3 = 36.0
    block = '101'
  [../]
  [./strain_101]
    type = ComputeSmallStrain
    block = '101'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_101]
    type = ComputeLinearElasticStress
    block = '101'
  [../]

  [./elasticity_tensor_102]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 106.0
    euler_angle_2 = -10.0
    euler_angle_3 = 109.0
    block = '102'
  [../]
  [./strain_102]
    type = ComputeSmallStrain
    block = '102'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_102]
    type = ComputeLinearElasticStress
    block = '102'
  [../]

  [./elasticity_tensor_103]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 17.0
    euler_angle_2 = 18.0
    euler_angle_3 = -19.0
    block = '103'
  [../]
  [./strain_103]
    type = ComputeSmallStrain
    block = '103'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_103]
    type = ComputeLinearElasticStress
    block = '103'
  [../]

  [./elasticity_tensor_104]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -27.2
    euler_angle_2 = 125.0
    euler_angle_3 = -9.0
    block = '104'
  [../]
  [./strain_104]
    type = ComputeSmallStrain
    block = '104'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_104]
    type = ComputeLinearElasticStress
    block = '104'
  [../]

  [./elasticity_tensor_105]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 19.0
    euler_angle_2 = 19.0
    euler_angle_3 = -49.0
    block = '105'
  [../]
  [./strain_105]
    type = ComputeSmallStrain
    block = '105'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_105]
    type = ComputeLinearElasticStress
    block = '105'
  [../]

  [./elasticity_tensor_106]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 67.0
    euler_angle_2 = 8.0
    euler_angle_3 = -81.0
    block = '106'
  [../]
  [./strain_106]
    type = ComputeSmallStrain
    block = '106'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_106]
    type = ComputeLinearElasticStress
    block = '106'
  [../]

  [./elasticity_tensor_107]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -16.0
    euler_angle_2 = 28.0
    euler_angle_3 = 0.0
    block = '107'
  [../]
  [./strain_107]
    type = ComputeSmallStrain
    block = '107'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_107]
    type = ComputeLinearElasticStress
    block = '107'
  [../]

  [./elasticity_tensor_108]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -6.0
    euler_angle_2 = -8.0
    euler_angle_3 = 19.0
    block = '108'
  [../]
  [./strain_108]
    type = ComputeSmallStrain
    block = '108'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_108]
    type = ComputeLinearElasticStress
    block = '108'
  [../]

  [./elasticity_tensor_109]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 100.0
    euler_angle_2 = 100.0
    euler_angle_3 = 79.0
    block = '109'
  [../]
  [./strain_109]
    type = ComputeSmallStrain
    block = '109'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_109]
    type = ComputeLinearElasticStress
    block = '109'
  [../]


  [./elasticity_tensor_110]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 260.0
    euler_angle_2 = 84.0
    euler_angle_3 = 89.0
    block = '110'
  [../]
  [./strain_110]
    type = ComputeSmallStrain
    block = '110'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_110]
    type = ComputeLinearElasticStress
    block = '110'
  [../]


  [./elasticity_tensor_111]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -66.0
    euler_angle_2 = 18.0
    euler_angle_3 = -77.0
    block = '111'
  [../]
  [./strain_111]
    type = ComputeSmallStrain
    block = '111'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_111]
    type = ComputeLinearElasticStress
    block = '111'
  [../]

  [./elasticity_tensor_112]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 177.0
    euler_angle_2 = 177.0
    euler_angle_3 = 177.0
    block = '112'
  [../]
  [./strain_112]
    type = ComputeSmallStrain
    block = '112'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_112]
    type = ComputeLinearElasticStress
    block = '112'
  [../]

  [./elasticity_tensor_113]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 84.0
    euler_angle_2 = -4.0
    euler_angle_3 = 196.0
    block = '113'
  [../]
  [./strain_113]
    type = ComputeSmallStrain
    block = '113'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_113]
    type = ComputeLinearElasticStress
    block = '113'
  [../]

  [./elasticity_tensor_114]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 15.0
    euler_angle_2 = 105.0
    euler_angle_3 = 1.5
    block = '114'
  [../]
  [./strain_114]
    type = ComputeSmallStrain
    block = '114'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_114]
    type = ComputeLinearElasticStress
    block = '114'
  [../]


  [./elasticity_tensor_115]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 23.0
    euler_angle_2 = 235.0
    euler_angle_3 = 2.32
    block = '115'
  [../]
  [./strain_115]
    type = ComputeSmallStrain
    block = '115'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_115]
    type = ComputeLinearElasticStress
    block = '115'
  [../]

  [./elasticity_tensor_116]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 75.0
    euler_angle_2 = 67.0
    euler_angle_3 = 1.56
    block = '116'
  [../]
  [./strain_116]
    type = ComputeSmallStrain
    block = '116'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_116]
    type = ComputeLinearElasticStress
    block = '116'
  [../]

  [./elasticity_tensor_117]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -0.2
    euler_angle_2 = 186.3
    euler_angle_3 = 55.6
    block = '117'
  [../]
  [./strain_117]
    type = ComputeSmallStrain
    block = '117'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_117]
    type = ComputeLinearElasticStress
    block = '117'
  [../]

  [./elasticity_tensor_118]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -6.28
    euler_angle_2 = 77.5
    euler_angle_3 = 77.5
    block = '118'
  [../]
  [./strain_118]
    type = ComputeSmallStrain
    block = '118'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_118]
    type = ComputeLinearElasticStress
    block = '118'
  [../]

  [./elasticity_tensor_119]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -0.2
    euler_angle_2 = 186.3
    euler_angle_3 = 55.6
    block = '119'
  [../]
  [./strain_119]
    type = ComputeSmallStrain
    block = '119'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_119]
    type = ComputeLinearElasticStress
    block = '119'
  [../]

  [./elasticity_tensor_120]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 200.0
    euler_angle_2 = 200.0
    euler_angle_3 = 5.6
    block = '120'
  [../]
  [./strain_120]
    type = ComputeSmallStrain
    block = '120'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_120]
    type = ComputeLinearElasticStress
    block = '120'
  [../]

  [./elasticity_tensor_121]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 75.3
    euler_angle_2 = 76.8
    euler_angle_3 = -15.6
    block = '121'
  [../]
  [./strain_121]
    type = ComputeSmallStrain
    block = '121'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_121]
    type = ComputeLinearElasticStress
    block = '121'
  [../]

  [./elasticity_tensor_122]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
    block = '122'
  [../]
  [./strain_122]
    type = ComputeSmallStrain
    block = '122'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_122]
    type = ComputeLinearElasticStress
    block = '122'
  [../]

  [./elasticity_tensor_123]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 66.0
    euler_angle_2 = 87.0
    euler_angle_3 = -67.5
    block = '123'
  [../]
  [./strain_123]
    type = ComputeSmallStrain
    block = '123'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_123]
    type = ComputeLinearElasticStress
    block = '123'
  [../]

  [./elasticity_tensor_124]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '124'
  [../]
  [./strain_124]
    type = ComputeSmallStrain
    block = '124'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_124]
    type = ComputeLinearElasticStress
    block = '124'
  [../]


  [./elasticity_tensor_125]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 210.0
    euler_angle_2 = 215.0
    euler_angle_3 = -15.6
    block = '125'
  [../]
  [./strain_125]
    type = ComputeSmallStrain
    block = '125'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_125]
    type = ComputeLinearElasticStress
    block = '125'
  [../]


  [./elasticity_tensor_126]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '126'
  [../]
  [./strain_126]
    type = ComputeSmallStrain
    block = '126'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_126]
    type = ComputeLinearElasticStress
    block = '126'
  [../]


  [./elasticity_tensor_127]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '127'
  [../]
  [./strain_127]
    type = ComputeSmallStrain
    block = '127'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_127]
    type = ComputeLinearElasticStress
    block = '127'
  [../]

  [./elasticity_tensor_128]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '128'
  [../]
  [./strain_128]
    type = ComputeSmallStrain
    block = '128'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_128]
    type = ComputeLinearElasticStress
    block = '128'
  [../]

  [./elasticity_tensor_129]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '129'
  [../]
  [./strain_129]
    type = ComputeSmallStrain
    block = '129'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_129]
    type = ComputeLinearElasticStress
    block = '129'
  [../]

  [./elasticity_tensor_130]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '130'
  [../]
  [./strain_130]
    type = ComputeSmallStrain
    block = '130'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_130]
    type = ComputeLinearElasticStress
    block = '130'
  [../]

  [./elasticity_tensor_131]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '131'
  [../]
  [./strain_131]
    type = ComputeSmallStrain
    block = '131'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_131]
    type = ComputeLinearElasticStress
    block = '131'
  [../]

  [./elasticity_tensor_132]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '132'
  [../]
  [./strain_132]
    type = ComputeSmallStrain
    block = '132'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_132]
    type = ComputeLinearElasticStress
    block = '132'
  [../]

  [./elasticity_tensor_133]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '133'
  [../]
  [./strain_133]
    type = ComputeSmallStrain
    block = '133'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_133]
    type = ComputeLinearElasticStress
    block = '133'
  [../]

  [./elasticity_tensor_134]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '134'
  [../]
  [./strain_134]
    type = ComputeSmallStrain
    block = '134'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_134]
    type = ComputeLinearElasticStress
    block = '134'
  [../]

  [./elasticity_tensor_135]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '135'
  [../]
  [./strain_135]
    type = ComputeSmallStrain
    block = '135'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_135]
    type = ComputeLinearElasticStress
    block = '135'
  [../]

  [./elasticity_tensor_136]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '136'
  [../]
  [./strain_136]
    type = ComputeSmallStrain
    block = '136'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_136]
    type = ComputeLinearElasticStress
    block = '136'
  [../]

  [./elasticity_tensor_137]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '137'
  [../]
  [./strain_137]
    type = ComputeSmallStrain
    block = '137'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_137]
    type = ComputeLinearElasticStress
    block = '137'
  [../]

  [./elasticity_tensor_138]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '138'
  [../]
  [./strain_138]
    type = ComputeSmallStrain
    block = '138'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_138]
    type = ComputeLinearElasticStress
    block = '138'
  [../]

  [./elasticity_tensor_139]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '139'
  [../]
  [./strain_139]
    type = ComputeSmallStrain
    block = '139'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_139]
    type = ComputeLinearElasticStress
    block = '139'
  [../]

  [./elasticity_tensor_140]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '140'
  [../]
  [./strain_140]
    type = ComputeSmallStrain
    block = '140'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_140]
    type = ComputeLinearElasticStress
    block = '140'
  [../]

  [./elasticity_tensor_141]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '141'
  [../]
  [./strain_141]
    type = ComputeSmallStrain
    block = '141'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_141]
    type = ComputeLinearElasticStress
    block = '141'
  [../]

  [./elasticity_tensor_142]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '142'
  [../]
  [./strain_142]
    type = ComputeSmallStrain
    block = '142'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_142]
    type = ComputeLinearElasticStress
    block = '142'
  [../]

  [./elasticity_tensor_143]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '143'
  [../]
  [./strain_143]
    type = ComputeSmallStrain
    block = '143'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_143]
    type = ComputeLinearElasticStress
    block = '143'
  [../]

  [./elasticity_tensor_144]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '144'
  [../]
  [./strain_144]
    type = ComputeSmallStrain
    block = '144'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_144]
    type = ComputeLinearElasticStress
    block = '144'
  [../]

  [./elasticity_tensor_145]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '145'
  [../]
  [./strain_145]
    type = ComputeSmallStrain
    block = '145'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_145]
    type = ComputeLinearElasticStress
    block = '145'
  [../]

  [./elasticity_tensor_146]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '146'
  [../]
  [./strain_146]
    type = ComputeSmallStrain
    block = '146'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_146]
    type = ComputeLinearElasticStress
    block = '146'
  [../]

  [./elasticity_tensor_147]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '147'
  [../]
  [./strain_147]
    type = ComputeSmallStrain
    block = '147'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_147]
    type = ComputeLinearElasticStress
    block = '147'
  [../]

  [./elasticity_tensor_148]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '148'
  [../]
  [./strain_148]
    type = ComputeSmallStrain
    block = '148'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_148]
    type = ComputeLinearElasticStress
    block = '148'
  [../]

  [./elasticity_tensor_149]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '149'
  [../]
  [./strain_149]
    type = ComputeSmallStrain
    block = '149'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_149]
    type = ComputeLinearElasticStress
    block = '149'
  [../]

  [./elasticity_tensor_150]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '150'
  [../]
  [./strain_150]
    type = ComputeSmallStrain
    block = '150'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_150]
    type = ComputeLinearElasticStress
    block = '150'
  [../]

  [./elasticity_tensor_151]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '151'
  [../]
  [./strain_151]
    type = ComputeSmallStrain
    block = '151'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_151]
    type = ComputeLinearElasticStress
    block = '151'
  [../]

  [./elasticity_tensor_152]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '152'
  [../]
  [./strain_152]
    type = ComputeSmallStrain
    block = '152'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_152]
    type = ComputeLinearElasticStress
    block = '152'
  [../]

  [./elasticity_tensor_153]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '153'
  [../]
  [./strain_153]
    type = ComputeSmallStrain
    block = '153'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_153]
    type = ComputeLinearElasticStress
    block = '153'
  [../]

  [./elasticity_tensor_154]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '154'
  [../]
  [./strain_154]
    type = ComputeSmallStrain
    block = '154'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_154]
    type = ComputeLinearElasticStress
    block = '154'
  [../]

  [./elasticity_tensor_155]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '155'
  [../]
  [./strain_155]
    type = ComputeSmallStrain
    block = '155'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_155]
    type = ComputeLinearElasticStress
    block = '155'
  [../]

  [./elasticity_tensor_156]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '156'
  [../]
  [./strain_156]
    type = ComputeSmallStrain
    block = '156'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_156]
    type = ComputeLinearElasticStress
    block = '156'
  [../]

  [./elasticity_tensor_157]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '157'
  [../]
  [./strain_157]
    type = ComputeSmallStrain
    block = '157'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_157]
    type = ComputeLinearElasticStress
    block = '157'
  [../]

  [./elasticity_tensor_158]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '158'
  [../]
  [./strain_158]
    type = ComputeSmallStrain
    block = '158'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_158]
    type = ComputeLinearElasticStress
    block = '158'
  [../]

  [./elasticity_tensor_159]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '159'
  [../]
  [./strain_159]
    type = ComputeSmallStrain
    block = '159'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_159]
    type = ComputeLinearElasticStress
    block = '159'
  [../]

  [./elasticity_tensor_160]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '160'
  [../]
  [./strain_160]
    type = ComputeSmallStrain
    block = '160'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_160]
    type = ComputeLinearElasticStress
    block = '160'
  [../]

  [./elasticity_tensor_161]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '161'
  [../]
  [./strain_161]
    type = ComputeSmallStrain
    block = '161'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_161]
    type = ComputeLinearElasticStress
    block = '161'
  [../]

  [./elasticity_tensor_162]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '162'
  [../]
  [./strain_162]
    type = ComputeSmallStrain
    block = '162'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_162]
    type = ComputeLinearElasticStress
    block = '162'
  [../]

  [./elasticity_tensor_163]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '163'
  [../]
  [./strain_163]
    type = ComputeSmallStrain
    block = '163'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_163]
    type = ComputeLinearElasticStress
    block = '163'
  [../]

  [./elasticity_tensor_164]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '164'
  [../]
  [./strain_164]
    type = ComputeSmallStrain
    block = '164'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_164]
    type = ComputeLinearElasticStress
    block = '164'
  [../]

  [./elasticity_tensor_165]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '165'
  [../]
  [./strain_165]
    type = ComputeSmallStrain
    block = '165'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_165]
    type = ComputeLinearElasticStress
    block = '165'
  [../]

  [./elasticity_tensor_166]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '166'
  [../]
  [./strain_166]
    type = ComputeSmallStrain
    block = '166'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_166]
    type = ComputeLinearElasticStress
    block = '166'
  [../]


  [./elasticity_tensor_167]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '167'
  [../]
  [./strain_167]
    type = ComputeSmallStrain
    block = '167'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_167]
    type = ComputeLinearElasticStress
    block = '167'
  [../]

  [./elasticity_tensor_168]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '168'
  [../]
  [./strain_168]
    type = ComputeSmallStrain
    block = '168'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_168]
    type = ComputeLinearElasticStress
    block = '168'
  [../]

  [./elasticity_tensor_169]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '169'
  [../]
  [./strain_169]
    type = ComputeSmallStrain
    block = '169'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_169]
    type = ComputeLinearElasticStress
    block = '169'
  [../]

  [./elasticity_tensor_170]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '170'
  [../]
  [./strain_170]
    type = ComputeSmallStrain
    block = '170'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_170]
    type = ComputeLinearElasticStress
    block = '170'
  [../]

  [./elasticity_tensor_171]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '171'
  [../]
  [./strain_171]
    type = ComputeSmallStrain
    block = '171'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_171]
    type = ComputeLinearElasticStress
    block = '171'
  [../]

  [./elasticity_tensor_172]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '172'
  [../]
  [./strain_172]
    type = ComputeSmallStrain
    block = '172'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_172]
    type = ComputeLinearElasticStress
    block = '172'
  [../]

  [./elasticity_tensor_173]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '173'
  [../]
  [./strain_173]
    type = ComputeSmallStrain
    block = '173'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_173]
    type = ComputeLinearElasticStress
    block = '173'
  [../]

  [./elasticity_tensor_174]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '174'
  [../]
  [./strain_174]
    type = ComputeSmallStrain
    block = '174'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_174]
    type = ComputeLinearElasticStress
    block = '174'
  [../]

  [./elasticity_tensor_175]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '175'
  [../]
  [./strain_175]
    type = ComputeSmallStrain
    block = '175'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_175]
    type = ComputeLinearElasticStress
    block = '175'
  [../]

  [./elasticity_tensor_176]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '176'
  [../]
  [./strain_176]
    type = ComputeSmallStrain
    block = '176'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_176]
    type = ComputeLinearElasticStress
    block = '176'
  [../]

  [./elasticity_tensor_177]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 20.0
    euler_angle_2 = 57.0
    euler_angle_3 = 35.6
    block = '177'
  [../]
  [./strain_177]
    type = ComputeSmallStrain
    block = '177'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_177]
    type = ComputeLinearElasticStress
    block = '177'
  [../]

  [./elasticity_tensor_178]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 28.0
    euler_angle_2 = 156.0
    euler_angle_3 = -45.7
    block = '178'
  [../]
  [./strain_178]
    type = ComputeSmallStrain
    block = '178'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_178]
    type = ComputeLinearElasticStress
    block = '178'
  [../]

  [./elasticity_tensor_179]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '179'
  [../]
  [./strain_179]
    type = ComputeSmallStrain
    block = '179'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_179]
    type = ComputeLinearElasticStress
    block = '179'
  [../]

  [./elasticity_tensor_180]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '180'
  [../]
  [./strain_180]
    type = ComputeSmallStrain
    block = '180'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_180]
    type = ComputeLinearElasticStress
    block = '180'
  [../]

  [./elasticity_tensor_181]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '181'
  [../]
  [./strain_181]
    type = ComputeSmallStrain
    block = '181'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_181]
    type = ComputeLinearElasticStress
    block = '181'
  [../]

  [./elasticity_tensor_182]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '182'
  [../]
  [./strain_182]
    type = ComputeSmallStrain
    block = '182'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_182]
    type = ComputeLinearElasticStress
    block = '182'
  [../]

  [./elasticity_tensor_183]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '183'
  [../]
  [./strain_183]
    type = ComputeSmallStrain
    block = '183'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_183]
    type = ComputeLinearElasticStress
    block = '183'
  [../]

  [./elasticity_tensor_184]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '184'
  [../]
  [./strain_184]
    type = ComputeSmallStrain
    block = '184'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_184]
    type = ComputeLinearElasticStress
    block = '184'
  [../]

  [./elasticity_tensor_185]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '185'
  [../]
  [./strain_185]
    type = ComputeSmallStrain
    block = '185'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_185]
    type = ComputeLinearElasticStress
    block = '185'
  [../]

  [./elasticity_tensor_186]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '186'
  [../]
  [./strain_186]
    type = ComputeSmallStrain
    block = '186'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_186]
    type = ComputeLinearElasticStress
    block = '186'
  [../]

  [./elasticity_tensor_187]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '187'
  [../]
  [./strain_187]
    type = ComputeSmallStrain
    block = '187'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_187]
    type = ComputeLinearElasticStress
    block = '187'
  [../]

  [./elasticity_tensor_188]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '188'
  [../]
  [./strain_188]
    type = ComputeSmallStrain
    block = '188'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_188]
    type = ComputeLinearElasticStress
    block = '188'
  [../]

  [./elasticity_tensor_189]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '189'
  [../]
  [./strain_189]
    type = ComputeSmallStrain
    block = '189'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_189]
    type = ComputeLinearElasticStress
    block = '189'
  [../]

  [./elasticity_tensor_190]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '190'
  [../]
  [./strain_190]
    type = ComputeSmallStrain
    block = '190'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_190]
    type = ComputeLinearElasticStress
    block = '190'
  [../]

  [./elasticity_tensor_191]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 38.0
    euler_angle_2 = 39.0
    euler_angle_3 = 59.0
    block = '191'
  [../]
  [./strain_191]
    type = ComputeSmallStrain
    block = '191'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_191]
    type = ComputeLinearElasticStress
    block = '191'
  [../]

  [./elasticity_tensor_192]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '192'
  [../]
  [./strain_192]
    type = ComputeSmallStrain
    block = '192'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_192]
    type = ComputeLinearElasticStress
    block = '192'
  [../]

  [./elasticity_tensor_193]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '193'
  [../]
  [./strain_193]
    type = ComputeSmallStrain
    block = '193'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_193]
    type = ComputeLinearElasticStress
    block = '193'
  [../]

  [./elasticity_tensor_194]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '194'
  [../]
  [./strain_194]
    type = ComputeSmallStrain
    block = '194'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_194]
    type = ComputeLinearElasticStress
    block = '194'
  [../]

  [./elasticity_tensor_195]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '195'
  [../]
  [./strain_195]
    type = ComputeSmallStrain
    block = '195'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_195]
    type = ComputeLinearElasticStress
    block = '195'
  [../]

  [./elasticity_tensor_196]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '196'
  [../]
  [./strain_196]
    type = ComputeSmallStrain
    block = '196'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_196]
    type = ComputeLinearElasticStress
    block = '196'
  [../]

  [./elasticity_tensor_197]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '197'
  [../]
  [./strain_197]
    type = ComputeSmallStrain
    block = '197'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_197]
    type = ComputeLinearElasticStress
    block = '197'
  [../]

  [./elasticity_tensor_198]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '198'
  [../]
  [./strain_198]
    type = ComputeSmallStrain
    block = '198'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_198]
    type = ComputeLinearElasticStress
    block = '198'
  [../]

  [./elasticity_tensor_199]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '199'
  [../]
  [./strain_199]
    type = ComputeSmallStrain
    block = '199'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_199]
    type = ComputeLinearElasticStress
    block = '199'
  [../]

  [./elasticity_tensor_200]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '200'
  [../]
  [./strain_200]
    type = ComputeSmallStrain
    block = '200'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_200]
    type = ComputeLinearElasticStress
    block = '200'
  [../]

  [./elasticity_tensor_201]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '201'
  [../]
  [./strain_201]
    type = ComputeSmallStrain
    block = '201'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_201]
    type = ComputeLinearElasticStress
    block = '201'
  [../]

  [./elasticity_tensor_202]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '202'
  [../]
  [./strain_202]
    type = ComputeSmallStrain
    block = '202'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_202]
    type = ComputeLinearElasticStress
    block = '202'
  [../]

  [./elasticity_tensor_203]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '203'
  [../]
  [./strain_203]
    type = ComputeSmallStrain
    block = '203'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_203]
    type = ComputeLinearElasticStress
    block = '203'
  [../]


  [./elasticity_tensor_204]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 76.0
    euler_angle_2 = 48.0
    euler_angle_3 = -29.0
    block = '204'
  [../]
  [./strain_204]
    type = ComputeSmallStrain
    block = '204'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_204]
    type = ComputeLinearElasticStress
    block = '204'
  [../]

  [./elasticity_tensor_205]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = -21.0
    euler_angle_2 = 208.0
    euler_angle_3 = 24.0
    block = '205'
  [../]
  [./strain_205]
    type = ComputeSmallStrain
    block = '205'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_205]
    type = ComputeLinearElasticStress
    block = '205'
  [../]

  [./elasticity_tensor_206]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '206'
  [../]
  [./strain_206]
    type = ComputeSmallStrain
    block = '206'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_206]
    type = ComputeLinearElasticStress
    block = '206'
  [../]

  [./elasticity_tensor_207]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '207'
  [../]
  [./strain_207]
    type = ComputeSmallStrain
    block = '207'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_207]
    type = ComputeLinearElasticStress
    block = '207'
  [../]

  [./elasticity_tensor_208]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '208'
  [../]
  [./strain_208]
    type = ComputeSmallStrain
    block = '208'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_208]
    type = ComputeLinearElasticStress
    block = '208'
  [../]

  [./elasticity_tensor_209]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '209'
  [../]
  [./strain_209]
    type = ComputeSmallStrain
    block = '209'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_209]
    type = ComputeLinearElasticStress
    block = '209'
  [../]

  [./elasticity_tensor_210]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '210'
  [../]
  [./strain_210]
    type = ComputeSmallStrain
    block = '210'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_210]
    type = ComputeLinearElasticStress
    block = '210'
  [../]

  [./elasticity_tensor_211]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '211'
  [../]
  [./strain_211]
    type = ComputeSmallStrain
    block = '211'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_211]
    type = ComputeLinearElasticStress
    block = '211'
  [../]

  [./elasticity_tensor_212]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '212'
  [../]
  [./strain_212]
    type = ComputeSmallStrain
    block = '212'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_212]
    type = ComputeLinearElasticStress
    block = '212'
  [../]

  [./elasticity_tensor_213]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '213'
  [../]
  [./strain_213]
    type = ComputeSmallStrain
    block = '213'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_213]
    type = ComputeLinearElasticStress
    block = '213'
  [../]

  [./elasticity_tensor_214]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '214'
  [../]
  [./strain_214]
    type = ComputeSmallStrain
    block = '214'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_214]
    type = ComputeLinearElasticStress
    block = '214'
  [../]

  [./elasticity_tensor_215]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '215'
  [../]
  [./strain_215]
    type = ComputeSmallStrain
    block = '215'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_215]
    type = ComputeLinearElasticStress
    block = '215'
  [../]

  [./elasticity_tensor_216]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '216'
  [../]
  [./strain_216]
    type = ComputeSmallStrain
    block = '216'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_216]
    type = ComputeLinearElasticStress
    block = '216'
  [../]

  [./elasticity_tensor_217]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '217'
  [../]
  [./strain_217]
    type = ComputeSmallStrain
    block = '217'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_217]
    type = ComputeLinearElasticStress
    block = '217'
  [../]

  [./elasticity_tensor_218]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '218'
  [../]
  [./strain_218]
    type = ComputeSmallStrain
    block = '218'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_218]
    type = ComputeLinearElasticStress
    block = '218'
  [../]

  [./elasticity_tensor_219]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '219'
  [../]
  [./strain_219]
    type = ComputeSmallStrain
    block = '219'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_219]
    type = ComputeLinearElasticStress
    block = '219'
  [../]

  [./elasticity_tensor_220]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '220'
  [../]
  [./strain_220]
    type = ComputeSmallStrain
    block = '220'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_220]
    type = ComputeLinearElasticStress
    block = '220'
  [../]

  [./elasticity_tensor_221]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '221'
  [../]
  [./strain_221]
    type = ComputeSmallStrain
    block = '221'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_221]
    type = ComputeLinearElasticStress
    block = '221'
  [../]

  [./elasticity_tensor_222]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '222'
  [../]
  [./strain_222]
    type = ComputeSmallStrain
    block = '222'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_222]
    type = ComputeLinearElasticStress
    block = '222'
  [../]

  [./elasticity_tensor_223]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '223'
  [../]
  [./strain_223]
    type = ComputeSmallStrain
    block = '223'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_223]
    type = ComputeLinearElasticStress
    block = '223'
  [../]

  [./elasticity_tensor_224]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '224'
  [../]
  [./strain_224]
    type = ComputeSmallStrain
    block = '224'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_224]
    type = ComputeLinearElasticStress
    block = '224'
  [../]

  [./elasticity_tensor_225]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '225'
  [../]
  [./strain_225]
    type = ComputeSmallStrain
    block = '225'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_225]
    type = ComputeLinearElasticStress
    block = '225'
  [../]

  [./elasticity_tensor_226]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '226'
  [../]
  [./strain_226]
    type = ComputeSmallStrain
    block = '226'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_226]
    type = ComputeLinearElasticStress
    block = '226'
  [../]

  [./elasticity_tensor_227]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '227'
  [../]
  [./strain_227]
    type = ComputeSmallStrain
    block = '227'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_227]
    type = ComputeLinearElasticStress
    block = '227'
  [../]

  [./elasticity_tensor_228]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '228'
  [../]
  [./strain_228]
    type = ComputeSmallStrain
    block = '228'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_228]
    type = ComputeLinearElasticStress
    block = '228'
  [../]

  [./elasticity_tensor_229]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '229'
  [../]
  [./strain_229]
    type = ComputeSmallStrain
    block = '229'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_229]
    type = ComputeLinearElasticStress
    block = '229'
  [../]

  [./elasticity_tensor_230]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '230'
  [../]
  [./strain_230]
    type = ComputeSmallStrain
    block = '230'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_230]
    type = ComputeLinearElasticStress
    block = '230'
  [../]

  [./elasticity_tensor_231]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '231'
  [../]
  [./strain_231]
    type = ComputeSmallStrain
    block = '231'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_231]
    type = ComputeLinearElasticStress
    block = '231'
  [../]

  [./elasticity_tensor_232]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '232'
  [../]
  [./strain_232]
    type = ComputeSmallStrain
    block = '232'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_232]
    type = ComputeLinearElasticStress
    block = '232'
  [../]

  [./elasticity_tensor_233]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '233'
  [../]
  [./strain_233]
    type = ComputeSmallStrain
    block = '233'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_233]
    type = ComputeLinearElasticStress
    block = '233'
  [../]

  [./elasticity_tensor_234]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '234'
  [../]
  [./strain_234]
    type = ComputeSmallStrain
    block = '234'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_234]
    type = ComputeLinearElasticStress
    block = '234'
  [../]

  [./elasticity_tensor_235]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '235'
  [../]
  [./strain_235]
    type = ComputeSmallStrain
    block = '235'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_235]
    type = ComputeLinearElasticStress
    block = '235'
  [../]

  [./elasticity_tensor_236]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '236'
  [../]
  [./strain_236]
    type = ComputeSmallStrain
    block = '236'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_236]
    type = ComputeLinearElasticStress
    block = '236'
  [../]

  [./elasticity_tensor_237]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '237'
  [../]
  [./strain_237]
    type = ComputeSmallStrain
    block = '237'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_237]
    type = ComputeLinearElasticStress
    block = '237'
  [../]

  [./elasticity_tensor_238]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '238'
  [../]
  [./strain_238]
    type = ComputeSmallStrain
    block = '238'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_238]
    type = ComputeLinearElasticStress
    block = '238'
  [../]

  [./elasticity_tensor_239]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '239'
  [../]
  [./strain_239]
    type = ComputeSmallStrain
    block = '239'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_239]
    type = ComputeLinearElasticStress
    block = '239'
  [../]

  [./elasticity_tensor_240]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '240'
  [../]
  [./strain_240]
    type = ComputeSmallStrain
    block = '240'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_240]
    type = ComputeLinearElasticStress
    block = '240'
  [../]

  [./elasticity_tensor_241]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '241'
  [../]
  [./strain_241]
    type = ComputeSmallStrain
    block = '241'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_241]
    type = ComputeLinearElasticStress
    block = '241'
  [../]

  [./elasticity_tensor_242]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '242'
  [../]
  [./strain_242]
    type = ComputeSmallStrain
    block = '242'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_242]
    type = ComputeLinearElasticStress
    block = '242'
  [../]

  [./elasticity_tensor_243]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    block = '243'
  [../]
  [./strain_243]
    type = ComputeSmallStrain
    block = '243'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_243]
    type = ComputeLinearElasticStress
    block = '243'
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
    file_base = out_484grain_structure
    elemental_as_nodal = true
  [../]
[]
