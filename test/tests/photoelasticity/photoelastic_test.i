
[Mesh]
  file = 6grains.e
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
    index_j = 0
    index_k = 0
    var1 = dn_1
    execute_on = 'timestep_end'
  [../]

  [./n_2_c]
    type = RefractiveIndex
    variable = n_2
    index_j = 1
    index_k = 1
    var1 = dn_2
    execute_on = 'timestep_end'
  [../]

  [./n_3_c]
    type = RefractiveIndex
    variable = n_3
    index_j = 2
    index_k = 2
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
    block = '1 2 3 4 5 6'
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
    type = ComputeElastoopticTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 0.0
    euler_angle_2 = 15.0
    euler_angle_3 = 65.0
    block = '1'
  [../]
  [./beta_tensor_1]
    type = ComputeIndicatrix
    block = '1'
    euler_angle_1 = 0.0
    euler_angle_2 = 15.0
    euler_angle_3 = 65.0
  [../]
  [./delta_beta_tensor_1]
    type = ComputeDeltaIndicatrix
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
    type = ComputeElastoopticTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 12.0
    euler_angle_2 = -15.0
    euler_angle_3 = 65.0
    block = '2'
  [../]
  [./beta_tensor_2]
    type = ComputeIndicatrix
    block = '2'
    euler_angle_1 = 12.0
    euler_angle_2 = -15.0
    euler_angle_3 = 65.0
  [../]
  [./delta_beta_tensor_2]
    type = ComputeDeltaIndicatrix
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
    type = ComputeElastoopticTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 1.0
    euler_angle_2 = 50.0
    euler_angle_3 = -132.0
    block = '3'
  [../]
  [./beta_tensor_3]
    type = ComputeIndicatrix
    block = '3'
    euler_angle_1 = 1.0
    euler_angle_2 = 50.0
    euler_angle_3 = -132.0
  [../]
  [./delta_beta_tensor_3]
    type = ComputeDeltaIndicatrix
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
    type = ComputeElastoopticTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 100.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
    block = '4'
  [../]
  [./beta_tensor_4]
    type = ComputeIndicatrix
    block = '4'
    euler_angle_1 = 100.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]
  [./delta_beta_tensor_4]
    type = ComputeDeltaIndicatrix
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
    type = ComputeElastoopticTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 200.0
    euler_angle_2 = 33.0
    euler_angle_3 = 5.0
    block = '5'
  [../]
  [./beta_tensor_5]
    type = ComputeIndicatrix
    block = '5'
    euler_angle_1 = 200.0
    euler_angle_2 = 33.0
    euler_angle_3 = 5.0
  [../]
  [./delta_beta_tensor_5]
    type = ComputeDeltaIndicatrix
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
    type = ComputeElastoopticTensor
    fill_method = symmetric21 #BTO is not symmetric21 FIX!! 
    # Use BaTiO3, crystal symmetry P4mm.
    P_mnkl = '0.5 0.106 0.2 0.0 0.0 0.0 0.5 0.2 0.0 0.0 0.0 0.77 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.1'
    euler_angle_1 = 110.0
    euler_angle_2 = 0.0
    euler_angle_3 = 100.0
    block = '6'
  [../]
  [./beta_tensor_6]
    type = ComputeIndicatrix
    block = '6'
    euler_angle_1 = 110.0
    euler_angle_2 = 0.0
    euler_angle_3 = 100.0
  [../]
  [./delta_beta_tensor_6]
    type = ComputeDeltaIndicatrix
    block = '6'
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
    value = -0.1
    boundary = 'side1'
  [../]

  [./center_disp_z_bottom]
    type = DirichletBC
    variable = 'disp_z'
    value = 0.0
    boundary = 'side2'
  [../]

  [./center_disp_x_bottom]
    type = DirichletBC
    variable = 'disp_x'
    value = 0.0
    boundary = 'side2'
  [../]

  [./center_disp_y_bottom]
    type = DirichletBC
    variable = 'disp_y'
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
    file_base = out_6grain_photo
    elemental_as_nodal = true
  [../]
[]
