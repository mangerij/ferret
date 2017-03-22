
[Mesh]
  file = exodus_disk_r8_h1_photo.e
[]

[MeshModifiers]
  [./centernodeset_1]
    type = AddExtraNodeset
    new_boundary = 'center_node_1'
    coord = '0.0 0.0 -0.5'
  [../]
  [./centernodeset_2]
    type = AddExtraNodeset
    new_boundary = 'center_node_2'
    coord = '0.0 0.0 0.5'
  [../]
  [./centernodeset_3]
    type = AddExtraNodeset
    new_boundary = 'center_node_3'
    coord = '0.0 0.0 0.166667'
  [../]
  [./centernodeset_4]
    type = AddExtraNodeset
    new_boundary = 'center_node_4'
    coord = '0.0 0.0 -0.166667'
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

  [./beta_11_impermeability]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./beta_22_impermeability]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./beta_12_impermeability]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./beta_13_impermeability]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./beta_33_impermeability]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./beta_23_impermeability]
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
    block  = '1 2 3 4'
  [../]
  [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    variable = strain_xy_elastic
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
  [./matl_e13]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 2
    variable = strain_xz_elastic
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    variable = strain_yy_elastic
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
  [./matl_e23]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 2
    variable = strain_yz_elastic
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
  [./matl_e33]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 2
    index_j = 2
    variable = strain_zz_elastic
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
  [./matl_s11]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = stress_xx_elastic
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
  [./matl_s12]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = stress_xy_elastic
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
  [./matl_s13]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
    variable = stress_xz_elastic
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
 [./matl_s22]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = stress_yy_elastic
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
  [./matl_s23]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    variable = stress_yz_elastic
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
  [./matl_s33]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    variable = stress_zz_elastic
    block  = '1 2 3 4'
    execute_on = 'timestep_end'
  [../]

  [./beta_e11]
    type = Impermeability
    index_i = 0
    index_j = 0
    variable = beta_11_impermeability
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
  [./beta_e12]
    type = Impermeability
    index_i = 0
    index_j = 1
    variable = beta_12_impermeability
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
  [./beta_e13]
    type = Impermeability
    index_i = 0
    index_j = 2
    variable = beta_13_impermeability
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
  [./beta_e22]
    type = Impermeability
    index_i = 1
    index_j = 1
    variable = beta_22_impermeability
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
  [./beta_e23]
    type = Impermeability    index_i = 1
    index_j = 2
    variable = beta_23_impermeability
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
  [./beta_e33]
    type = Impermeability
    index_i = 2
    index_j = 2
    variable = beta_33_impermeability
    execute_on = 'timestep_end'
    block  = '1 2 3 4'
  [../]
[]


[Materials]
  [./eigen_strain_zz] #Use for stress-free strain (ie epitaxial)
    type = ComputeEigenstrain
    block = '1 2 3 4'
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
    block = '2 3'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    block = '2 3'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
    block = '2 3'
  [../]
  [./photoelastic_tensor_1]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric9 #symmetric21
    #C1111, C1122, C1133, C1123, C1113, C1112, C2222, C2233, C2223, C2213, C2212, C3333, C3323, C3313, C3312, C2323, C2313, C2312, C1313, C1312, C1212
    P_mnkl = '0.5 0.106 0.20 0.50 0.106 0.5 0.77 0.77 0.77'
    block = '2 3'
  [../]
  [./delta_beta_tensor_1]
    type = ComputeDeltaBetaTensor
    block = '2 3'
  [../]
  [./beta_tensor_1]
    type = ComputeBetaTensor
    n_o = 2.437
    n_e = 2.365
    block = '2 3'
  [../]

  [./elasticity_tensor_2]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 45.0
    euler_angle_2 = 37.0
    euler_angle_3 = -76.0 #ne = -2.162e-3, no = -4.279e-4, bf = -1.734e-3
    block = '1 4'
  [../]
  [./strain_2]
    type = ComputeSmallStrain
    block = '1 4'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_2]
    type = ComputeLinearElasticStress
    block = '1 4'
  [../]
  [./photoelastic_tensor_2]
    type = ComputePhotostrictiveTensor
    fill_method = symmetric9 #symmetric21
    #C1111, C1122, C1133, C1123, C1113, C1112, C2222, C2233,

 C2223, C2213, C2212, C3333, C3323, C3313, C3312, C2323, C2313, C2312, C1313, C1312, C1212
    P_mnkl = '0.5 0.106 0.20 0.0 0.0 0.0 0.5 0.20 0.0 #FIX



0.50 0.106 0.5 0.77 0.77 0.77'
    block = '1 4'
  [../]
  [./delta_beta_tensor_2]
    type = ComputeDeltaBetaTensor
    block = '1 4'
  [../]
  [./beta_tensor_2]
    type = ComputeBetaTensor
    n_o = 2.437
    n_e = 2.365
    block = '1 4'
    euler_angle_1 = 45.0
    euler_angle_2 = 37.0
    euler_angle_3 = -76.0
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
    value = -0.0025
    boundary = '5'
  [../]

  [./center_disp_z_bottom]
    type = DirichletBC
    variable = 'disp_z'
    value = 0.0025
    boundary = '6'
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
    file_base = out_disk
    elemental_as_nodal = true
  [../]
[]
