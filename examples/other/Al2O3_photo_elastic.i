
[Mesh]
  file = exodus_disk_r8_h1_photo.e
  uniform_refine = 1
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

  [./n_o_ij]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./n_e_ij]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dn]
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

  [./ref_deps_00_ij]
    type = RefractiveIndex
    index_one = 0
    index_two = 0
    variable = n_o_ij
    execute_on = 'timestep_end'
  [../]
  [./ref_deps_22_ij]
    type = RefractiveIndex
    index_one = 2
    index_two = 2
    variable = n_e_ij
    execute_on = 'timestep_end'
  [../]

  [./dn_r]
    type = Birefringence
    variable = dn
    n_o = n_o_ij
    n_e = n_e_ij
    execute_on = 'timestep_end'
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
    fill_method = symmetric9
    #S11 S21 S31 S12 S22 S32 S13 S23 S33
    P_mnkl = '0.5 0.106 0.20 0.50 0.106 0.5 0.77 0.77 0.77'
    block = '2 3'
  [../]
  [./beta_tensor_1]
    type = ComputeBetaTensor
    block = '2 3'
  [../]
  [./indicatrix_1]
    type = ComputeIndicatrix
    block = '2 3'
    #from Handbook of Optics Chp. 17
    refractive_index_bulk_ordinary = 2.437
    refractive_index_bulk_extraordinary = 2.365
  [../]


  [./elasticity_tensor_2]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #BaTiO3 from MaterialsProject
    C_ijkl = '260.06 105.79 76.90 260.06 105.79 260.06 81.57 81.57 116.28'
    euler_angle_1 = 0.0 # 45.0
    euler_angle_2 = 54.74 #37.0
    euler_angle_3 = 45.0 #-76.0  #ne = -2.162e-3, no = -4.279e-4, bf = -1.734e-3
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
    fill_method = symmetric9
    P_mnkl = '0.5 0.106 0.20 0.50 0.106 0.5 0.77 0.77 0.77'
    block = '1 4'
  [../]
  [./beta_tensor_2]
    type = ComputeBetaTensor
    block = '1 4'
  [../]
  [./indicatrix_2]
    type = ComputeIndicatrix
    block = '1 4'
    #from BTO paper
    refractive_index_bulk_ordinary = 2.412
    refractive_index_bulk_extraordinary = 2.360
    euler_angle_1 = 0.0 # 45.0
    euler_angle_2 = 54.74 #37.0
    euler_angle_3 = 45.0 #-76.0
  [../]
[]

[Kernels]
  #Elastic problem
  [./SolidMechanics]
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
