[Mesh]
  file = 6grains.e
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
    block = '1 2 3'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    block = '1 2 3'
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
    block = '1 2 3'
  [../]

  [./elasticity_tensor_2]
    type = ComputeElasticityTensor
    C_ijkl = '319 99.6 99.6 319 99.6 319 109.53 109.53 109.53'
    fill_method = symmetric9
    block = '4 5'
    euler_angle_1 = 0.0
    euler_angle_2 = 83.0
    euler_angle_3 = 0.0
  [../]
  [./strain_2]
    type = ComputeSmallStrain
    block = '4 5'
  [../]
  [./stress_2]
    type = ComputeLinearElasticStress
    block = '4 5'
  [../]

  [./elasticity_tensor_3]
    type = ComputeElasticityTensor
    C_ijkl = '539 295.6 295.6 539 295.6 539 349.53 349.53 349.53'
    fill_method = symmetric9
    block = '6'
    euler_angle_1 = 55.0
    euler_angle_2 = 26.0
    euler_angle_3 = -10.0
  [../]
  [./strain_3]
    type = ComputeSmallStrain
    block = '6'
  [../]
  [./stress_3]
    type = ComputeLinearElasticStress
    block = '6'
  [../]
[]

[Kernels]
  [./TensorMechanics]
  [../]
[]



[BCs]
   active = 'anchor_up_Z anchor_dn_Z anchor_up_X anchor_dn_X anchor_up_Y anchor_dn_Y'
  [./anchor_up_X]
    type = DirichletBC
    variable = disp_x
    boundary = '2 13 15 21 27 28'
    value = 0.0
  [../]

  [./anchor_up_Y]
    type = DirichletBC
    variable = disp_y
    boundary = '2 13 15 21 27 28'
    value = 0.0
  [../]

  [./anchor_up_Z]
    type = DirichletBC
    variable = disp_z
    boundary = '2 13 15 21 27 28'
    value = 0.025
  [../]
 
  [./anchor_dn_X]
    type = DirichletBC
    variable = disp_x
    boundary = '6 9 18 22 24 29'
    value = 0.0
  [../]

  [./anchor_dn_Y]
    type = DirichletBC
    variable = disp_y
    boundary = '6 9 18 22 24 29'
    value = 0.0
  [../]

  [./anchor_dn_Z]
    type = DirichletBC
    variable = disp_z
    boundary = '6 9 18 22 24 29'
    value = -0.025
  [../]

[]


[Postprocessors]
    [./Felastic]
      type = ElasticEnergy
      block = '1 2 3 4 5 6'
      execute_on = 'timestep_end'
    [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -ksp_snes_ew'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type -pc_hypre_type'
    petsc_options_value = '    121                1e-8      1e-8    hypre    boomeramg'
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
[]

[Outputs]
  print_linear_residuals = true
  [./out]
    type = Exodus
    file_base = out_grains_elastic
    elemental_as_nodal = true
  [../]
[]


