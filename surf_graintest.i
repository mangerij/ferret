[Mesh]
  file = 6grains.e
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
    C_ijkl = '309 195.6 195.6 309 309 539 349.53 349.53 349.53'
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


 [./surface_elasticity_X_surf3]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_x
    boundary = '24'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
# Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
# Intrinsic surface stress
    taus = '-1.7e-06'
    component = 0
  [../]

  [./surface_elasticity_Y_surf3]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_y
    boundary = '24'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
# Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
# Intrinsic surface stress
    taus = '-1.7e-06'
    component = 1
  [../]

  [./surface_elasticity_Z_surf3]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_z
    boundary = '24'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
# Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
# Intrinsic surface stress
    taus = '-1.7e-06'
    component = 2
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
  solve_type = 'PJFNK'       #"PJFNK, JFNK, NEWTON"
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_grains_surf_elastic
    elemental_as_nodal = true
  [../]
[]


