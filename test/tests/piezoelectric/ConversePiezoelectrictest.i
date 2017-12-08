
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 7
  ny = 7
  nz = 5
  xmin = -3
  xmax = 3
  ymin = -3
  ymax = 3
  zmin = -2
  zmax = 2
  elem_type = HEX8
[]

[GlobalParams]
  potential_int = potential_int
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
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
  [./potential_int]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
    [./stress_xx]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./stress_yy]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./stress_zz]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./stress_xy]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./stress_yz]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./stress_zx]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./strain_xx]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./strain_yy]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./strain_zz]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./strain_xy]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./strain_yz]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./strain_zx]
      order = CONSTANT
      family = MONOMIAL
    [../]
[]

[Kernels]
  #Elastic problem
  [./TensorMechanics]
  #This is an action block
  [../]
  [./piezocouple_0]
    type = ConversePiezoelectricStrain
    variable = disp_x
    component = 0
  [../]
  [./piezocouple_1]
    type = ConversePiezoelectricStrain
    variable = disp_y
    component = 1
  [../]
  [./piezocouple_2]
    type = ConversePiezoelectricStrain
    variable = disp_z
    component = 2
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     permittivity = 0.0721616
  [../]
  [./strain_charge]
     type = PiezoelectricStrainCharge
     variable = potential_int
  [../]
[]


[AuxKernels]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  [../]
  [./stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 1
    index_j = 2
  [../]
  [./stress_zx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zx
    index_i = 2
    index_j = 0
  [../]
  [./strain_xx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_xx
    index_i = 0
    index_j = 0
  [../]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
  [../]
  [./strain_zz]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zz
    index_i = 2
    index_j = 2
  [../]
  [./strain_xy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_xy
    index_i = 0
    index_j = 1
  [../]
  [./strain_yz]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_yz
    index_i = 1
    index_j = 2
  [../]
  [./strain_zx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zx
    index_i = 2
    index_j = 0
  [../]
[]

[Materials]
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '209.7 121.1 105.1 209.7 105.1 210.9 42.47 42.47 44.29'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
  [../]
  [./d333]
    type = ComputePiezoTensor
    fill_method2 = symmetric9
    fill_method = general
    compute_piezostrictive_coeff = true
    C_ijkl = '209.7 121.1 105.1 209.7 105.1 210.9 42.47 42.47 44.29'
    d_kpq = '0 0 -0.00415 0 0 0 -0.00415 0 0 0 0 0 0 0 -0.00415 0 -0.00415 0 -0.005 0 0 0 -0.005 0 0 0 0.0124'
    d_pqkT = '0 0 -0.005 0 0 0 -0.00415 0 0 0 0 0 0 0 -0.005 0 -0.00415 0 -0.00415 0 0 0 -0.00415 0 0 0 0.0124'
  [../]
[]


[BCs]
  # Boundary Condition System
  [./back_pot]
    type =DirichletBC
    variable = potential_int
    boundary = 'back'
    value = -0.25
  [../]
  [./front_pot]
    type =DirichletBC
    variable = potential_int
    boundary = 'front'
    value = 0.25
  [../]

  #[./stablizer_x]
  #  type = DirichletBC
  #  variable = 'disp_x'
  #  boundary = 'front'
  #  value = 0.0
  #[../]
  #[./stablizer_y]
  #  type = DirichletBC
  #  variable = 'disp_y'
  #  boundary = 'front'
  #  value = 0.0
  #[../]
  #[./stablizer_z]
  #  type = DirichletBC
  #  variable = 'disp_z'
  #  boundary = 'front'
  #  value = 0.0
  [../]
[]

[Postprocessors]
  [./Felastic]
    type = ElasticEnergy
    execute_on = 'timestep_end'
  [../]
[]

[Problem]
  null_space_dimension = 6
[]

[UserObjects]
  [./rigidbodymodes_x]
     type = RigidBodyModes3D
     subspace_name = NullSpace
     subspace_indices = '0 1 2 3 4 5'
     modes = 'trans_x trans_y trans_z rot_x rot_y rot_z'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart  -snes_atol -ksp_rtol -pc_type'
    petsc_options_value = '    121                1e-10      1e-8     bjacobi'
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
    file_base = out_test_conversepiezoelectric
    elemental_as_nodal = true
  [../]
[]
