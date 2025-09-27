
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 20
  ny = 10
  nz = 10
  xmin = -10
  xmax = 10
  ymin = -5
  ymax = 5
  zmin = -5
  zmax = 5
  elem_type = HEX8
[]

[GlobalParams]
  potential_E_int = potential_E_int
  disp_x = u_x
  disp_y = u_y
  disp_z = u_z
  displacements = 'u_x u_y u_z'
[]

[Variables]
  [./u_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./u_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./u_z]
    order = FIRST
    family = LAGRANGE
  [../]
  [./potential_E_int]
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
    variable = u_x
    component = 0
  [../]
  [./piezocouple_1]
    type = ConversePiezoelectricStrain
    variable = u_y
    component = 1
  [../]
  [./piezocouple_2]
    type = ConversePiezoelectricStrain
    variable = u_z
    component = 2
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_E_int
  [../]
  [./strain_charge]
     type = PiezoelectricStrainCharge
     variable = potential_E_int
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
    type = ComputePiezostrictiveTensor
    fill_method = general
    e_ijk = '0 0 -0.00415 0 0 0 -0.00415 0 0 0 0 0 0 0 -0.00415 0 -0.00415 0 -0.005 0 0 0 -0.005 0 0 0 0.0124'
  [../]
  [./permitivitty_1]
    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '0.0721616'
  [../]
[]

[Functions]
  [./bc_func_1]
    type = ParsedFunction
    expression = -100.0*sin(0.05*t)
  [../]
[]

[BCs]
  # Boundary Condition System
  [./back_pot]
    type = DirichletBC
    variable = potential_E_int
    boundary = 'back'
    value = 0.0
  [../]
  [./front_pot]
    type = FunctionDirichletBC
    variable = potential_E_int
    boundary = 'front'
    function = bc_func_1
  [../]

  [./fixed_end_x]
    type = DirichletBC
    variable = u_x
    boundary = 'left'
    value = 0.0
  [../]
  [./fixed_end_y]
    type = DirichletBC
    variable = u_y
    boundary = 'left'
    value = 0.0
  [../]
  [./fixed_end_z]
    type = DirichletBC
    variable = u_z
    boundary = 'left'
    value = 0.0
  [../]
[]

[Postprocessors]
  [./sxx]
    type = ElementAverageValue
    variable = stress_xx
    execute_on = 'initial timestep_end final'
  [../]
  [./sxz]
    type = ElementAverageValue
    variable = stress_zx
    execute_on = 'initial timestep_end final'
  [../]
  [./szz]
    type = ElementAverageValue
    variable = stress_zz
    execute_on = 'initial timestep_end final'
  [../]
  [./Felastic]
    type = ElasticEnergy
    execute_on = 'timestep_end'
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
  type = Transient
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"

  [./TimeIntegrator]
    type = NewmarkBeta
  [../]

  dtmin = 1e-10
  dtmax = 5.0

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 25  #usually 10
    linear_iteration_ratio = 100
    dt = 0.1
    growth_factor = 1.5
  [../]


  num_steps = 200

[]


[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_piezoelectric
    elemental_as_nodal = true
  [../]
[]
