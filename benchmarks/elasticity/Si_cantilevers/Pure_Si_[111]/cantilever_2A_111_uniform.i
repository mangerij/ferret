# This is an example for a cantilever bending calculation

[Mesh]
  file = cantilever_2A.e
  uniform_refine = 1
  displacements = 'disp_x disp_y disp_z'
[] # Mesh

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

[] # Variables

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

[] # AuxVariables

[TensorMechanics]
  [./solid]
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
[]

[AuxKernels]

  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 1
    index_j = 1
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 2
    index_j = 2
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 3
    index_j = 3
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 1
    index_j = 2
  [../]
  [./stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 2
    index_j = 3
  [../]
  [./stress_zx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zx
    index_i = 1
    index_j = 3
  [../]
  [./strain_xx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_xx
    index_i = 1
    index_j = 1
  [../]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_yy
    index_i = 2
    index_j = 2
  [../]
  [./strain_zz]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zz
    index_i = 3
    index_j = 3
  [../]
  [./strain_xy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_xy
    index_i = 1
    index_j = 2
  [../]
  [./strain_yz]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_yz
    index_i = 2
    index_j = 3
  [../]
  [./strain_zx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zx
    index_i = 1
    index_j = 3
  [../]

[] # AuxKernels


[BCs]

  active = 'anchor3_X anchor3_Y anchor3_Z anchor5_Z'

  #  Cantilever side with the (-1,0,0) normal; this one is clamped: 

  [./anchor3_X]
    type = DirichletBC
    variable = disp_x
    boundary = '3'
    value = 0.0
  [../]

  [./anchor3_Y]
    type = DirichletBC
    variable = disp_y
    boundary = '3'
    value = 0.0
  [../]

  [./anchor3_Z]
    type = DirichletBC
    variable = disp_z
    boundary = '3'
    value = 0.0
  [../]

  #  Cantilever side with the (1,0,0) normal; this one is pulled down [um]: 

  [./anchor5_X]
    type = DirichletBC
    variable = disp_x
    boundary = '5'
    value = 0.0
  [../]

  [./anchor5_Y]
    type = DirichletBC
    variable = disp_y
    boundary = '5'
    value = 0.0
  [../]

  [./anchor5_Z]
    type = DirichletBC
    variable = disp_z
    boundary = '5'
    value = -0.1
  [../]

[] # BCs

[Materials]

  [./grain1]
    type = LinearElasticMaterial
    block = '1'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# Elastic constants for silicon (cubic) in GPa at room T
# C_11 = 166.0 GPa, C_12 = 64.0 GPa, C_44 = 79.5 GPa.
    C_ijkl = '166.0 64.0 64.0 166.0 64.0 166.0 79.5 79.5 79.5'
# This setup orients one of the cubic axes along (111)
    euler_angle_1 =   0.0
    euler_angle_2 =  54.73561
    euler_angle_3 =  45.0
  [../]

[] # Materials

[Preconditioning]
   active = 'smp'
   [./smp]
   type = SMP
   full = true
   [../]
[]

[Executioner]

  type = Steady
#  petsc_options = '-ksp_monitor'
#  petsc_options_iname = '-ksp_type -pc_type'
#  petsc_options_value = 'gmres lu'

#  type = Transient
  petsc_options = '-snes_mf_operator -snes_monitor -snes_converged_reason -ksp_monitor -ksp_converged_reason -ksp_view -snes_view'
  petsc_options_iname = '-ksp_type -pc_type -pc_asm_overlap -sub_pc_type'
  petsc_options_value = 'gmres asm 20 lu'
##  petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -pc_hypre_type'
##  petsc_options_value = 'gmres hypre hypre boomeramg'

#  petsc_options = '-snes_mf_operator -ksp_monitor -ksp_view -snes_view -pc_svd_monitor'
#  petsc_options_iname = '-ksp_type -pc_type'
#  petsc_options_value = '    gmres      svd'

  nl_abs_tol = 1e-12
  #l_abs_tol  = 1e-5

  l_max_its = 150

#  start_time = 0.0
#  dt = 1.0
#  num_steps = 2
#  end_time = 2.0
[] # Executioner

[Output]
  file_base = out_cantilever_2A_uniform_asm20_meshref1
#  file_base = out_cantilever_2A_uniform_asm
  interval = 1
  output_initial = true
  elemental_as_nodal = true
  exodus = true
  perf_log = true
[] # Output
