
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
  len_scale = 1.0
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
[]

[Materials]
  [./eigen_strain_zz] #Use for stress-free strain (ie epitaxial)
    type = ComputeEigenstrain
    block = '1'
    # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
    eigen_base = '1 0 0 0 1 0 0 0 0'
    eigenstrain_name = eigenstrain
    prefactor = 0.0
  [../]
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #from MaterialsProject
    C_ijkl = '281 115.74 115.74 281 115.74 281 97.18 97.18 97.18'
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

  [./slab_ferroelectric]
    block = '1'
    type = ComputeElectrostrictiveTensor
    Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.03375 0.03375 0.03375'
    C_ijkl = '281 115.74 115.74 281 115.74 281 97.18 97.18 97.18'
  [../]
[]




[Kernels]
  #Elastic problem
  [./TensorMechanics]
  #This is an action block
  [../]
[]


[BCs]
  [./center_disp_x]
    type = DirichletBC
    variable = 'disp_x'
    value = 0.0
    boundary = 'center_node_1 center_node_2 center_node_3 center_node_4'
  [../]
  [./center_disp_y]
    type = DirichletBC
    variable = 'polar_y'
    value = 0.0
    boundary = 'center_node_1 center_node_2 center_node_3 center_node_4'
  [../]



  [./center_disp_z]
    type = DirichletBC
    variable = 'disp_z'
    value = -0.01
    boundary = '2'
  [../]

  [./center_disp_z]
    type = DirichletBC
    variable = 'disp_z'
    value = 0.01
    boundary = '3'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    250              1e-10      1e-8      1e-6      bjacobi   '
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    execute_on = 'timestep_end'
    file_base = out_disk
    elemental_as_nodal = true
  [../]
[]
