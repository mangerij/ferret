
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 40
  ny = 40
  nz = 12
  ymax = 20
  ymin = -20
  xmin = -20
  xmax = 20
  zmin = -6
  zmax = 6
[]


[GlobalParams]
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  displacements = 'disp_x disp_y disp_z'
  #use_displaced_mesh = false
  C_ijkl = '178 96.4 96.4 178 96.4 178 122 122 122'
  prefactor = -0.01
[]



[Variables]
  [./disp_x]
    block = '0'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./disp_y]
    block = '0'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./disp_z]
    block = '0'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
[]

[Kernels]
  #Elastic problem
  [./TensorMechanics]
  #This is an action block
  [../]
[]

[BCs]
#  [./top_disp_x]
#    variable = disp_x
#    type = DirichletBC
#    value = 0
#    boundary = 'front'
#  [../]
#  [./top_disp_y]
#    variable = disp_y
#    type = DirichletBC
#    value = 0
#    boundary = 'front'
#  [../]
#  [./top_disp_z]
#    variable = disp_z
#    type = DirichletBC
#    value = 0
#    boundary = 'front'
#  [../]

  #[./bot_disp_x]
  #  variable = disp_x
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'back'
  #[../]
  #[./bot_disp_y]
  #  variable = disp_y
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'back'
  #[../]
  #[./bot_disp_z]
  #  variable = disp_z
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'back'
  #[../]
  #
  #[./top_disp_x]
  #  variable = disp_x
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'top'
  #[../]
  #[./top_disp_y]
  #  variable = disp_y
  #  type = DirichletBC
  #  value = -0.5
  #  boundary = 'top'
  #[../]
  #[./top_disp_z]
  #  variable = disp_z
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'top'
  #[../]
  #
  #[./bottom_disp_x]
  #  variable = disp_x
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'bottom'
  #[../]
  #[./bottom_disp_y]
  #  variable = disp_y
  #  type = DirichletBC
  #  value = 0.5
  #  boundary = 'bottom'
  #[../]
  #[./bottom_disp_z]
  #  variable = disp_z
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'bottom'
  #[../]
  #
  #
  #
  #
  #[./right_disp_x]
  #  variable = disp_x
  #  type = DirichletBC
  #  value = -0.5
  #  boundary = 'right'
  #[../]
  #[./right_disp_y]
  #  variable = disp_y
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'right'
  #[../]
  #[./right_disp_z]
  #  variable = disp_z
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'right'
  #[../]
  #
  #[./left_disp_x]
  #  variable = disp_x
  #  type = DirichletBC
  #  value = 0.5
  #  boundary = 'left'
  #[../]
  #[./left_disp_y]
  #  variable = disp_y
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'left'
  #[../]
  #[./left_disp_z]
  #  variable = disp_z
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'left'
  #[../]
[]

[Materials]
#   [./eigen_strain_xx_yy] #Use for stress-free strain (ie epitaxial)
#    type = ComputeEigenstrain
#    block = '0'
#   # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
#    eigen_base = '1 0 0 0 1 0 0 0 0'
#  [../]
  [./BTO]
    type=LinearElasticMaterial
    block = '0'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = '0'
    fill_method = symmetric9
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = '0'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '0'
  [../]
[]

[Postprocessors]
  [./elastic_energy]
    type = ElasticEnergy
    execute_on = 'timestep_end'
    block = '0'
  [../]
  [./volume]
    type = VolumePostprocessor
    use_displaced_mesh = true
  [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -options_left'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type    -sub_pc_factor_zeropivot -pc_factor_zeropivot -pc_side '
    petsc_options_value = '    121                1e-8      1e-8    gamg        1e-50    1e-50      left        '
  [../]
[]

#Limits exist on -snes_rtol =< 1e-10.

[Executioner]
  type = Steady
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_island_comp
    output_initial = true
    elemental_as_nodal = true
  [../]
[]
