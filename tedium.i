[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 2
  ny = 2
  nz = 2
  ymax = 1
  ymin = 0
  xmin = 0
  xmax = 1
  zmin = 0
  zmax = 1
[]

[GlobalParams]
   #length scale
   len_scale = 1.0
   #BulkEnergy coefficients
  # alpha1 = -1.8202e8 # 3.8(T-479)*10^5 C^{-2}m^2 (T = 0 K)
  # alpha11 = -7.3e7
  # alpha111 = 2.6e8
  # alpha12 = 7.5e8
  # alpha112 = 6.1e8
  # alpha123 = -3.7e9
  # #Wall stuff
  # G110 = 0.6e-10
  # G11/G110 = 0.6
  # G12/G110 = 0.0
  # G44/G110 = 0.3
  # G44P/G110 = 0.3
   permittivity = 8.854187e-12
   polar_x = polar_x
   polar_y = polar_y
   polar_z = polar_z
   potential_int = potential_int
   potential_ext = potential_ext
[]

[Variables]
  #[./polar_x]
  #  order = FIRST
  #  family = LAGRANGE
  #  block='0'
  #[../]
  #[./polar_y]
  #  order = FIRST
  #  family = LAGRANGE
  #  block='0'
  #[../]
  #[./polar_z]
  #  order = FIRST
  #  family = LAGRANGE
  #  block='0'
  #[../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
    block='0'
  [../]
  #[./potential_ext]
  #  order = FIRST
  #  family = LAGRANGE
  #  block='0'
  #[../]
[]

[Kernels]
  #Bulk energy density
  #[./polar_electric_E]
  #   type=PolarElectricEStrong
  #   variable=potential_int
  #   block='0'
  #[../]
  [./FE_E_int]
     type=Electrostatics
     variable=potential_int
     block='0'
  [../]
  #[./FE_E_ext]
  #   type=Electrostatics
  #   variable=potential_ext
  #   block='0'
  #[../]
  #[./polar_electric_px]
  #   type=PolarElectricPStrong
  #   variable = polar_x
  #   component = 0
  #   block='0'
  #[../]
  #[./polar_electric_py]
  #   type=PolarElectricPStrong
  #   variable = polar_y
  #   component = 1
  #   block='0'
  #[../]
  #[./polar_electric_pz]
  #   type=PolarElectricPStrong
  #   variable = polar_z
  #   component = 2
  #   block='0'
  #[../]
[]

#[ICs]
#  [./polar_x_randic]
#     type=RandomIC
#     variable=polar_x
#     min = 0.6
#     max = 0.75
#  [../]
#  [./polar_y_randic]
#     type=RandomIC
#     variable=polar_y
#     min = 0.6
#     max = 0.75
#  [../]
#  [./polar_z_randic]
#     type=RandomIC
#     variable=polar_z
#     min = 0.1
#     max = 0.15
#  [../]
#  #[./disp_x_randic]
#  #   type=RandomIC
#  #   variable=disp_x
#  #   min = 1
#  #   max = 1.1
#  #[../]
#  #[./disp_y_randic]
#  #   type=RandomIC
#  #   variable=disp_y
#  #   min = 1
#  #   max = 1.1
#  #[../]
#  #[./disp_z_randic]
#  #   type=RandomIC
#  #   variable=disp_z
#  #   min = 1
#  #   max = 1.1
#  #[../]
#[]

[BCs]
  [./top_potential]
    type = DirichletBC
    value = 1.0
    boundary = 'top'
     variable = potential_int
  [../]
  [./bottom_potential]
    type = DirichletBC
    value = 0.0
    boundary = 'bottom'
    variable = potential_int
  [../]
[]


#[BCs]
#   [./potential_int_upz]
#     type = DirichletBC
#     variable = potential_int
#     boundary = '1'
#     value = 0.0
#   [../]
#   [./potential_int_downz]
#     type = DirichletBC
#     variable = potential_int
#     boundary = '2'
#     value = 0.0
#   [../]
#
#   #Displacements [5 6 7 8] [x y -x -y] boundaries, [3 4] [top bottom]
#  # [./disp_x_slab3]
#  #   type = PresetBC
#  #   variable = disp_x
#  #   boundary = '3'
#  #   value = 0.0
#  # [../]
#  # [./disp_y_slab3]
#  #   type = PresetBC
#  #   variable = disp_y
#  #   boundary = '3'
#  #   value = 0.0
#  # [../]
#  # [./disp_z_slab3]
#  #   type = PresetBC
#  #   variable = disp_z
#  #   boundary = '3'
#  #   value = 0.0
#  # [../]
#   #
#  # [./disp_x_slab4]
#  #   type = PresetBC
#  #   variable = disp_x
#  #   boundary = '4'
#  #   value = 0.0
#  # [../]
#  # [./disp_y_slab4]
#  #   type = PresetBC
#  #   variable = disp_y
#  #   boundary = '4'
#  #   value = 0.0
#  # [../]
#  # [./disp_z_slab4]
#  #   type = PresetBC
#  #   variable = disp_z
#  #   boundary = '4'
#  #   value = 0.0
#  # [../]
#
#  # [./disp_x_slab5]
#  #   type = DirichletBC
#  #   variable = disp_x
#  #   boundary = '5'
#  #  # value = 1.0e-1
#  #   value = 0.5
#  # [../]
#   [./disp_y_slab5]
#     type = DirichletBC
#     variable = disp_y
#     boundary = '5'
#     value = 0.5
#   [../]
#  # [./disp_z_slab5]
#  #   type = PresetBC
#  #   variable = disp_z
#  #   boundary = '5'
#  #  # value = 1.0e-1
#  #   value = 0.0
#  # [../]
#
#  # [./disp_x_slab6]
#  #   type = PresetBC
#  #   variable = disp_x
#  #   boundary = '6'
#  #  # value = 1.0e-1
#  #  value = 0.0
#  # [../]
#  # [./disp_y_slab6]
#  #   type = PresetBC
#  #   variable = disp_y
#  #   boundary = '6'
#  #  # value = 1.0e-1
#  #  value = 0.0
#  # [../]
#  # [./disp_z_slab6]
#  #   type = PresetBC
#  #   variable = disp_z
#  #   boundary = '6'
#  #  # value = 1.0e-1
#  #  value = 0.0
#  # [../]
#
#  # [./disp_x_slab7]
#  #   type = DirichletBC
#  #   variable = disp_x
#  #   boundary = '7'
#  #  # value = -1.0e-1
#  #   value = -0.5
#  # [../]
#   [./disp_y_slab7]
#     type = DirichletBC
#     variable = disp_y
#     boundary = '7'
#     value = -0.5
#   [../]
#  # [./disp_z_slab7]
#  #   type = PresetBC
#  #   variable = disp_z
#  #   boundary = '7'
#  #  # value = -1.0e-1
#  #   value = 0.0
#  # [../]
#
# #  [./disp_x_slab8]
# #    type = PresetBC
# #    variable = disp_x
# #    boundary = '8'
# #   # value = -1.0e-1
# #   value = 0.0
# #  [../]
# #  [./disp_y_slab8]
# #    type = PresetBC
# #    variable = disp_y
# #    boundary = '8'
# #   # value = -1.0e-1
# #   value = 0.0
# #  [../]
# #  [./disp_z_slab8]
# #    type = PresetBC
# #    variable = disp_z
# #    boundary = '8'
# #   # value = -1.0e-1
# #   value = 0.0
# #[../]
#  [./potential_ext_upz]
#    type = DirichletBC
#    variable = potential_ext
#    boundary = '1'
#    value = 0.0 #this is near-zero
#  [../]
#  [./potential_ext_downz]
#    type = DirichletBC
#    variable = potential_ext
#    boundary = '2'
#    value = 0.0
#  [../]
#[]

#[Materials]
#  [./slab_ferroelectric]
#    type=LinearFerroelectricMaterial
#    block = '2'
#    disp_x = disp_x
#    disp_y = disp_y
#    disp_z = disp_z
#    #in GPA. from N. Pandech et al. Ceramic. Internat., times 1e9 to convert to N/m^2
#    # C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl = '380.0e9 150.0e9 150.0e9 380.0e9 150.0e9 380.0e9 110.0e9 110.0e9 110.0e9'
#    #in m^4/C^2. from http://arxiv.org/pdf/1205.5640.pdf
#    # Q11 Q12 Q13 Q22 Q23 Q33 Q44 Q55 Q66
#    Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.034 0.034 0.034'
#    euler_angle_1 = 0.0 #currently will only rotate C_ijkl
#    euler_angle_2 = 0.0
#    euler_angle_3 = 0.0
#  [../]
#  #[./slab_elastic]
#  #  type=LinearElasticMaterial
#  #  block = '2'
#  #  disp_x = disp_x
#  #  disp_y = disp_y
#  #  disp_z = disp_z
#  #  #in GPA. from N. Pandech et al. Ceramic. Internat., times 1e9 to convert to N/m^2
#  #  # C11 C12 C13 C22 C23 C33 C44 C55 C66
#  #  C_ijkl = '380.0e9 150.0e9 150.0e9 380.0e9 150.0e9 380.0e9 110.0e9 110.0e9 110.0e9'
#  #  euler_angle_1 = 0.0 #currently will only rotate C_ijkl
#  #  euler_angle_2 = 0.0
#  #  euler_angle_3 = 0.0
#  #[../]
#  [./elasticity_tensor]
#    type = ComputeElasticityTensor
#    block = '2'
#    C_ijkl = '380.0e9 150.0e9 150.0e9 380.0e9 150.0e9 380.0e9 110.0e9 110.0e9 110.0e9'
#    fill_method = symmetric9
#  [../]
#
#  [./vacuum]
#    type=GenericConstantMaterial
#    block = '1'
#  [../]
#[]

[Executioner]
  type=Steady
  nl_max_its = 350
  #[./TimeStepper]
  #  type = IterationAdaptiveDT
  #  dt = 1.15e-15
  #  optimal_iterations = 10
  #  growth_factor = 1.01
  #  cutback_factor =  0.99
  #[../]
  full = 'true'
  solve_type = 'newton'
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1.0e-25
  dtmax = 1.81e15
  num_steps = 15000
  dt = 1.15e-15
  #splitting = 'ferretsplit'
  petsc_options = '-snes_test_display'
  petsc_options_iname = '-snes_type'
  petsc_options_value = 'test'
  #petsc_options ='-snes_linesearch_monitor -options_left '
  #petsc_options_iname = '-snes_rtol -ksp_rtol -pc_type  -sub_pc_type -pc_asm_overlap -sub_pc_factor_zeropivot -pc_factor_zeropivot'
  #petsc_options_value = '  1e-8       1e-14      asm         lu            2               1e-50                    1e-50'
[]

#
#[Splits]
#  [./ferretsplit]
#    type = Split
#    splitting = 'ferroelectric elastic' #split to two subproblems
#    splitting_type = schur #schur split somewhat bugged right now (only serial)
#    schur_type = full
#    schur_pre = A11
#    petsc_options ='-snes_linesearch_monitor'
#  [../]
#
#  [./ferroelectric]
#    vars = 'polar_x polar_y polar_z potential_int potential_ext'
#    petsc_options='-dm_view ' #'-ksp_monitor -inner_ksp_monitor'
#    petsc_options_iname='-ksp_type   -ksp_gmres_restart  -ksp_rtol -inner_pc_type -pc_type   -pc_asm_overlap -inner_pc_factor_zeropivot  -pc_factor_zeropivot'
#    petsc_options_value=' gmres            350              1e-14     hypre         hypre         5         1e-50                1e-50  '
#  [../]
#  [./elastic]
#    vars = 'disp_x disp_y disp_z'
#    petsc_options='-dm_view'  # ''-ksp_monitor'
#    petsc_options_iname='-ksp_type  -ksp_gmres_restart -inner_pc_type -ksp_rtol -pc_type   -pc_asm_overlap -inner_pc_factor_zeropivot -pc_factor_zeropivot'
#    petsc_options_value = 'gmres       350                 hypre         1e-14     hypre          5           1e-50                        1e-50'
#  [../]
#[]


[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_PBP_test
    output_initial = true
    elemental_as_nodal = false
    interval = 1
  [../]
  #[./debug]
  #  type = VariableResidualNormsDebugOutput
  #[../]
[]
