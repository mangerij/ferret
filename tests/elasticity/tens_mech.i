[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 4
  ny = 4
  nz = 4
  ymax = 2
  ymin = -2
  xmin = -2
  xmax = 2
  zmin = -2
  zmax = 2
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Kernels]
  [./TensorMechanics]
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
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

[BCs]
  [./back_x_bc]
    type = DirichletBC
    variable = disp_x
    boundary = 'top'
     value = 0.03
  [../]
  [./back_x_bc]
    type = DirichletBC
    variable = disp_x
    boundary = 'bottom'
     value = -0.03
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '209.7e-09 121.1e-09 105.1e-09 209.7e-09 105.1e-09 210.9e-09 42.47e-09 42.47e-09 44.29e-09'
    fill_method = symmetric9
    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    block = '0'
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = '0'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '0'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-info -snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -options_left'
    petsc_options_iname = '-ksp_gmres_restart -snes_rtol -ksp_rtol -pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_zeropivot -pc_factor_zeropivot '
    petsc_options_value = '    675              1e-8      1e-10      asm        1               lu              1e-50                    1e-50    '
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'       #"PJNK, JFNK, NEWTON"
[]


[Outputs]
  file_base = out_test_tens_mech
  #output_initial = true
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    elemental_as_nodal = true
  [../]
[]
