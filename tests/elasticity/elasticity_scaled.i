[Mesh]
  type=GeneratedMesh
  dim=3
  nx=5
  ny=5
  nz=5
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  zmin=0.0
  zmax=1.0
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

[Kernels]
  [./stressdiv_0]
    type = StressDivergenceTensorsScaled
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_x
    component = 0
    len_scale = 1e-9
  [../]
  [./stressdiv_1]
    type = StressDivergenceTensorsScaled
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_y
    component = 1
    len_scale = 1e-9
  [../]
  [./stressdiv_2]
    type = StressDivergenceTensorsScaled
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_z
    component = 2
    len_scale = 1e-9
  [../]
[]

[BCs]
  [./back_x_bc]
    type = DirichletBC
    variable = disp_x
    boundary = 'top'
    value = 0
  [../]
  [./back_y_bc]
    type = DirichletBC
    variable = disp_y
    boundary = 'top'
    value = 0
  [../]
  [./back_z_bc]
    type = DirichletBC
    variable = disp_z
    boundary = 'top'
    value = 0
  [../]


  [./front_x_bc]
    type = DirichletBC
    variable = disp_x
    boundary = 'bottom'
    value = 0.0
  [../]
  [./front_y_bc]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0.0
  [../]
  [./front_z_bc]
    type = DirichletBC
    variable = disp_z
    boundary = 'bottom'
    value = 0.1
  [../]
[]

[Materials]
  [./cube]
    type = LinearElasticMaterial
    block = '0'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
# C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '209.7e9 121.1e9 105.1e9 209.7e9 105.1e9 210.9e9 42.47e9 42.47e9 44.29e9'
    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]
[]

[Executioner]
  type = Steady
  petsc_options = '-snes_converged_reason -ksp_converged_reason -ksp_view -snes_view'
  petsc_options_iname = '-ksp_type -snes_type   -pc_type '
  petsc_options_value = 'gmres      newtonls      hypre'
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_elasticity_scaled
  [../]
[]
