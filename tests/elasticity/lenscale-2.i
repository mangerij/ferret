[Mesh]
  type=GeneratedMesh
  dim=3
  nx=3
  ny=3
  nz=3
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
  #[./TensorMechanics]
  #  disp_x = disp_x
  #  disp_y = disp_y
  #  disp_z = disp_z
  #[../]
  [./stressdiv_0]
    type = StressDivergenceTensorsScaled
    variable = disp_x
    component = 0
    len_scale = 1
  [../]
  [./stressdiv_1]
    type = StressDivergenceTensorsScaled
    variable = disp_y
    component = 1
    len_scale = 1
  [../]
  [./stressdiv_2]
    type = StressDivergenceTensorsScaled
    variable = disp_z
    component = 2
    len_scale = 1
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
  active = 'cube'
  [./cube]
    type = LinearElasticMaterial
    block = '0'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
# C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '209.7e-09 121.1e-09 105.1e-09 209.7e-09 105.1e-09 210.9e-09 42.47e-09 42.47e-09 44.29e-09'
    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options = '-snes_check_jacobian'
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_elasticity_scaled
  [../]
[]
