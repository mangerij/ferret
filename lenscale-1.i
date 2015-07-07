[Mesh]
  type=GeneratedMesh
  dim=3
  nx=2
  ny=2
  nz=2
  xmin=0.0
  xmax=10.0
  ymin=0.0
  ymax=10.0
  zmin=0.0
  zmax=10.0
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
    block = '0'
    len_scale = 1.0e-9
  [../]
  [./stressdiv_1]
    type = StressDivergenceTensorsScaled
    variable = disp_y
    component = 1
    block = '0'
    len_scale = 1.0e-9
  [../]
  [./stressdiv_2]
    type = StressDivergenceTensorsScaled
    variable = disp_z
    component = 2
    block = '0'
    len_scale = 1.0e-9
  [../]
[]

#[BCs]
#  [./back_x_bc]
#    type = DirichletBC
#    variable = disp_x
#    boundary = 'top'
#    value = 0
#  [../]
#  [./back_y_bc]
#    type = DirichletBC
#    variable = disp_y
#    boundary = 'top'
#    value = 0
#  [../]
#  [./back_z_bc]
#    type = DirichletBC
#    variable = disp_z
#    boundary = 'top'
#    value = 0
#  [../]
#
#
#  [./front_x_bc]
#    type = DirichletBC
#    variable = disp_x
#    boundary = 'bottom'
#    value = 0.0
#  [../]
#  [./front_y_bc]
#    type = DirichletBC
#    variable = disp_y
#    boundary = 'bottom'
#    value = 0.0
#  [../]
#  [./front_z_bc]
#    type = DirichletBC
#    variable = disp_z
#    boundary = 'bottom'
#    value = 0.1
#  [../]
#[]

[Materials]
  active = 'cube'
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
  solve_type = 'NEWTON'
  #petsc_options = '-snes_test_display'
  petsc_options_iname = '-snes_type'
  petsc_options_value = 'test'
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_elasticity_scaled
  [../]
[]
