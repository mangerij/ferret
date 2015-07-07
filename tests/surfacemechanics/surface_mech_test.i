[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 3
  ny = 3
  nz = 3
  ymax = 2
  ymin = -2
  xmin = -2
  xmax = 2
  zmin = -2
  zmax = 2
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
  [./TensorMechanics]
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


[BCs]
 [./surface_elasticity_X_surf1]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_x
    boundary = 'top'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
# Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
# Intrinsic surface stress
    taus = '-1.7e-09'
    component = 0
  [../]
  [./surface_elasticity_Y_surf1]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_y
    boundary = 'top'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
# Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
# Intrinsic surface stress
    taus = '-1.7e-09'
    component = 1
  [../]
  [./surface_elasticity_Z_surf1]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_z
    boundary = 'top'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
# Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
# Intrinsic surface stress
    taus = '-1.7e-09'
    component = 2
  [../]
  [./surface_elasticity_X_surf1]
     type = SurfaceMechanicsBC
     disp_x = disp_x
     disp_y = disp_y
     disp_z = disp_z
     variable = disp_x
     boundary = 'bottom'
     surface_euler_angle_1 = 0.0
     surface_euler_angle_2 = 0.0
     surface_euler_angle_3 = 0.0
 # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
     Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
 # Intrinsic surface stress
     taus = '-1.7e-09'
     component = 0
   [../]
   [./surface_elasticity_Y_surf1]
     type = SurfaceMechanicsBC
     disp_x = disp_x
     disp_y = disp_y
     disp_z = disp_z
     variable = disp_y
     boundary = 'bottom'
     surface_euler_angle_1 = 0.0
     surface_euler_angle_2 = 0.0
     surface_euler_angle_3 = 0.0
 # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
     Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
 # Intrinsic surface stress
     taus = '-1.7e-09'
     component = 1
   [../]
   [./surface_elasticity_Z_surf1]
     type = SurfaceMechanicsBC
     disp_x = disp_x
     disp_y = disp_y
     disp_z = disp_z
     variable = disp_z
     boundary = 'bottom'
     surface_euler_angle_1 = 0.0
     surface_euler_angle_2 = 0.0
     surface_euler_angle_3 = 0.0
 # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
     Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
 # Intrinsic surface stress
     taus = '-1.7e-09'
     component = 2
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
    C_ijkl = '209.7e-09 121.1e-09 105.1e-09 209.7e-09 105.1e-09 210.9e-09 42.47e-09 42.47e-09 44.29e-09'
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
    file_base = out_ZnO_NW_NotScaled
  [../]
[]
