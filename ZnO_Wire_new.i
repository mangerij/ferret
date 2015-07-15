[Mesh]
  file = wire_100.e
#  uniform_refine = 1
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
  [./EgZnO]
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
  ### Here we use the results from Phys Rev B 88, 235210 (2013) Wagner et al
  # which computes the strain-induced bandgap change of wurtzite ZnO using
  # the HSE approach. Here E0, db, du, Rb and nu are material properties.
  # Units will be in eV. The biaxial stress in the wire may be different than the core-shell case.
  [./bandgap]
    type = BandGapAuxZnO
    variable = EgZnO
    relaxed_energy = 3.200
    uniaxial_strain_rate = -3.800
    biaxial_strain_rate = -0.450
    biaxial_relaxation_coeff = 0.929
    poisson_ratio = 0.31
  [../]
[] # AuxKernels


[BCs]
  active = 'back_x_bc back_y_bc back_z_bc top_z_bc'
  [./back_x_bc]
    type = DirichletBC
    variable = disp_x
    boundary = 3
     value = 0
   [../]
   [./back_y_bc]
     type = DirichletBC
     variable = disp_y
     boundary = 3
     value = 0
   [../]
  [./back_z_bc]
    type = DirichletBC
    variable = disp_z
    boundary = 3
    value = 0
  [../]

  #Load boundary condition (-.1 um displacement along z)

  #  [./top_x_bc]
  #   type = DirichletBC
  #    variable = disp_x
  #    boundary = 2
  #    value = 0
  #  [../]
  #
  #  [./top_y_bc]
  #    type = DirichletBC
  #   variable = disp_y
  #    boundary = 2
  #    value = 0
  #  [../]
  [./top_z_bc]
    type = DirichletBC
    variable = disp_z
    boundary = 2
    value = -66                     #Edit this value#
 [../]

 #-------------Surface--------#

 [./surface_elasticity_X_surf1]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_x
    boundary = '1'
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
    boundary = '1'
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
    boundary = '1'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
  # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
  # Intrinsic surface stress
    taus = '-1.7e-09'
    component = 2
  [../]

 [./surface_elasticity_X_surf2]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_x
    boundary = '2'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
  # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
  # Intrinsic surface stress
    taus = '-1.7e-09'
    component = 0
  [../]

  [./surface_elasticity_Y_surf2]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_y
    boundary = '2'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
  # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
  # Intrinsic surface stress
    taus = '-1.7e-09'
    component = 1
  [../]
  [./surface_elasticity_Z_surf2]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_z
    boundary = '2'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
  # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
  # Intrinsic surface stress
    taus = '-1.7e-09'
    component = 2
  [../]
  [./surface_elasticity_X_surf3]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_x
    boundary = '3'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
  # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
  # Intrinsic surface stress
    taus = '-1.7e-09'
    component = 0
  [../]
  [./surface_elasticity_Y_surf3]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_y
    boundary = '3'
    surface_euler_angle_1 = 0.0
    surface_euler_angle_2 = 0.0
    surface_euler_angle_3 = 0.0
  # Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
    Cs_ijkl = '49.1e-09 0.0e-09 15.1e-09 13.7e-09 0.0e-09 0.0e-09 15.1e-09 0.0e-09 34.9e-09'
  # Intrinsic surface stress
    taus = '-1.7e-09'
    component = 1
  [../]
  [./surface_elasticity_Z_surf3]
    type = SurfaceMechanicsBC
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    variable = disp_z
    boundary = '3'
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
  #active = 'cube'
  #[./cube]
  #  type = LinearElasticMaterial
  #  block = '1'
  #  disp_x = disp_x
  #  disp_y = disp_y
  #  disp_z = disp_z
  ## C11 C12 C13 C22 C23 C33 C44 C55 C66
  #  C_ijkl = '209.7e-09 121.1e-09 105.1e-09 209.7e-09 105.1e-09 210.9e-09 42.47e-09 42.47e-09 44.29e-09'
  #  euler_angle_1 = 0.0
  #  euler_angle_2 = 0.0
  #  euler_angle_3 = 0.0
  #[../]
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
    block = '1'
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = '1'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '1'
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
    petsc_options_iname = '-ksp_gmres_restart -snes_rtol -ksp_rtol -pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_zeropivot -pc_factor_zeropivot -pc_hypre_type'
    petsc_options_value = '    675              1e-8      1e-12      hypre        5               lu              1e-50                    1e-50     boomeramg'
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'       #"PJNK, JFNK, NEWTON"
[]

[Debug]
  show_var_residual_norms = true
[]

[Outputs]
  file_base = wire100test
  #output_initial = true
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    elemental_as_nodal = true
  [../]
[]
