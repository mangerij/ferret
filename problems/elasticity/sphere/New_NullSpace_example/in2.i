# Patch Test

# This tests calculates the stress in a sheet with a cylindrical hole, with
# uniforn stress boundary conditions sigma_yy applied to the edges perpendicular to
# the y-axis. There exists an analytical solution for an infinite sheet with
# a cylindrical hole in it.

#

[Mesh]#Comment
  file = sphere_id7_od13_meshed_v2.e
#  uniform_refine = 1
  displacements = 'disp_x disp_y disp_z'
[] # Mesh

[Problem]
  dimNullSpace     = 6
[]
#
[UserObjects]
  [./RigidModes3DNullSpace]
     type=RigidBodyModes3D
     variable = disp_x
     subspace_name = NullSpace
     subspace_indices = '0 1 2 3 4 5 '
     modes = 'trans_x trans_y trans_z rot_x rot_y rot_z'
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
  [./elastic_energy]
    order = CONSTANT
    family = MONOMIAL
  [../]

[] #AuxVariables

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
    index_j=1
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 2
    index_j=2
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
    index_i = 3
    index_j = 1
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
    index_i = 3
    index_j = 1
  [../]

  [./elastic_energy]
    type = TensorElasticEnergyAux
    variable = elastic_energy
  [../]


[] # AuxKernels

[BCs]
#
  [./outsidex]
    type = HydrostaticBC
    variable = disp_x
    boundary = '1001'
    pressure = 1.e-08
    component = 0
  [../]
  [./outsidey]
    type = HydrostaticBC
    variable = disp_y
    boundary = '1001'
    pressure = 1.e-08
    component = 1
  [../]  
  [./outsidez]
    type = HydrostaticBC
    variable = disp_z
    boundary = '1001'
    pressure = 1.e-08
    component = 2
  [../]
#  [./surface_elasticity_x]
#    type = SurfaceMechanicsBC
#    variable = 'disp_x'
#    disp_x = disp_x
#    disp_y = disp_y
#    disp_z = disp_z
#    boundary = '1001'
#    surface_euler_angle_1 = '0.0'
#    surface_euler_angle_2 = '0.0'
#    surface_euler_angle_3 = '0.0'
## Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
#    Cs_ijkl = '200e-09 0.e-09 0.e-09 0.e-09 0.e-09 0.e-09 0.e-09 0.e-09 200e-09'
## Intrinsic surface stress
#    taus_ij = '0.0 0.0 0.0 0.0'
#    component = 0
#  [../]
#  [./surface_elasticity_y]
#    type = SurfaceMechanicsBC
#   variable = 'disp_y'
#   disp_x = disp_x
#   disp_y = disp_y
#   disp_z = disp_z 
#   boundary = '1001'
#    surface_euler_angle_1 = '0.0'
#    surface_euler_angle_2 = '0.0'
#    surface_euler_angle_3 = '0.0'
## Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
#    Cs_ijkl = '200e-09 0.e-09 0.e-09 0.e-09 0.e-09 0.e-09 0.e-09 0.e-09 200e-09'
## Intrinsic surface stress
#    taus_ij = '0.0 0.0 0.0 0.0'
#    component = 1
#3  [../]
#  [./surface_elasticity_z]
#    type = SurfaceMechanicsBC
#    variable = 'disp_z'
#    boundary = '1001'
#    disp_x = disp_x
#    disp_y = disp_y
#    disp_z = disp_z
#    surface_euler_angle_1 = '0.0'
#    surface_euler_angle_2 = '0.0'
#    surface_euler_angle_3 = '0.0'
## Surface elastic tensor C_1111, C_1112, C_1122, C_1212, C_1222, C_1211, C_2211, C_2212, C_2222
#    Cs_ijkl = '200e-09 0.e-09 0.e-09 0.e-09 0.e-09 0.e-09 0.e-09 0.e-09 200e-09'
## Intrinsic surface stress
#    taus_ij = '0.0 0.0 0.0 0.0'
#    component = 2
#  [../]

#
#  [./outsidex]
#    type = HydrostaticDirichletBC
#    variable = disp_x
#    boundary = '1001'
#    center_x = 0.
#    center_y = 0.
#    center_z = 0.
#    displacement = -0.1
#    component = 0
#  [../]
#  [./outsidey]
#    type = HydrostaticDirichletBC
#    variable = disp_y
#    boundary = '1001'
#    center_x = 0.
#    center_y = 0.
#    center_z = 0.
#    displacement = -0.1
#    component = 1
#  [../]  
#  [./outsidez]
#    type = HydrostaticDirichletBC
#    variable = disp_z
#    boundary = '1001'
#    center_x = 0.
#    center_y = 0.
#    center_z = 0.
#    displacement = -0.1
#    component = 2 
#  [../]

[] # BCs

[Materials]

  [./Goo1]
    type = LinearElasticMaterial
    block = '101'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./Goo2]
    type = LinearElasticMaterial
    block = '102'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./Goo3]
    type = LinearElasticMaterial
    block = '103'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./Goo4]
    type = LinearElasticMaterial
    block = '104'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]
  
  [./Goo5]
    type = LinearElasticMaterial
    block = '105'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./Goo6]
    type = LinearElasticMaterial
    block = '106'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./Goo7]
    type = LinearElasticMaterial
    block = '107'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./Goo8]
    type = LinearElasticMaterial
    block = '108'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./Goo9]
    type = LinearElasticMaterial
    block = '109'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./Goo10]
    type = LinearElasticMaterial
    block = '110'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./Goo11]
    type = LinearElasticMaterial
    block = '111'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./Goo12]
    type = LinearElasticMaterial
    block = '112'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./Goo13]
    type = LinearElasticMaterial
    block = '113'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./Goo14]
    type = LinearElasticMaterial
    block = '114'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./Goo15]
    type = LinearElasticMaterial
    block = '115'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]
  
  [./Goo0]
    type = LinearElasticMaterial
    block = '11'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
# cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
#    C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
# isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='266e-09 266e-09 266e-09 266e-09 266e-09 266e-09 60e-09 60e-09 60e-09'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

[] # Materials

[Preconditioning]
  active = 'smp'
#  active = 'FDP'
  [./smp]
  type = SMP
  full = true
#  [../]
#  [./FDP]
#   type = FDP
#  [../]
[]
[Postprocessors]
  [./integrated_elastic_energy]
    type = ElementIntegralVariablePostprocessor
    variable = elastic_energy
    block = 11
  [../]
  [./volume]
    type = VolumePostprocessor
    use_displaced_mesh = true
#    variable = disp_x
    block = 11
  [../]
[]
[Executioner]

   type = Steady
#  petsc_options = '-ksp_monitor'
#  petsc_options_iname = '-ksp_type -pc_type'
#  petsc_options_value = 'gmres lu'
#  petsc_options = `snes snes_view -ksp_view -snes_monitor -ksp_monitor -pc_asm_print_subdomains'
#  petsc_options_iname = `-ksp_type -pc_type -pc_asm_decomposition -pc_asm_sub_pc_sype`
#  petsc_options_value = `gmres asm block lu'
#  type = Transient
  petsc_options = '-snes_monitor -snes_converged_reason -ksp_monitor -ksp_converged_reason'
  petsc_options_iname = '-ksp_type -pc_type -pc_asm_overlap -ksp_max_it -snes_max_it'
  petsc_options_value = '    gmres     asm  16              100        3'
#  petsc_options_iname = '-ksp_type -pc_type '
#  petsc_options_value = '    gmres      svd'
#  petsc_options = '-ksp_monitor -ksp_view -snes_view'
#  petsc_options_iname = '-ksp_type -pc_type -pc_asm_overlap -sub_pc_type'
#  petsc_options_value = 'gmres asm 8 lu'
#  petsc_options_iname = 'ksp_type -pc_type'
#  petsc_options_value = 'gmres asm'
#  petsc_options_value = 'gmres ilu'
#  petsc_options_iname = '-ksp_type -pc_type'
#  petsc_options_value = 'gmres lu'

  nl_abs_tol = 1e-10
#  l_rel_tol = 1e-8
#  l_abs_tol  = 1e-10

  l_max_its = 50

#  start_time = 0.0
#  dt = 1.0
#  num_steps = 2
#  end_time = 2.0
[] # Executioner

[Output]
  file_base = out_sphere_amorphous_HBC_lu
  interval = 1
  output_initial = true
  elemental_as_nodal = true
  exodus = true
  tecplot = true
  perf_log = true
[] # Output
