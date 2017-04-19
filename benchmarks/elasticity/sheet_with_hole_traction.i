# Patch Test

# This tests calculates the stress in a sheet with a cylindrical hole, with
# uniforn stress boundary conditions sigma_yy applied to the edges perpendicular to
# the y-axis. There exists an analytical solution for an infinite sheet with
# a cylindrical hole in it.

#

# With Young's modulus at 1e6 and Poisson's ratio at 0, the shear
#  modulus is 5e5 (G=E/2/(1+nu)).  Therefore,
#
#  stress xx = 1e6 * 1e-6 = 1
#  stress yy = 1e6 * 2e-6 = 2
#  stress zz = 1e6 * 3e-6 = 3
#  stress xy = 2 * 5e5 * 1e-6 / 2 = 0.5
#             (2 * G   * gamma_xy / 2 = 2 * G * epsilon_xy)
#  stress yz = 2 * 5e5 * 2e-6 / 2 = 1
#  stress zx = 2 * 5e5 * 3e-6 / 2 = 1.5


[Mesh]
  file = sheet_with_hole_thick.e
[]

[GlobalParams]
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

[] #AuxVariables

[Kernels]
  [./TensorMechanics]
  #This is an action block
  [../]
[]

#[Functions]
#  [./f_disp_x]
#      type=SolutionFunction
#      file_type=exodusII
#      mesh=in.e
#      variable=disp_x
#      timestep=2
#  [../]
#  [./f_disp_y]
#      type=SolutionFunction
#      file_type=exodusII
#      mesh=in.e
#      variable=disp_y
#      timestep=2
#  [../]
#  [./f_disp_z]
#      type=SolutionFunction
#      file_type = exodusII
#      mesh = in.e
#      variable=disp_z
#      timestep=2
#  [../]
#  [./f_stress_xx]
#      type=SolutionFunction
#      file_type=exodusII
#      mesh=in.e
#      variable=stress_xx
#      timestep=2
#  [../]
#  [./f_stress_xy]
#      type=SolutionFunction
#      file_type=exodusII
#      mesh=in.e
#      variable=stress_xy
#      timestep=2
#  [../]
#  [./f_stress_yy]
#      type=SolutionFunction
#      file_type=exodusII
#      mesh=in.e
#      variable=stress_yy
#      timestep=2
#  [../]
#  [./f_stress_yz]
#      type=SolutionFunction
#      file_type=exodusII
#      mesh=in.e
#      variable=stress_yz
#      timestep=2
#  [../]
#  [./f_stress_zx]
#      type=SolutionFunction
#      file_type=exodusII
#      mesh=in.e
#      variable=stress_zx
#      timestep=2
#  [../]
#  [./f_stress_zz]
#      type=SolutionFunction
#      file_type=exodusII
#      mesh=in.e
#      variable=stress_zz
#      timestep=2
#  [../]
#[]


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

#hacks
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
# End hack


[] # AuxKernels

[BCs]

  #[./anchor_up_X_pos]
  #  type = StressFunctionBC
  #  variable = disp_x
  #  boundary = '2'
  #  component=0
  #  stress_xx = f_stress_xx
  #  stress_xy = f_stress_xy
  #  stress_yy = f_stress_yy
  #  stress_yz = f_stress_yz
  #  stress_zx = f_stress_zx
  #  stress_zz = f_stress_zz
  #[../]
  #
  #[./anchor_up_Y_pos]
  #  type = StressFunctionBC
  #  variable = disp_y
  #  boundary = '2'
  #  component = 1
  #  stress_xx = f_stress_xx
  #  stress_xy = f_stress_xy
  #  stress_yy = f_stress_yy
  #  stress_yz = f_stress_yz
  #  stress_zx = f_stress_zx
  #  stress_zz = f_stress_zz
  #[../]
  #
  #[./anchor_up_Z_pos]
  #  type = StressFunctionBC
  #  variable = disp_z
  #  boundary = '2'
  #  component = 2
  #  stress_xx = f_stress_xx
  #  stress_xy = f_stress_xy
  #  stress_yy = f_stress_yy
  #  stress_yz = f_stress_yz
  #  stress_zx = f_stress_zx
  #  stress_zz = f_stress_zz
  #[../]

  [./anchor_up_X_neg]
    type = DirichletBC
    variable = disp_x
    boundary = '1'
    value = 0.0
  [../]

  [./anchor_up_Z_neg]
    type = DirichletBC
    variable = disp_z
    boundary = '1'
    value = 0.0
  [../]


  [./anchor_up_Y_neg]
    type = DirichletBC
    variable = disp_y
    boundary = '1'
    value = 1.e-6
  [../]


[] # BCs

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
    # cubic symmetry C11 C12 C13 C22 C23 C33 C44 C55 C66
    # C_ijkl ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 0.5e6 1.5e6 1.5e6 1.5e6'
    # isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl = '3.e6 1.e6 1.e6 3.e6 1.e6 3.e6 1.e6 1.e6 1.e6'
    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
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
[] # Materials

[Preconditioning]
  active = 'smp'
  [./smp]
  type = SMP
  full = true
  [../]
#   type = FDP
[]
[Executioner]

   type = Steady
#  petsc_options = '-ksp_monitor'
   petsc_options_iname = '-ksp_type -pc_type'
   petsc_options_value = 'gmres lu'
#  petsc_options = `snes snes_view -ksp_view -snes_monitor -ksp_monitor -pc_asm_print_subdomains'
#  petsc_options_iname = `-ksp_type -pc_type -pc_asm_decomposition -pc_asm_sub_pc_sype`
#  petsc_options_value = `gmres asm block lu'
#  type = Transient

  petsc_options = '-snes_monitor -snes_converged_reason -ksp_converged_reason'
#  petsc_options_iname = '-ksp_type -pc_type '
#  petsc_options_value = '    gmres      svd'
#  petsc_options = '-ksp_monitor -ksp_view -snes_view'
#  petsc_options_iname = '-ksp_type -pc_type -pc_asm_overlap -sub_pc_type'
#  petsc_options_value = 'gmres asm 10 lu'
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

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    execute_on = 'timestep_end'
    file_base = out_hole
    elemental_as_nodal = true
  [../]
[]
