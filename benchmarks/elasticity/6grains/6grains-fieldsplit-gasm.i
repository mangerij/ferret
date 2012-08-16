# Patch Test

# This test is designed to compute constant xx, yy, zz, xy, yz, and zx
#  stress on a set of irregular hexes.  The mesh is composed of one
#  block with seven elements.  The elements form a unit cube with one
#  internal element.  There is a nodeset for each exterior node.

# The cube is displaced by 1e-6 units in x, 2e-6 in y, and 3e-6 in z.
#  The faces are sheared as well (1e-6, 2e-6, and 3e-6 for xy, yz, and
#  zx).  This gives a uniform strain/stress state for all six unique
#  tensor components.

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


[Mesh]#Comment
  file = 6grains.e
  displacements = 'disp_x disp_y disp_z'
[] # Mesh

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
  [./elastic_energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hydrostatic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./firstinv]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./secondinv]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./thirdinv]
    order = CONSTANT
    family = MONOMIAL
  [../]

[] # AuxVariables

[SolidMechanics]
  [./solid]
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
[]

[AuxKernels]

  [./stress_xx]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_xx
    index = 0
  [../]
  [./stress_yy]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_yy
    index = 1
  [../]
  [./stress_zz]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_zz
    index = 2
  [../]
  [./stress_xy]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_xy
    index = 3
  [../]
  [./stress_yz]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_yz
    index = 4
  [../]
  [./stress_zx]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_zx
    index = 5
  [../]
  [./elastic_energy]
    type = ElasticEnergyAux
    variable = elastic_energy
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]

[] # AuxKernels

[BCs]

  [./anchor_up_X]
    type = DirichletBC
    variable = disp_x
    boundary = '2 13 15 21 27 28'
    value = 0.0
  [../]

  [./anchor_up_Y]
    type = DirichletBC
    variable = disp_y
    boundary = '2 13 15 21 27 28'
    value = 0.0
  [../]

  [./anchor_up_Z]
    type = DirichletBC
    variable = disp_z
    boundary = '2 13 15 21 27 28'
    value = 2e-6
  [../]
 
  [./anchor_dn_X]
    type = DirichletBC
    variable = disp_x
    boundary = '6 9 18 22 24 29'
    value = 0.0
  [../]

  [./anchor_dn_Y]
    type = DirichletBC
    variable = disp_y
    boundary = '6 9 18 22 24 29'
    value = 0.0
  [../]

  [./anchor_dn_Z]
    type = DirichletBC
    variable = disp_z
    boundary = '6 9 18 22 24 29'
    value = -2e-6
  [../]

[] # BCs

[Materials]

  [./Goo1]
    type = LinearGeneralAnisotropicMaterial
    block = '1'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
    C_matrix ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 6.0e6 0.5e6 0.5e6 1.5e6'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./Goo2]
    type = LinearGeneralAnisotropicMaterial
    block = '2'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
    C_matrix ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 6.0e6 0.5e6 0.5e6 1.5e6'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./Goo3]
    type = LinearGeneralAnisotropicMaterial
    block = '3'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
    C_matrix ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 6.0e6 0.5e6 0.5e6 1.5e6'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]  

  [./Goo4]
    type = LinearGeneralAnisotropicMaterial
    block = '4'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
    C_matrix ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 6.0e6 0.5e6 0.5e6 1.5e6'

    euler_angle_1 = 0.0
    euler_angle_2 = 83.0
    euler_angle_3 = 0.0
  [../]  

  [./Goo5]
    type = LinearGeneralAnisotropicMaterial
    block = '5'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
    C_matrix ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 6.0e6 0.5e6 0.5e6 1.5e6'

    euler_angle_1 = 0.0
    euler_angle_2 = 83.0
    euler_angle_3 = 0.0
  [../]  

  [./Goo6]
    type = LinearGeneralAnisotropicMaterial
    block = '6'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false
    C_matrix ='1.0e6 0.5e6 0.5e6 1.0e6 0.5e6 6.0e6 0.5e6 0.5e6 1.5e6'

    euler_angle_1 = 0.0
    euler_angle_2 = 83.0
    euler_angle_3 = 0.0
  [../]
  

[] # Materials

#[Preconditioning]
#   type = FDP
#[]
[Executioner]

  type = Steady
#  petsc_options = '-ksp_monitor'
#  petsc_options_iname = '-ksp_type -pc_type'
#  petsc_options_value = 'gmres lu'

#  type = Transient
petsc_options = '-snes_mf_operator -snes_view -snes_monitor -ksp_monitor -dm_view -fieldsplit_disp_x_ksp_monitor -fieldsplit_disp_y_disp_z_fieldsplit_y_ksp_monitor'
petsc_options_iname = '-pc_type    -pc_fieldsplit_decomposition -pc_fieldsplit_type -fieldsplit_disp_x_disp_y_pc_type  -fieldsplit_disp_x_disp_y_fieldsplit_type    -fieldsplit_disp_x_disp_y_pc_fieldsplit_decomposition -fieldsplit_disp_x_disp_y_fieldsplit_disp_x_pc_type -fieldsplit_disp_x_disp_y_fieldsplit_disp_x_pc_asm_blocks -fieldsplit_disp_x_disp_y_fieldsplit_disp_x_sub_pc_type -fieldsplit_disp_x_disp_y_fieldsplit_disp_y_pc_type -fieldsplit_disp_x_disp_y_fieldsplit_disp_y_pc_asm_blocks -fieldsplit_disp_x_disp_y_fieldsplit_disp_y_sub_pc_type -fieldsplit_disp_z_pc_type -fieldsplit_disp_z_pc_asm_blocks -fieldsplit_disp_z_sub_pc_type'

petsc_options_value = 'fieldsplit     var:disp_x,disp_y;disp_z;            schur                     fieldsplit                                     schur                                                 var                                                   asm                                                               10                                                  lu                                                      asm                                                      10                                                       lu                          asm                              10                             lu'

  nl_abs_tol = 1e-10
#  l_abs_tol  = 1e-10

  l_max_its = 30

#  start_time = 0.0
#  dt = 1.0
#  num_steps = 2
#  end_time = 2.0
[] # Executioner

[Output]
  file_base = out_6grain_example_fine_uniax_BC_mc
  interval = 1
  output_initial = true
  elemental_as_nodal = true
  exodus = true
  perf_log = true
[] # Output
