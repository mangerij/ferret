###########################
# This example demonstrates the problem due to  interaction of the scaling of C_ijkl and the snes linear search.
# We consider two different C_ijkl, see the [Materials] block. These two C_ijkl are only differ by a scale of 1e3.
# For the first  case, the default snes line search (which is cubit) doesn't converge,  but '-snes_linesearch_type basic' converges to the correct result.
# For the second case, both default and basic line search converges to the correct results.
# For both case, '-snes_type test' indicates the Jacobian is correctly calculated.
# The above numerics might suggests the following possible causes.
# (1) The cubit line search has a bug. (2) A terrible conditioning of the Jacobians causes numerical instability which leads to the failure of the linear search. This terrible scaling might due to the fact that MOOSE handles the boundary conditions in a naive way: they simply set the Dirichlet BC nodes to be one which theorectically might lead to a very worse conditioning of the system even if the discrete differential operator is not ill-condition. This might be what we see in this example.
###########################
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

[TensorMechanics]
  [./solid]
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
[]

[BCs]
  [./anchor_upz_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'upz_outter'
    value = 0
  [../]
  [./anchor_upz_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'upz_outter'
    value = 0
  [../]
  [./anchor_upz_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'upz_outter'
    value = 1.0
  [../]
  [./anchor_downz_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'downz_outter'
    value = 0
  [../]
  [./anchor_downz_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'downz_outter'
    value = 0
  [../]
  [./anchor_downz_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'downz_outter'
    value = -1.0
  [../]
[]

[Materials]
  [./whole]
    type = LinearElasticMaterial
    block = 'exterior interior'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false

    # isotropic bulk modulus=1.5, shear modulus=1.5
    C_ijkl='3.e6 1.e6 1.e6 3.e6 1.e6 3.e6 1.e6 1.e6 1.e6'
    #C_ijkl='3.e3 1.e3 1.e3 3.e3 1.e3 3.e3 1.e3 1.e3 1.e3'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]
[] # Materials

[Preconditioning]
   [./smp]
     type=SMP
     full=true
     pc_side=left
   [../]
[]

[Executioner]
  type = Steady
  solve_type=newton
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  #petsc_options_iname='-pc_type -snes_linesearch_type'
  #petsc_options_value='jacobi     basic'
  petsc_options_iname='-pc_type'
  petsc_options_value='jacobi'
[]

[Outputs]
  output_initial = true
  exodus = true
  [./console]
    type = Console
    perf_log = true
  [../]
[]