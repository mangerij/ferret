[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  nx = 10
  ny = 10
  elem_type = QUAD4
[]


[GlobalParams]
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  antiferrodis_A_x = antiferrodis_A_x
  antiferrodis_A_y = antiferrodis_A_y
  antiferrodis_A_z = antiferrodis_A_z

  Ap = -1.09028
  Bp = 6.68020e-1
  Cp = -4.67264e-1

  Ar = -5.33957e-3
  Br = 1.35373e-5
  Cr = -7.73973e-6

  Bpr = 6.99618e-3
  Cpr = -2.14508e-3
  Cprp = -1.40134e-2

[]


[Variables]
  [./diffused]
    order = FIRST
    family = LAGRANGE
  [../]


  # ODE variables
  [./polar_x]
    family = SCALAR
    order = FIRST
    initial_condition = 0.01
  [../]
  [./polar_y]
    family = SCALAR
    order = FIRST
    initial_condition = 0.01
  [../]
  [./polar_z]
    family = SCALAR
    order = FIRST
    initial_condition = 0.01
  [../]
  [./antiferrodis_A_x]
    family = SCALAR
    order = FIRST
    initial_condition = 0.1
  [../]
  [./antiferrodis_A_y]
    family = SCALAR
    order = FIRST
    initial_condition = 0.1
  [../]
  [./antiferrodis_A_z]
    family = SCALAR
    order = FIRST
    initial_condition = 0.1
  [../]



[]

[Kernels]
  [./td]
    type = TimeDerivative
    variable = diffused
  [../]
  [./diff]
    type = Diffusion
    variable = diffused
  [../]
[]

[ScalarKernels]
  [./td_px]
    type = ODETimeDerivative
    variable = polar_x
  [../]
  [./td_py]
    type = ODETimeDerivative
    variable = polar_y
  [../]
  [./td_pz]
    type = ODETimeDerivative
    variable = polar_z
  [../]

  [./td_Ax]
    type = ODETimeDerivative
    variable = antiferrodis_A_x
  [../]
  [./td_Ay]
    type = ODETimeDerivative
    variable = antiferrodis_A_y
  [../]
  [./td_Az]
    type = ODETimeDerivative
    variable = antiferrodis_A_z
  [../]


  [./bulk_px]
    type = ScalarBulkEnergyP
    variable = polar_x
    component = 0
  [../]
  [./bulk_py]
    type = ScalarBulkEnergyP
    variable = polar_y
    component = 1
  [../]
  [./bulk_pz]
    type = ScalarBulkEnergyP
    variable = polar_z
    component = 2
  [../]

  [./bulk_Ax]
    type = ScalarBulkEnergyA
    variable = antiferrodis_A_x
    component = 0
  [../]
  [./bulk_Ay]
    type = ScalarBulkEnergyA
    variable = antiferrodis_A_y
    component = 1
  [../]
  [./bulk_Az]
    type = ScalarBulkEnergyA
    variable = antiferrodis_A_z
    component = 2
  [../]


  [./roto_Px]
    type = ScalarRotopolarEnergy
    variable = polar_x
    component = 0
  [../]
  [./roto_Py]
    type = ScalarRotopolarEnergy
    variable = polar_y
    component = 1
  [../]
  [./roto_Pz]
    type = ScalarRotopolarEnergy
    variable = polar_z
    component = 2
  [../]
  [./roto_Ax]
    type = ScalarRotopolarEnergy
    variable = antiferrodis_A_x
    component = 3
  [../]
  [./roto_Ay]
    type = ScalarRotopolarEnergy
    variable = antiferrodis_A_y
    component = 4
  [../]
  [./roto_Az]
    type = ScalarRotopolarEnergy
    variable = antiferrodis_A_z
    component = 5
  [../]

[]


[BCs]
  [./right]
    type = DirichletBC
    variable = diffused
    boundary = 1
    value = 0.0
  [../]

  [./left]
    type = DirichletBC
    variable = diffused
    boundary = 0
    value = 1.0
  [../]
[]

[Postprocessors]
  [./Px]
    type = ScalarVariable
    variable = polar_x
    execute_on = timestep_end
  [../]
  [./Py]
    type = ScalarVariable
    variable = polar_y
    execute_on = timestep_end
  [../]
  [./Pz]
    type = ScalarVariable
    variable = polar_z
    execute_on = timestep_end
  [../]

  [./Ax]
    type = ScalarVariable
    variable = antiferrodis_A_x
    execute_on = timestep_end
  [../]
  [./Ay]
    type = ScalarVariable
    variable = antiferrodis_A_y
    execute_on = timestep_end
  [../]
  [./Az]
    type = ScalarVariable
    variable = antiferrodis_A_z
    execute_on = timestep_end
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121            1e-8          1e-8       1e-6     bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  start_time = 0
  dt = 1.0

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'
  num_steps = 50
[]

[Outputs]
  print_linear_residuals = true
  [./out]
    type = Exodus
    file_base = out_LK
  [../]
[]
