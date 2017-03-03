[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 8
  ny = 8
  nz = 8
  xmin = -8
  xmax = 8
  ymin = -8
  ymax = 8
  zmin = -8
  zmax = 8
  elem_type = HEX8
[]

[GlobalParams]
  alpha1 = -0.1524
  alpha2 = 1.76  
  alpha3 = 3.73
  alpha4 = -.591
  alpha5 = -.651
  x1 = -4.18
  x2 = -36.9
  x3 = -251
  x4 = 121
  x5 = 2410
  x6 = 6490
  len_scale = 1.0
  G110 = 0.173
  G11/G110 = 2.0
  G12/G110 = 0
  G44/G110 = 1.0
  G44P/G110 = 1.0
  T = 0.0
  Tc = 120.0
  epsilon = 0.02
  permittivity = 0.5843763
  polar_x = polar_x
  polar_y = polar_y

[]



[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = ConstantIC
      value = 0.05
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = ConstantIC
      value = 0.05
    [../]
  [../]
[]

[Kernels]
  #Bulk energy density
  [./bed_x]
    type = BulkEnergyDerivativePSTO
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = BulkEnergyDerivativePSTO
    variable = polar_y
    component = 1
  [../]
 

[]

[BCs]

   [./polar_x]
    type = DirichletBC
    boundary = 'back'
    variable = polar_x
    value = 0
  [../]
   [./polar_y]
    type = DirichletBC
    boundary = 'back'
    variable = polar_y
    value = 0
  [../]

[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason '
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type    -pc_factor_zeropivot'
    petsc_options_value = '    121            1e-8      1e-8    gamg           1e-50     '
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_wall_test
    output_initial = true
    elemental_as_nodal = true
  [../]
[]
