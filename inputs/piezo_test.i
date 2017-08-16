

# This file just evolves from the PE state an open circuit BC
# chunk of LNO using coefficients from Chen's appendix 
# (and assuming gradient terms are == PTO).

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 12
  ny = 12
  nz = 8
  xmin = -5
  xmax = 5
  ymin = -5
  ymax = 5
  zmin = -3
  zmax = 3
  elem_type = HEX8
[]

[GlobalParams]
  potential_int = potential_int
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
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
  [./potential_int]
    order = THIRD
    family = HERMITE
  [../]
[]


[Kernels]
  #Elastic problem
  [./TensorMechanics]
  #This is an action block
  [../]
  [./piezocouple_0]
    type = ConversePiezoelectricStrain
    variable = disp_x
    component = 0
  [../]
  [./piezocouple_1]
    type = ConversePiezoelectricStrain
    variable = disp_y
    component = 1
  [../]
  [./piezocouple_2]
    type = ConversePiezoelectricStrain
    variable = disp_z
    component = 2
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     permittivity = 0.08854187
  [../]
  [./disp_x_time]
     type = TimeDerivativeScaled
     variable = disp_x
     time_scale = 1.0
  [../]
  [./disp_y_time]
     type = TimeDerivativeScaled
     variable = disp_y
     time_scale = 1.0
  [../]
  [./disp_z_time]
     type = TimeDerivativeScaled
     variable = disp_z
     time_scale = 1.0
  [../]
[]

[Materials]
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
  [../]
  [./d333]
    type = ComputePiezoTensor
    fill_method2 = symmetric9
    fill_method = general
    compute_piezostrictive_coeff = true
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
    d_ijk = '0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1'
  [../]
[]

[Functions]
  [./bc_func_1]
    type = ParsedFunction
    value = '0.005*sin(omega*t)'
    vars = 'omega'
    vals = '0.1'
  [../]
  [./bc_func_2]
    type = ParsedFunction
    value = '-0.005*sin(omega*t)'
    vars = 'omega'
    vals = '0.1'
  [../]
[]

[BCs]
  # Boundary Condition System
  [./front_pot]
    type = FunctionDirichletBC
    variable = potential_int
    boundary = 'front'
    function = bc_func_1
  [../]
  [./back_pot]
    type = FunctionDirichletBC
    variable = potential_int
    boundary = 'back'
    function = bc_func_2
  [../]
[]




[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart  -snes_atol -ksp_rtol -pc_type'
    petsc_options_value = '    121                1e-10      1e-8    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
    [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    optimal_iterations = 6
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.8
[../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  #dt = 0.5
  dtmin = 1e-13
  dtmax = 0.1
  num_steps = 1000
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_piezo_test
    elemental_as_nodal = true
    interval = 1
  [../]
[]
