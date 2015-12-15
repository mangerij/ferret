[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 12
  ny = 12
  nz = 6
  xmin = -6
  xmax = 6
  ymin = -6
  ymax = 6
  zmin = -2.5
  zmax = 2.5
  elem_type = HEX8
[]

[GlobalParams]
  len_scale = 1.0
  alpha1 = -0.1722883 # (3.766(T-765.1)*10^5) C^{-2} nm^2 (T = 300 K)
  alpha11 = -0.07253
  alpha111 = 0.26
  alpha12 = 0.75
  alpha112 = 0.61
  alpha123 = -3.67

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
[]


[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    #[./InitialCondition]
    #  type = ConstantIC
    #  value = 0.1
    #[../]
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
  [../]
[]


[Kernels]
  #Bulk energy density
  [./bed_x]
    type = BulkEnergyDerivativeSixth
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = BulkEnergyDerivativeSixth
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = BulkEnergyDerivativeSixth
    variable = polar_z
    component = 2
  [../]
  ##Time dependence
  [./polar_x_time]
    type=TimeDerivativeScaled
    variable = polar_x
    time_scale = 1.0
  [../]
  [./polar_y_time]
    type=TimeDerivativeScaled
    variable = polar_y
    time_scale = 1.0
  [../]
  [./polar_z_time]
    type=TimeDerivativeScaled
    variable = polar_z
    time_scale = 1.0
  [../]
[]


[BCs]

[]

[Postprocessors]
  [./bulk_energy]
   type = BulkEnergy
   execute_on = 'timestep_end'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_rtol -ksp_rtol -pc_type    -sub_pc_factor_zeropivot -pc_factor_zeropivot -pc_side '
    petsc_options_value = '    121            1e-8	 1e-12      gamg                   1e-50    1e-50	  left        '
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dt = 0.85
  num_steps = 90
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_bulk6
    output_initial = true
    elemental_as_nodal = true
    interval = 15
  [../]
[]
