
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 10
  ny = 10
  nz = 7
  xmin = -4
  xmax = 4
  ymin = -4
  ymax = 4
  zmin = -3
  zmax = 3
  elem_type = HEX8
[]

[GlobalParams]
  len_scale = 1.0
  alpha1 = -0.1722883
  alpha11 = -0.07253
  alpha111 = 0.26
  alpha12 = 0.75
  alpha112 = 0.61
  alpha123 = -3.67
  G110 = 0.173
  G11/G110 = 0.6
  G12/G110 = 0
  G44/G110 = 0.3
  G44P/G110 = 0.3
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  displacements = 'disp_x disp_y disp_z'
  prefactor = 0.01 #negative = tension, positive = compression
[]

[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FluctuationsIC
      epsilon = 1.0e-5
      q1 = '354 1005 645'
      q2 = '715 1065 1132'
      q3 = '391 305 1106'
      q4 = '1053 1116 627'
      h = 0.22
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FluctuationsIC
      epsilon = 1.0e-5
      q1 = '354 850 10'
      q2 = '715 28 5'
      q3 = '391 305 1106'
      q4 = '653 1116 627'
      h = 0.15
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FluctuationsIC
      epsilon = 1.0e-5
      q1 = '860 165 645'
      q2 = '715 665 1332'
      q3 = '361 15 706'
      q4 = '253 1116 627'
      h = 0.05
    [../]
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
  [../]
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


[Materials]
  [./eigen_strain_zz] #Use for stress-free strain (ie epitaxial)
   type = ComputeEigenstrain
   block = '0'
  # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
   eigen_base = '1 0 0 0 1 0 0 0 0'
   eigenstrain_name = eigenstrain
 [../]
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    eigenstrain_names = eigenstrain
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
  [../]

  [./slab_ferroelectric]
    type = ComputeElectrostrictiveTensor
    Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.03375 0.03375 0.03375'
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
  [../]

[]


[Kernels]
  [./TensorMechanics]
  [../]
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
  [./walled_x]
    type = WallEnergyDerivative
    variable = polar_x
    component = 0
  [../]
  [./walled_y]
    type = WallEnergyDerivative
    variable = polar_y
    component = 1
  [../]
  [./walled_z]
     type = WallEnergyDerivative
     variable = polar_z
     component = 2
  [../]
  [./ferroelectriccouplingp_xx]
    type = FerroelectricCouplingP
    variable = polar_x
    component = 0
  [../]
  [./ferroelectriccouplingp_yy]
    type = FerroelectricCouplingP
    variable = polar_y
    component = 1
  [../]
  [./ferroelectriccouplingp_zz]
    type = FerroelectricCouplingP
    variable = polar_z
    component = 2
  [../]
  [./ferroelectriccouplingX_xx]
    type = FerroelectricCouplingX
    variable = disp_x
    component = 0
  [../]
  [./ferroelectriccouplingX_yy]
    type = FerroelectricCouplingX
    variable = disp_y
    component = 1
  [../]
  [./ferroelectriccouplingX_zz]
    type = FerroelectricCouplingX
    variable = disp_z
    component = 2
  [../]
  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     permittivity = 0.08854187
  [../]
  [./polar_electric_px]
     type = PolarElectricPStrong
     variable = polar_x
     component = 0
  [../]
  [./polar_electric_py]
     type = PolarElectricPStrong
     variable = polar_y
     component = 1
  [../]
  [./polar_electric_pz]
     type = PolarElectricPStrong
     variable = polar_z
     component = 2
  [../]
  [./polar_x_time]
     type = TimeDerivativeScaled
     variable=polar_x
    time_scale = 1.0
  [../]
  [./polar_y_time]
     type = TimeDerivativeScaled
     variable=polar_y
    time_scale = 1.0
  [../]
  [./polar_z_time]
     type = TimeDerivativeScaled
     variable = polar_z
    time_scale = 1.0
  [../]
[]


[BCs]
  [./disp_x_front]
    type = DirichletBC
    boundary = 'front'
    value = 0.0
    variable = disp_x
  [../]
  [./disp_y_front]
    type = DirichletBC
    boundary = 'front'
    value = 0.0
    variable = disp_y
  [../]
  [./disp_z_front]
    type = DirichletBC
    boundary = 'front'
    value = 0.0
    variable = disp_z
  [../]
  [./potential_int_front]
    type = DirichletBC
    boundary = 'front'
    value = 0.0001
    variable = potential_int
  [../]

  [./potential_int_back]
    type = DirichletBC
    boundary = 'back'
    value = 0.0001
    variable = potential_int
  [../]
  [./disp_x_back]
    type = DirichletBC
    boundary = 'back'
    value = 0.0
    variable = disp_x
  [../]
  [./disp_y_back]
    type = DirichletBC
    boundary = 'back'
    value = 0.0
    variable = disp_y
  [../]
  [./disp_z_back]
    type = DirichletBC
    boundary = 'back'
    value = 0.0
    variable = disp_z
  [../]
[]



[Postprocessors]
   [./Fbulk]
      type = BulkEnergy
      execute_on = 'initial timestep_end'
    [../]
    [./Fwall]
      type = WallEnergy
      execute_on = 'initial timestep_end'
    [../]
    [./Felastic]
      type = ElasticEnergy
      execute_on = 'initial timestep_end'
    [../]
    [./Fcoupled]
      type = CoupledEnergy
      execute_on = 'initial timestep_end'
    [../]
    [./Felec]
      type = ElectrostaticEnergy
      execute_on = 'initial timestep_end'
    [../]
    [./Ftotal]
      type = TotalEnergyFlow
      Fbulk = Fbulk
      Fwall = Fwall
      Fcoupled = Fcoupled
      Felec = Felec
      execute_on = 'initial timestep_end'
    [../]
    [./perc_change]
     type = PercentChangePostprocessor
     postprocessor = Ftotal
     execute_on = 'initial timestep_end'
    [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121                1e-6      1e-8    bjacobi'
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
  dtmin = 1e-13
  dtmax = 0.1
  num_steps = 5
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = outPTOchunk_test_fluct
    elemental_as_nodal = true
    interval = 1
  [../]
[]
