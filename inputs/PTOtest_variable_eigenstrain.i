
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 26
  ny = 26
  nz = 11
  xmin = -12
  xmax = 12
  ymin = -12
  ymax = 12
  zmin = -5
  zmax = 5
  elem_type = HEX8
[]

[GlobalParams]
  len_scale = 1.0
  alpha1 = -0.1722883 # (3.766(T-765.1)*10^5) C^{-2} nm^2 (T = 293 K)
  alpha11 = -0.07253
  alpha111 = 0.26
  alpha12 = 0.75
  alpha112 = 0.61
  alpha123 = -3.67
  G110 = 0.173
  G11_G110 = 0.6
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
  #prefactor = 0.01 #negative = tension, positive = compression
[]


[Functions]
  [./parsed_function]
    type = ParsedFunction
    value = '0.01 * sin(0.05 * x) * sin(0.05 * y)'
  [../]
  [./parsed_grad_function]
    type =ParsedGradFunction
    value = '0.02 * sin(x) * sin(y)'
    grad_x = '0.02 *  cos(x) * sin(y)'
    grad_y = '0.02 *  sin(x) * cos(y)'
  [../]
[]

[AuxVariables]
  [./a]
   order = THIRD
   family = HERMITE
  [../]
[]

[ICs]
  [./a_ic]
    type = FunctionIC
    variable = 'a'
    function = parsed_function
  [../]
[]

[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -1.0e-5
      max = 1.0e-5
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -1.0e-5
      max = 1.0e-5
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -1.0e-5
      max = 1.0e-5
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

  [./var_dependence]
    type = DerivativeParsedMaterial
    block = 0
    function = a
    args = a
    output_properties = 'var_dep'
    f_name = var_dep
    enable_jit = true
    derivative_order = 2
  [../]


  [./eigen_strain]
    type = ComputeVariableEigenstrain
    block = 0
   eigen_base = '1 1 1 1 1 1 1 1 1'

    prefactor = var_dep
    args = a
  [../]


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

  [./slab_ferroelectric]
    type = ComputeElectrostrictiveTensor
    Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.03375 0.03375 0.03375'
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
  [../]

[]


[Kernels]
  #Elastic problem
  [./TensorMechanics]
  #This is an action block
  [../]
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
  ##Wall energy penalty
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
##Polarization-strain coupling

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
  ##Electrostatics
  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
     permittivity = 0.08854187
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
  ##Time dependence
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
#  [./disp_x_front]
#    type = DirichletBC
#    boundary = 'front'
#    value = 0.0
#    variable = disp_x
#  [../]
#  [./disp_y_front]
#    type = DirichletBC
#    boundary = 'front'
#    value = 0.0
#    variable = disp_y
#  [../]
#  [./disp_z_front]
#    type = DirichletBC
#    boundary = 'front'
#    value = 0.0
#    variable = disp_z
#  [../]
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
      execute_on = 'timestep_end'
    [../]
    [./Fwall]
      type = WallEnergy
      execute_on = 'timestep_end'
    [../]
    [./Felastic]
      type = ElasticEnergy
      execute_on = 'timestep_end'
    [../]
    [./Fcoupled]
      type = CoupledEnergy
      execute_on = 'timestep_end'
    [../]
    [./Felec]
      type = ElectrostaticEnergy
      permittivity = 0.08854187
      execute_on = 'timestep_end'
    [../]
    [./Ftotal]
      type = TotalEnergyFlow
      Fbulk = Fbulk
      Fwall = Fwall
      Fcoupled = Fcoupled
      Felec = Felec
      execute_on = 'timestep_end'
    [../]
    [./perc_change]
     type = PercentChangePostprocessor
     postprocessor = Ftotal
   [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121                1e-6      1e-8    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
    [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    #iteration_window = 3
    optimal_iterations = 6 #should be 5 probably
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.8
[../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  #dt = 0.5
  dtmin = 1e-13
  dtmax = 0.1
  num_steps = 50
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = outPTOchunk_test_var_eigen
    elemental_as_nodal = true
    interval = 1
  [../]
[]

