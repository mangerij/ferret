
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 32
  ny = 32
  xmin = -8
  xmax = 8
  ymin = -8
  ymax = 8
[]

[GlobalParams]
  len_scale = 1.0
  alpha1 = -0.1722883 # (3.766(T-765.1)*10^5) C^{-2} nm^2 (T = 293 K)
  alpha11 = -0.07253
  alpha111 = 0.26
  alpha12 = 0.75
  alpha112 = 0.61
  alpha123 = -3.67

  polar_x = polar_x
  polar_y = polar_y
  potential_E_int = potential_E_int
  disp_x = disp_x
  disp_y = disp_y
  displacements = 'disp_x disp_y'
  prefactor = 0.0
[]



[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-4
      max = 0.01e-4
      seed = 5
      legacy_generator = true
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-4
      max = 0.01e-4
      seed = 5
      legacy_generator = true
    [../]
  [../]
  [./potential_E_int]
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
[]


[Materials]

  [./Landau_G]
    type = GenericConstantMaterial
    prop_names = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '0.173 0.6 0.0 0.3 0.3'
  [../]

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
    Q_mnkl = '-0.089 0.026 0.026 -0.089 0.026 -0.089 -0.03375 -0.03375 -0.03375'
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
  [../]

  [./permitivitty_1]
    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '0.08854187'
  [../]

[]


[Kernels]
  #Elastic problem
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
  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_E_int
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_E_int
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
[]


[BCs]

  [./potential_E_int_front]
    type = DirichletBC
    boundary = 'top'
    value = 0.0001
    variable = potential_E_int
  [../]

  [./potential_E_int_back]
    type = DirichletBC
    boundary = 'bottom'
    value = 0.0001
    variable = potential_E_int
  [../]

  [./disp_x_back]
    type = DirichletBC
    boundary = 'top'
    value = 0.0
    variable = disp_x
  [../]
  [./disp_y_back]
    type = DirichletBC
    boundary = 'top'
    value = 0.0
    variable = disp_y
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
      type = ElectrostrictiveEnergy
      execute_on = 'initial timestep_end'
    [../]
    [./Felec]
      type = ElectrostaticEnergy
      execute_on = 'initial timestep_end'
    [../]
    [./Ftotal]
      type = LinearCombinationPostprocessor
      pp_names = 'Fbulk Fwall Fcoupled Felec'
      pp_coefs = ' 1 1 1 1'
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
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121               1e-10      1e-8      1e-6    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  scheme = 'implicit-euler'
  #dt = 0.5
  dtmin = 1e-13
  dtmax = 0.1
  num_steps = 10
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = outPTOchunk_test_2D
    elemental_as_nodal = true
    interval = 1
  [../]
[]
