
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 15
  ny = 15
  nz = 10
  xmin = -5
  xmax = 5
  ymin = -5
  ymax = 5
  zmin = -3
  zmax = 3
  elem_type = HEX8
[]



[GlobalParams]
  len_scale = 1.0

  # -0.1722883 #room temp PTO
  alpha1 = -0.5
  alpha11 = 0.25
  alpha111 = 0.0
  alpha12 = 0.0
  alpha112 = 0.0
  alpha123 = 0.0

  G110 = 0.173
  G11_G110 = 1.6
  G12_G110 = 0
  G44_G110 = 0.8
  G44P_G110 = 0.8

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z

  displacements = 'disp_x disp_y disp_z'

  potential_int = potential_int
[]


[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.0001
      max = 0.0001
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.0001
      max = 0.0001
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.0001
      max = 0.0001
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
  [./eigen_strain_zz]
    type = ComputeEigenstrain
    # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
    prefactor = 0.0
    eigen_base = '0 0 0 0 0 0 0 0 0'
    eigenstrain_name = eigenstrain
  [../]
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '200.0 100.0 100.0 200.0 100.0 200.0 85.0 85.0 85.0'
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
    Q_mnkl = '0.09 -0.02 -0.02 0.09 -0.02 0.09 0.03 0.03 0.03'
    C_ijkl = '200.0 100.0 100.0 200.0 100.0 200.0 85.0 85.0 85.0'
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

[]

[Postprocessors]
  [./avePz]
    type = ElementAverageValue
    variable = polar_z
    execute_on = 'timestep_end'
  [../]
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
  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = Fbulk
  [../]
[]



[UserObjects]
 [./kill]
  type = Terminator
  expression = 'perc_change <= 1.0e-3'
 [../]
[]



[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    120               1e-10      1e-8     1e-4    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.12
    optimal_iterations = 6
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.8
  [../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 0.4
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    execute_on = 'timestep_end'
    file_base = out_prob3a
    elemental_as_nodal = true
    interval = 1
  [../]
  [./outcsv]
    type = CSV
    file_base = out_prob3a
    execute_on = 'timestep_end'
  [../]
[]
