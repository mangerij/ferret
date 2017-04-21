
[Mesh]
  file = exodus_disk_r8_h1.e
[]


[GlobalParams]
  len_scale = 1.0

  alpha11 = -0.07253
  alpha111 = 0.26
  alpha12 = 0.75
  alpha112 = 0.61
  alpha123 = -3.67

  G110 = 0.173
  G11/G110 = 1.6
  G12/G110 = 0
  G44/G110 = 0.8
  G44P/G110 = 0.8

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  potential_int = potential_int

  # Elastic interactions are turned off for now.
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z

  displacements = 'disp_x disp_y disp_z'

  # -0.1722883 #room temp PTO
  alpha1 = -0.1722883

  # negative = tension, positive = compression
  prefactor = 0.01

  # value runs from 0 to 0.005 I believe
  lambda = 0.01

[]


[Functions]
  #-----------------------------------------------------#
  # This is a skyrmion solution from Svitlana and Igor  #
  # It was fit in Mathematica                           #
  #-----------------------------------------------------#


  [./parsed_function_x_skyrm]
    type = ParsedFunction
    value = '-(0.738217-0.00686984*(x^2+y^2)^(0.5)+0.00644497*(x^2+y^2)-0.0188174*(x^2+y^2)^(1.5)+0.00441745*(x^2+y^2)^2-0.000274842*(x^2+y^2)^(5/2))*sin(-0.028395+0.267482*(x^2+y^2)^(0.5)-0.146762*(x^2+y^2)+0.0632932*(x^2+y^2)^(1.5)-0.00790942*(x^2+y^2)^(2)+0.000294936*(x^2+y^2)^(5/2))*sin(atan(y/x))'
  [../]
  [./parsed_function_y_skyrm]
    type = ParsedFunction
    value = '(0.738217-0.00686984*(x^2+y^2)^(0.5)+0.00644497*(x^2+y^2)-0.0188174*(x^2+y^2)^(1.5)+0.00441745*(x^2+y^2)^2-0.000274842*(x^2+y^2)^(5/2))*sin(-0.028395+0.267482*(x^2+y^2)^(0.5)-0.146762*(x^2+y^2)+0.0632932*(x^2+y^2)^(1.5)-0.00790942*(x^2+y^2)^(2)+0.000294936*(x^2+y^2)^(5/2))*cos(atan(y/x))'
  [../]
  [./parsed_function_z_skyrm]
    type = ParsedFunction
    value = '(0.738217-0.00686984*(x^2+y^2)^(0.5)+0.00644497*(x^2+y^2)-0.0188174*(x^2+y^2)^(1.5)+0.00441745*(x^2+y^2)^2-0.000274842*(x^2+y^2)^(5/2))*cos(-0.028395+0.267482*(x^2+y^2)^(0.5)-0.146762*(x^2+y^2)+0.0632932*(x^2+y^2)^(1.5)-0.00790942*(x^2+y^2)^(2)+0.000294936*(x^2+y^2)^(5/2))'
  [../]

  #-----------------------------------------------------#
  # This is a cylindrical domain solution from Svitlana #
  # and Igor. It was fit in Mathematica.                #
  #-----------------------------------------------------#


  [./parsed_function_z_cd]
    type = ParsedFunction
    value = '0.719527-0.0061793*(x^2+y^2)^(0.5)-0.00641062*(x^2+y^2)+0.00508983*(x^2+y^2)^(1.5)-0.0020986*(x^2+y^2)^2+0.000171088*(x^2+y^2)^(5/2)'
  [../]
[]




[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = parsed_function_x_skyrm
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = parsed_function_y_skyrm
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = parsed_function_z_skyrm
    [../]
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
    block = '1'
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


[AuxVariables]
  [./chern]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./chernMag]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pMag]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./divP]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./cherndens]
    type = ChernSimonsDensity  #Need to code in skyrmion number instead of this
    variable = chern
  [../]
  [./chernMagdens]
    type = ChernSimonsDensityMag  #Need to code in skyrmion number instead of this
    variable = chernMag
  [../]
  [./pmag]
    type = PolarMag
    variable = pMag
  [../]
  [./divPd]
    type = DivP
    variable = divP
  [../]
[]


[Materials]
  [./eigen_strain_zz] #Use for stress-free strain (ie epitaxial)
    type = ComputeEigenstrain
    block = '1'
    # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
    eigen_base = '1 0 0 0 1 0 0 0 0'
    eigenstrain_name = eigenstrain
  [../]
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #from MaterialsProject
    C_ijkl = '281 115.74 115.74 281 115.74 281 97.18 97.18 97.18'
    block = '1'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    block = '1'
    eigenstrain_names = eigenstrain
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
    block = '1'
  [../]

  [./slab_ferroelectric]
    block = '1'
    type = ComputeElectrostrictiveTensor
    Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.03375 0.03375 0.03375'
    C_ijkl = '281 115.74 115.74 281 115.74 281 97.18 97.18 97.18'
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
      block = '1'
      variable = disp_x
      component = 0
    [../]
    [./ferroelectriccouplingX_yy]
      type = FerroelectricCouplingX
      block = '1'
      variable = disp_y
      component = 1
    [../]
    [./ferroelectriccouplingX_zz]
      type = FerroelectricCouplingX
      block = '1'
      variable = disp_z
      component = 2
    [../]

    [./depol_z]
      type = DepolEnergy
      permitivitty = 0.00885
      variable = polar_z
      avePz = avePz
    [../]


    ##Electrostatics
    [./polar_x_electric_E]
       type = PolarElectricEStrong
       variable = potential_int
    [../]
    [./FE_E_int]
       type = Electrostatics
       variable = potential_int
       block = '1'
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
   [./top_electrode_top]
      type = DirichletBC
      variable = 'potential_int'
      value = 0.0001
      boundary = '2'
   [../]

   [./top_electrode_bottom]
      type = DirichletBC
      variable = 'potential_int'
      value = 0.0001
      boundary = '3'
   [../]

   [./radial_side]
     type = NeumannBC
     variable = 'polar_x polar_y polar_z'
     value = 0
     boundary = '1'
   [../]
[]

[Postprocessors]
   [./pmagave]
     type = ElementAverageValue
     variable = pMag
     execute_on = 'timestep_end'
   [../]
   [./avePz]
     type = ElementAverageValue
     variable = polar_z
     execute_on = 'initial linear nonlinear timestep_begin timestep_end'
   [../]
   [./avgChern]
     block = '1'
     type = ElementAverageValue
     variable = chern
   [../]
   [./avgChernMag]
     block = '1'
     type = ElementAverageValue
     variable = chernMag
   [../]
   [./Fbulk]
      type = BulkEnergy
      block = '1'
      execute_on = 'initial timestep_end'
    [../]
    [./Fwall]
      type = WallEnergy
      block = '1'
      execute_on = 'initial timestep_end'
    [../]
    [./Felastic]
      type = ElasticEnergy
      block = '1'
      execute_on = 'initial timestep_end'
    [../]
    [./Fcoupled]
      block = '1'
      type = CoupledEnergy
      execute_on = 'initial timestep_end'
    [../]
    [./Fdepol]
      block = '1'
      type = DepolarizationEnergy
      execute_on = 'initial timestep_end'
    [../]
    [./Felec]
      block = '1'
      type = ElectrostaticEnergy
      execute_on = 'initial timestep_end'
    [../]
    [./Ftotal]
      type = TotalEnergyFlow
      Fbulk = Fbulk
      Fwall = Fwall
      Fcoupled = Fcoupled
      Felec = Felec
      #Fdepol = Fdepol
      execute_on = 'initial timestep_end'
    [../]
    [./perc_change]
     type = PercentChangePostprocessor
     postprocessor = Ftotal
   [../]
[]


[UserObjects]
 [./kill]
  type = Terminator
  expression = 'perc_change <= 2.0e-3'
 [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    250               1e-10      1e-8     1e-8    bjacobi'
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
  dtmin = 1e-13
  dtmax = 0.7
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    execute_on = 'timestep_end'
    file_base = out_nanodisk_skyrm_0_21st
    elemental_as_nodal = true
    interval = 1
  [../]
  [./outcsv]
    type = CSV
    file_base = out_nanodisk_skyrm_0_21st
    execute_on = 'timestep_end'
  [../]
[]
