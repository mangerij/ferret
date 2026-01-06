
##see J. Appl. Phys. 126, 151101 (2019);

Nx = 2
Ny = 2
Nz = 2

xMax = 0.01  # in micrometers = 10 nm
yMax = 0.01
zMax = 0.01


[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${Nx}
    ny = ${Ny}
    nz = ${Nz}
    xmin = 0.0
    xmax = ${xMax}
    ymin = 0.0
    ymax = ${yMax}
    zmin = 0.0
    zmax = ${zMax}
    elem_type = HEX8
  []
[]

[GlobalParams]
  mag1_x = mag1_x
  mag1_y = mag1_y
  mag1_z = mag1_z

  mag2_x = mag2_x
  mag2_y = mag2_y
  mag2_z = mag2_z
[]

[Materials]
  [./constants_kOe] # Constants used in other material properties
    type = GenericConstantMaterial
    prop_names = ' H0         Ms        g0          He        Ha      '
    prop_values = '8.0e-5    1.0      2.8e9      0.000526     8.2e-6      '
  [../]
  [./constants] # Constants used in other material properties
    type = GenericConstantMaterial
    prop_names = ' alpha     mu0   nx ny nz   long_susc t'
    prop_values = '0.0       1256   0  0  1          1.0     0'
  [../]

  [./a_long]
    type = GenericFunctionMaterial
    prop_names = 'alpha_long'
    prop_values = 'bc_func_1'
  [../]
[]

[Functions]
  [./bc_func_1]
    type = ParsedFunction
    expression = 'st'
    symbol_names = 'st'
    symbol_values = '0.5'
  [../]
[]

[Variables]
  [./mag1_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi1
      theta = polar_theta1
      M0s = 1.0 #amplitude of the RandomConstrainedVectorFieldIC
      component  = 0
    [../]
  [../]
  [./mag1_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi1
      theta = polar_theta1
      M0s = 1.0
      component  = 1
    [../]
  [../]
  [./mag1_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi1
      theta = polar_theta1
      M0s = 1.0
      component  = 2
    [../]
  [../]

  [./mag2_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi2
      theta = polar_theta2
      M0s = 1.0
      component  = 0
    [../]
  [../]
  [./mag2_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi2
      theta = polar_theta2
      M0s = 1.0
      component  = 1
    [../]
  [../]
  [./mag2_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi2
      theta = polar_theta2
      M0s = 1.0
      component  = 2
    [../]
  [../]

[]

[AuxVariables]

  #--------------------------------------------#
  #                                            #
  #  field to seed IC that obeys constraint    #
  #                                            #
  #--------------------------------------------#

  [./azimuth_phi1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.3
      max = 0.35
      seed = 2
    [../]
  [../]
  [./polar_theta1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.0001
      max = 0.0002
      seed = 37
    [../]
  [../]

  [./azimuth_phi2]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.3
      max = 0.35
      seed = 2
    [../]
  [../]
  [./polar_theta2]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 3.1415
      max = 3.1416
      seed = 37
    [../]
  [../]

  [./mag1_s]
    order = FIRST
    family = LAGRANGE
  [../]
  [./mag2_s]
    order = FIRST
    family = LAGRANGE
  [../]

  [./Neel_L_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Neel_L_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Neel_L_z]
    order = FIRST
    family = LAGRANGE
  [../]

  [./SSMag_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./SSMag_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./SSMag_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  #---------------------------------------#
  #                                       #
  #          Time dependence              #
  #                                       #
  #---------------------------------------#

  [./mag1_x_time]
    type = TimeDerivative
    variable = mag1_x
  [../]
  [./mag1_y_time]
    type = TimeDerivative
    variable = mag1_y
  [../]
  [./mag1_z_time]
    type = TimeDerivative
    variable = mag1_z
  [../]

  [./mag2_x_time]
    type = TimeDerivative
    variable = mag2_x
  [../]
  [./mag2_y_time]
    type = TimeDerivative
    variable = mag2_y
  [../]
  [./mag2_z_time]
    type = TimeDerivative
    variable = mag2_z
  [../]

  #---------------------------------------#
  #                                       #
  #     AFM resonance kernel terms        #
  #                                       #
  #---------------------------------------#

  [./afmr1_x]
    type = UniaxialAFMSublattice
    variable = mag1_x
    mag_sub = 0
    component = 0
  [../]
  [./afmr1_y]
    type = UniaxialAFMSublattice
    variable = mag1_y
    mag_sub = 0
    component = 1
  [../]
  [./afmr1_z]
    type = UniaxialAFMSublattice
    variable = mag1_z
    mag_sub = 0
    component = 2
  [../]

  [./afmr2_x]
    type = UniaxialAFMSublattice
    variable = mag2_x
    mag_sub = 1
    component = 0
  [../]
  [./afmr2_y]
    type = UniaxialAFMSublattice
    variable = mag2_y
    mag_sub = 1
    component = 1
  [../]
  [./afmr2_z]
    type = UniaxialAFMSublattice
    variable = mag2_z
    mag_sub = 1
    component = 2
  [../]

  #---------------------------------------#
  #                                       #
  #          LLB constraint terms         #
  #                                       #
  #---------------------------------------#

  [./llb1_x]
    type = LongitudinalLLB
    variable = mag1_x
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    component = 0
  [../]
  [./llb1_y]
    type = LongitudinalLLB
    variable = mag1_y
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    component = 1
  [../]

  [./llb1_z]
    type = LongitudinalLLB
    variable = mag1_z
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    component = 2
  [../]

  [./llb2_x]
    type = LongitudinalLLB
    variable = mag2_x
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    component = 0
  [../]
  [./llb2_y]
    type = LongitudinalLLB
    variable = mag2_y
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    component = 1
  [../]

  [./llb2_z]
    type = LongitudinalLLB
    variable = mag2_z
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    component = 2
  [../]
[]


[AuxKernels]
  [./mag1_mag]
    type = VectorMag
    variable = mag1_s
    vector_x = mag1_x
    vector_y = mag1_y
    vector_z = mag1_z
    execute_on = 'initial timestep_end final'
  [../]

  [./mag2_mag]
    type = VectorMag
    variable = mag2_s
    vector_x = mag2_x
    vector_y = mag2_y
    vector_z = mag2_z
    execute_on = 'initial timestep_end final'
  [../]


  [./Neel_Lx]
    type = VectorDiffOrSum
    variable = Neel_L_x
    var1 = mag1_x
    var2 = mag2_x
    diffOrSum = 0
    execute_on = 'initial timestep_end final'
  [../]
  [./Neel_Ly]
    type = VectorDiffOrSum
    variable = Neel_L_y
    var1 = mag1_y
    var2 = mag2_y
    diffOrSum = 0
    execute_on = 'initial timestep_end final'
  [../]
  [./Neel_Lz]
    type = VectorDiffOrSum
    variable = Neel_L_z
    var1 = mag1_z
    var2 = mag2_z
    diffOrSum = 0
    execute_on = 'initial timestep_end final'
  [../]

  [./smallSignalMag_x]
    type = VectorDiffOrSum
    variable = SSMag_x
    var1 = mag1_x
    var2 = mag2_x
    diffOrSum = 1
    execute_on = 'initial timestep_end final'
  [../]
  [./smallSignalMag_y]
    type = VectorDiffOrSum
    variable = SSMag_y
    var1 = mag1_y
    var2 = mag2_y
    diffOrSum = 1
    execute_on = 'initial timestep_end final'
  [../]
  [./smallSignalMag_z]
    type = VectorDiffOrSum
    variable = SSMag_z
    var1 = mag1_z
    var2 = mag2_z
    diffOrSum = 1
    execute_on = 'initial timestep_end final'
  [../]

[]

[Kernels]
  #---------------------------------------#
  #                                       #
  #          Time dependence              #
  #                                       #
  #---------------------------------------#

  [./mag1_x_time]
    type = TimeDerivative
    variable = mag1_x
  [../]
  [./mag1_y_time]
    type = TimeDerivative
    variable = mag1_y
  [../]
  [./mag1_z_time]
    type = TimeDerivative
    variable = mag1_z
  [../]

  [./mag2_x_time]
    type = TimeDerivative
    variable = mag2_x
  [../]
  [./mag2_y_time]
    type = TimeDerivative
    variable = mag2_y
  [../]
  [./mag2_z_time]
    type = TimeDerivative
    variable = mag2_z
  [../]

  #---------------------------------------#
  #                                       #
  #     AFM resonance kernel terms        #
  #                                       #
  #---------------------------------------#

  [./afmr1_x]
    type = UniaxialAFMSublattice
    variable = mag1_x
    mag_sub = 0
    component = 0
  [../]
  [./afmr1_y]
    type = UniaxialAFMSublattice
    variable = mag1_y
    mag_sub = 0
    component = 1
  [../]
  [./afmr1_z]
    type = UniaxialAFMSublattice
    variable = mag1_z
    mag_sub = 0
    component = 2
  [../]

  [./afmr2_x]
    type = UniaxialAFMSublattice
    variable = mag2_x
    mag_sub = 1
    component = 0
  [../]
  [./afmr2_y]
    type = UniaxialAFMSublattice
    variable = mag2_y
    mag_sub = 1
    component = 1
  [../]
  [./afmr2_z]
    type = UniaxialAFMSublattice
    variable = mag2_z
    mag_sub = 1
    component = 2
  [../]

  #---------------------------------------#
  #                                       #
  #          LLB constraint terms         #
  #                                       #
  #---------------------------------------#

  [./llb1_x]
    type = LongitudinalLLB
    variable = mag1_x
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    component = 0
  [../]
  [./llb1_y]
    type = LongitudinalLLB
    variable = mag1_y
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    component = 1
  [../]

  [./llb1_z]
    type = LongitudinalLLB
    variable = mag1_z
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    component = 2
  [../]

  [./llb2_x]
    type = LongitudinalLLB
    variable = mag2_x
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    component = 0
  [../]
  [./llb2_y]
    type = LongitudinalLLB
    variable = mag2_y
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    component = 1
  [../]

  [./llb2_z]
    type = LongitudinalLLB
    variable = mag2_z
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    component = 2
  [../]
[]

[BCs]
  #---------------------------------------#
  #                                       #
  #  periodic magnetization distribution  #
  #                                       #
  #---------------------------------------#

  [./Periodic]
    [./xyz]
      auto_direction = 'x y z'
      variable = 'mag1_x mag1_y mag1_z mag2_x mag2_y mag2_z'
    [../]
  [../]
[]

[Postprocessors]
   [./dt]
     type = TimestepSize
   [../]

  #---------------------------------------#
  #                                       #
  #       Average |M| and along other     #
  #       directions                      #
  #                                       #
  #---------------------------------------#

  [./<M1>]
    type = ElementAverageValue
    variable = mag1_s
    execute_on = 'initial timestep_end final'
  [../]
  [./<M2>]
    type = ElementAverageValue
    variable = mag2_s
    execute_on = 'initial timestep_end final'
  [../]

  [./<m1x>]
    type = ElementAverageValue
    variable = mag1_x
    execute_on = 'initial timestep_end final'
  [../]
  [./<m1y>]
    type = ElementAverageValue
    variable = mag1_y
    execute_on = 'initial timestep_end final'
  [../]
  [./<m1z>]
    type = ElementAverageValue
    variable = mag1_z
    execute_on = 'initial timestep_end final'
  [../]

  [./<m2x>]
    type = ElementAverageValue
    variable = mag2_x
    execute_on = 'initial timestep_end final'
  [../]
  [./<m2y>]
    type = ElementAverageValue
    variable = mag2_y
    execute_on = 'initial timestep_end final'
  [../]
  [./<m2z>]
    type = ElementAverageValue
    variable = mag2_z
    execute_on = 'initial timestep_end final'
  [../]

  [./<Lx>]
    type = ElementAverageValue
    variable = Neel_L_x
    execute_on = 'initial timestep_end final'
  [../]
  [./<Ly>]
    type = ElementAverageValue
    variable = Neel_L_y
    execute_on = 'initial timestep_end final'
  [../]
  [./<Lz>]
    type = ElementAverageValue
    variable = Neel_L_z
    execute_on = 'initial timestep_end final'
  [../]

  [./<SSmx>]
    type = ElementAverageValue
    variable = SSMag_x
    execute_on = 'initial timestep_end final'
  [../]
  [./<SSmy>]
    type = ElementAverageValue
    variable = SSMag_y
    execute_on = 'initial timestep_end final'
  [../]
  [./<SSmz>]
    type = ElementAverageValue
    variable = SSMag_z
    execute_on = 'initial timestep_end final'
  [../]
[]


[Preconditioning]
  #---------------------------------------#
  #                                       #
  #            Solver options             #
  #                                       #
  #---------------------------------------#

  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121               1e-8      1e-6      1e-6      bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  [./TimeIntegrator]
    type = NewmarkBeta
  [../]
  dtmin = 1e-11
  dtmax = 1.0e-3
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 25
    iteration_window = 2000
    growth_factor = 1.3
    cutback_factor = 0.8
    dt = 12.0e-9
  [../]
  verbose = true

  num_steps = 18000
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_AFMR_ex
    elemental_as_nodal = true
    time_step_interval = 500
  [../]
  [./outCSV]
    type = CSV
    file_base = out_AFMR_ex
    time_step_interval = 1
  [../]
[]
