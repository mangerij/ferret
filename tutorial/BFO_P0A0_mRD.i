
Dedef = 3.7551
D0def = 0.003
K1def = -5.0068
K1cdef = -0.0550748
Ktdef = -0.00365997

alphadef = 0.1

[Mesh]
  [fileload]
    type = FileMeshGenerator
    file = BFO_P0A0.e
    use_for_exodus_restart = true
  []
[]


[GlobalParams]
  mag1_x = mag1_x
  mag1_y = mag1_y
  mag1_z = mag1_z

  mag2_x = mag2_x
  mag2_y = mag2_y
  mag2_z = mag2_z

  antiphase_A_x = antiphase_A_x
  antiphase_A_y = antiphase_A_y
  antiphase_A_z = antiphase_A_z

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z


[]

[Materials]
  ############################################################################
  ##
  ##       material constants used.
  ##
  ##
  ############################################################################

  [./constants]
    type = GenericConstantMaterial
    prop_names = ' alpha           De       D0             g0mu0Ms      g0          K1         K1c      Kt     permittivity '
    prop_values = '${alphadef}   ${Dedef} ${D0def}     48291.9      48291.9      ${K1def}  ${K1cdef} ${Ktdef}     1.0     '
  [../]

  [./a_long]
    type = GenericFunctionMaterial
    prop_names = 'alpha_long'
    prop_values = 'bc_func_1'
  [../]

[]

[Functions]

  ##############################
  ##
  ## Define the ramping function
  ## expression to be used
  ##
  ##############################

  [./bc_func_1]
    type = ParsedFunction
    expression = 'st'
    vars = 'st'
    vals = '1e3'
  [../]
[]

[Variables]
  [./mag1_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.703126
      max = -0.703125
    [../]
  [../]
  [./mag1_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.0079299
      max = -0.0079298
    [../]
  [../]
  [./mag1_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.711033
      max = 0.711034
    [../]
  [../]

  [./mag2_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.711033
      max = 0.711034
    [../]
  [../]
  [./mag2_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.0079299
      max = -0.0079298
    [../]
  [../]
  [./mag2_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.703126
      max = -0.703125
    [../]
  [../]
[]

[Kernels]

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

  [./afmex1_x]
    type = AFMSublatticeSuperexchange
    variable = mag1_x
    mag_sub = 0
    component = 0
  [../]
  [./afmex1_y]
    type = AFMSublatticeSuperexchange
    variable = mag1_y
    mag_sub = 0
    component = 1
  [../]
  [./afmex1_z]
    type = AFMSublatticeSuperexchange
    variable = mag1_z
    mag_sub = 0
    component = 2
  [../]

  [./afmex2_x]
    type = AFMSublatticeSuperexchange
    variable = mag2_x
    mag_sub = 1
    component = 0
  [../]
  [./afmex2_y]
    type = AFMSublatticeSuperexchange
    variable = mag2_y
    mag_sub = 1
    component = 1
  [../]
  [./afmex2_z]
    type = AFMSublatticeSuperexchange
    variable = mag2_z
    mag_sub = 1
    component = 2
  [../]

  [./afmdmi1_x]
    type = AFMSublatticeDMInteraction
    variable = mag1_x
    mag_sub = 0
    component = 0
  [../]
  [./afmdmi1_y]
    type = AFMSublatticeDMInteraction
    variable = mag1_y
    mag_sub = 0
    component = 1
  [../]
  [./afmdmi1_z]
    type = AFMSublatticeDMInteraction
    variable = mag1_z
    mag_sub = 0
    component = 2
  [../]

  [./afmdmi2_x]
    type = AFMSublatticeDMInteraction
    variable = mag2_x
    mag_sub = 1
    component = 0
  [../]
  [./afmdmi2_y]
    type = AFMSublatticeDMInteraction
    variable = mag2_y
    mag_sub = 1
    component = 1
  [../]
  [./afmdmi2_z]
    type = AFMSublatticeDMInteraction
    variable = mag2_z
    mag_sub = 1
    component = 2
  [../]

  [./afma1_x]
    type = AFMEasyPlaneAnisotropy
    variable = mag1_x
    mag_sub = 0
    component = 0
  [../]
  [./afma1_y]
    type = AFMEasyPlaneAnisotropy
    variable = mag1_y
    mag_sub = 0
    component = 1
  [../]
  [./afma1_z]
    type = AFMEasyPlaneAnisotropy
    variable = mag1_z
    mag_sub = 0
    component = 2
  [../]

  [./afma2_x]
    type = AFMEasyPlaneAnisotropy
    variable = mag2_x
    mag_sub = 1
    component = 0
  [../]
  [./afma2_y]
    type = AFMEasyPlaneAnisotropy
    variable = mag2_y
    mag_sub = 1
    component = 1
  [../]
  [./afma2_z]
    type = AFMEasyPlaneAnisotropy
    variable = mag2_z
    mag_sub = 1
    component = 2
  [../]

  [./afmsia1_x]
    type = AFMSingleIonCubicSixthAnisotropy
    variable = mag1_x
    mag_sub = 0
    component = 0
  [../]
  [./afmsia1_y]
    type = AFMSingleIonCubicSixthAnisotropy
    variable = mag1_y
    mag_sub = 0
    component = 1
  [../]
  [./afmsia1_z]
    type = AFMSingleIonCubicSixthAnisotropy
    variable = mag1_z
    mag_sub = 0
    component = 2
  [../]

  [./afmsia2_x]
    type = AFMSingleIonCubicSixthAnisotropy
    variable = mag2_x
    mag_sub = 1
    component = 0
  [../]
  [./afmsia2_y]
    type = AFMSingleIonCubicSixthAnisotropy
    variable = mag2_y
    mag_sub = 1
    component = 1
  [../]
  [./afmsia2_z]
    type = AFMSingleIonCubicSixthAnisotropy
    variable = mag2_z
    mag_sub = 1
    component = 2
  [../]

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


[AuxVariables]
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


  [./antiphase_A_x]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = antiphase_A_x
    initial_from_file_timestep = 'LATEST'
  [../]
  [./antiphase_A_y]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = antiphase_A_y
    initial_from_file_timestep = 'LATEST'
  [../]
  [./antiphase_A_z]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = antiphase_A_z
    initial_from_file_timestep = 'LATEST'
  [../]


  [./polar_x]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_x
    initial_from_file_timestep = 'LATEST'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_y
    initial_from_file_timestep = 'LATEST'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_z
    initial_from_file_timestep = 'LATEST'
  [../]

  [./ph]
    order = FIRST
    family = LAGRANGE
  [../]

  [./th1]
    order = FIRST
    family = LAGRANGE
  [../]
  [./th2]
    order = FIRST
    family = LAGRANGE
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

  [./phc]
    type = AngleBetweenTwoVectors
    variable = ph
    var1x = mag1_x
    var1y = mag1_y
    var1z = mag1_z
    var2x = mag2_x
    var2y = mag2_y
    var2z = mag2_z

    execute_on = 'initial timestep_end final'
  [../]

  [./th1c]
    type = AngleBetweenTwoVectors
    variable = th1
    var1x = mag1_x
    var1y = mag1_y
    var1z = mag1_z
    var2x = polar_x
    var2y = polar_y
    var2z = polar_z

    execute_on = 'initial timestep_end final'
  [../]

  [./th2c]
    type = AngleBetweenTwoVectors
    variable = th2
    var1x = mag2_x
    var1y = mag2_y
    var1z = mag2_z
    var2x = polar_x
    var2y = polar_y
    var2z = polar_z

    execute_on = 'initial timestep_end final'
  [../]

[]



[BCs]
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

  [./ph]
    type = ElementAverageValue
    execute_on = 'initial timestep_end final'
    variable = ph
  [../]
  [./th1]
    type = ElementAverageValue
    execute_on = 'initial timestep_end final'
    variable = th1
  [../]
  [./th2]
    type = ElementAverageValue
    execute_on = 'initial timestep_end final'
    variable = th2
  [../]

  [./FafmSLexch]
    type = AFMSublatticeSuperexchangeEnergy
    execute_on = 'initial timestep_end final'
    mag1_x = mag1_x
    mag1_y = mag1_y
    mag1_z = mag1_z
    mag2_x = mag2_x
    mag2_y = mag2_y
    mag2_z = mag2_z
    energy_scale = 1.0

  [../]
  [./FafmSLdmi]
    type = AFMSublatticeDMInteractionEnergy
    execute_on = 'initial timestep_end final'
    energy_scale = 1.0
  [../]

  [./Fllb1]
    type = MagneticExcessLLBEnergy
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    execute_on = 'initial timestep_end final'
  [../]
  [./Fllb2]
    type = MagneticExcessLLBEnergy
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    execute_on = 'initial timestep_end final'
  [../]

  [./Fa1]
    type = AFMEasyPlaneAnisotropyEnergy
    execute_on = 'initial timestep_end final'
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    energy_scale = 1.0
  [../]
  [./Fa2]
    type = AFMEasyPlaneAnisotropyEnergy
    execute_on = 'initial timestep_end final'
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    energy_scale = 1.0
  [../]

  [./Fsia1]
    type = AFMSingleIonCubicSixthAnisotropyEnergy
    execute_on = 'initial timestep_end final'
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    energy_scale = 1.0
  [../]
  [./Fsia2]
    type = AFMSingleIonCubicSixthAnisotropyEnergy
    execute_on = 'initial timestep_end final'
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    energy_scale = 1.0
  [../]

  [./Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'FafmSLexch FafmSLdmi Fa1 Fa2 Fsia1 Fsia2'
    pp_coefs = ' 1.0 1.0 1.0 1.0 1.0 1.0'
    execute_on = 'initial timestep_end final'
  [../]

  [./FtotLLB]
    type = LinearCombinationPostprocessor
    pp_names = 'Fllb1 Fllb2'
    pp_coefs = ' 1.0 1.0'
    execute_on = 'initial timestep_end final'
  [../]

  [./perc_change]
    type = EnergyRatePostprocessor
    postprocessor = Ftot
    dt = dt
    execute_on = 'timestep_end final'
  [../]


  [./elapsed]
    type = PerfGraphData
    section_name = "Root"  # for profiling the problem
    data_type = total
  [../]
[]

[UserObjects]
  [./kill]
    type = Terminator
    expression = 'perc_change <= 1.0e-5'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type '
    petsc_options_value = '    526               1e-8      1e-8      1e-8     bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'

  [./TimeIntegrator]
    type = NewmarkBeta
  [../]

  dtmin = 1e-18
  dtmax = 1.0e-8

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 10  #usually 10
    linear_iteration_ratio = 100
    dt = 1e-9
  [../]
  num_steps = 10000
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_BFO_P0A0_mRD
    elemental_as_nodal = true
    time_step_interval = 4
  [../]
[]

