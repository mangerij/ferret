
[Mesh]
  [./fmg]
    type = FileMeshGenerator
    file = sphere_tet_approx_size0_05.e
  []
[]


[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z
  potential_H_int = potential_H_int
[]

[Materials]
  [./constants] # Constants used in other material properties
    type = GenericConstantMaterial
    prop_names = ' alpha    Ae    Ms    g0     mu0   nx ny nz   long_susc t'
    prop_values = '0.01    0.013  1.2  176.086 1256.64   1  0  0          1.0     0'
  [../]
  [./a_long]
    type = GenericFunctionMaterial
    prop_names = 'alpha_long'
    prop_values = 'bc_func_1'
  [../]
  [./permitivitty_1]
    type = GenericConstantMaterial
    prop_names = 'permittivity'  #dummy variable at the moment since we use the "electrostatics" kernel
    prop_values = '1.0'
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
    value = 'st'   #*tanh(sl*t)+1.0'
    vars = 'st' #sl'
    vals = '100.0' # 795.775'
  [../]
[]

[Variables]
  [./mag_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi
      theta = polar_theta
      M0s = 1.0 #amplitude of the RandomConstrainedVectorFieldIC
      component  = 0
    [../]
  [../]
  [./mag_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi
      theta = polar_theta
      M0s = 1.0 #amplitude of the RandomConstrainedVectorFieldIC
      component  = 1
    [../]
  [../]
  [./mag_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi
      theta = polar_theta
      M0s = 1.0 #amplitude of the RandomConstrainedVectorFieldIC
      component  = 2
    [../]
  [../]
[]

[AuxVariables]


# Note that POTENTIAL is an AUX variable in the LLG master input file 
# but a VARIABLE in others

  [./potential_H_int]
    order = FIRST
    family = LAGRANGE
  [../]


  #--------------------------------------------#
  #                                            #
  #  field to seed IC that obeys constraint    #
  #                                            #
  #--------------------------------------------#
  [./azimuth_phi]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 1.1
      max = 1.12
      seed = 2
    [../]
  [../]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 1.4
      max = 1.41
      seed = 37
    [../]
  [../]

  [./mag_s]
    order = FIRST
    family = LAGRANGE
  [../]

  [./H_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./H_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./H_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./mag_mag]
    type = VectorMag
    variable = mag_s
    vector_x = mag_x
    vector_y = mag_y
    vector_z = mag_z
    execute_on = 'timestep_end final'
  [../]

  [./hxo]
    type = DemagFieldAux
    component = 0
    variable = H_x
    execute_on = 'timestep_end final'
  [../]
  [./hyo]
    type = DemagFieldAux
    component = 1
    variable = H_y
    execute_on = 'timestep_end final'
  [../]
  [./hzo]
    type = DemagFieldAux
    component = 2
    variable = H_z
    execute_on = 'timestep_end final'
  [../]
[]

[Kernels]
  #---------------------------------------#
  #                                       #
  #          Time dependence              #
  #                                       #
  #---------------------------------------#

  [./mag_x_time]
    type = TimeDerivative
    variable = mag_x
  [../]
  [./mag_y_time]
    type = TimeDerivative
    variable = mag_y
  [../]
  [./mag_z_time]
    type = TimeDerivative
    variable = mag_z
  [../]

  #---------------------------------------#
  #                                       #
  #          Magnetic exchange            #
  #                                       #
  #---------------------------------------#

  [./dllg_x_exch]
    type = ExchangeCartLL
    variable = mag_x
    component = 0
  [../]
  [./dllg_y_exch]
    type = ExchangeCartLL
    variable = mag_y
    component = 1
  [../]
  [./dllg_z_exch]
    type = ExchangeCartLL
    variable = mag_z
    component = 2
  [../]

  #---------------------------------------#
  #                                       #
  #         demagnetization field         #
  #                                       #
  #---------------------------------------#

  [./d_HM_x]
    type = InteractionCartLL
    variable = mag_x
    component = 0
  [../]
  [./d_HM_y]
    type = InteractionCartLL
    variable = mag_y
    component = 1
  [../]
  [./d_HM_z]
    type = InteractionCartLL
    variable = mag_z
    component = 2
  [../]

  #---------------------------------------#
  #                                       #
  #          LLB constraint terms         #
  #                                       #
  #---------------------------------------#

  [./llb_x]
    type = LongitudinalLLB
    variable = mag_x
    component = 0
  [../]
  [./llb_y]
    type = LongitudinalLLB
    variable = mag_y
    component = 1
  [../]

  [./llb_z]
    type = LongitudinalLLB
    variable = mag_z
    component = 2
  [../]
[]

[BCs]
  #---------------------------------------#
  #                                       #
  #  ground the magnetostatic potential   #
  #  at two boundaries                    #
  #                                       #
  #---------------------------------------#

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

  [./<M>]
    type = ElementAverageValue
    variable = mag_s
    execute_on = 'timestep_end final'
  [../]

  [./<mx>]
    type = ElementAverageValue
    variable = mag_x
    execute_on = 'timestep_end final'
  [../]
  [./<my>]
    type = ElementAverageValue
    variable = mag_y
    execute_on = 'timestep_end final'
  [../]
  [./<mz>]
    type = ElementAverageValue
    variable = mag_z
    execute_on = 'timestep_end final'
  [../]

  #---------------------------------------#
  #                                       #
  #   Calculate exchange energy of        #
  #   the magnetic body                   #
  #                                       #
  #---------------------------------------#

  [./Fexch]
    type = MagneticExchangeEnergy
    execute_on = 'timestep_end final'
  [../]

  #---------------------------------------#
  #                                       #
  #   Calculate demagnetization energy    #
  #   of the magnetic body                #
  #                                       #
  #---------------------------------------#

  [./Fdemag]
    type = MagnetostaticEnergyCart
    execute_on = 'timestep_end final'
  [../]

  #---------------------------------------#
  #                                       #
  #   Calculate excess energy from missed #
  #   LLB targets                         #
  #                                       #
  #---------------------------------------#

  [./Fllb]
    type = MagneticExcessLLBEnergy
    execute_on = 'timestep_end final'
  [../]


  #---------------------------------------#
  #                                       #
  #   add all the energy contributions    #
  #   and calculate their percent change  #
  #                                       #
  #---------------------------------------#

  [./Ftot]
    type = LinearCombinationPostprocessor 
    pp_names = 'Fexch Fdemag Fllb'
    pp_coefs = ' 1.0 1.0 1.0' 
    execute_on = 'timestep_end final'
  [../]

  [./perc_change]
    type = EnergyRatePostprocessor
    postprocessor = Ftot
    dt = dt 
    execute_on = 'timestep_end final'
  [../]
[]


[UserObjects]
  [./kill]
    type = Terminator
    expression = 'perc_change <= 1.0e-5'
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
    petsc_options = '-snes_ksp_ew -snes_converged_reason'
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121               1e-10      1e-8      1e-6      bjacobi'
  [../]
[]

#[Debug]
#   show_var_residual_norms = true
#[]

[MultiApps]
  [./poissonLLG]
    type = FullSolveMultiApp
    input_files = poissonLLG.i
    execute_on = initial
  [../]
  [./laplaceLLG]
    type = FullSolveMultiApp
    input_files = laplaceLLG.i
    execute_on = timestep_begin
  [../]
[]

[Transfers]



  [./M_from_sub]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = potential_H_int
    variable = potential_H_int
    multi_app = poissonLLG
  [../]


  [./pot_H_from_sub]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = potential_H_int
    variable = potential_H_int
    multi_app = poissonLLG
  [../]
  [./pot_H_to_sub]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    source_variable = potential_H_int
    variable = potential_H_int1
    multi_app = laplaceLLG
  [../]
  [./pot_H_from_sub2]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = potential_H_int
    variable = potential_H_int
    multi_app = laplaceLLG
  [../]
[]

[Executioner]
  type = Transient            
  solve_type = 'NEWTON'
  [./TimeIntegrator]
    type = ImplicitEuler
  [../]
  dtmin = 1e-9
  dtmax = 1.0e-3
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 18
    growth_factor = 1.3
    cutback_factor = 0.8
    dt = 3.0e-5
  [../]
  verbose = true
  num_steps = 800
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_LLG_ringdown
    elemental_as_nodal = true
    interval = 1
  [../]
  [./outCSV]
    type = CSV
    file_base = out_LLG_ringdown
  [../]
[]
