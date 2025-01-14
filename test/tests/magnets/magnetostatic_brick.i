
[Mesh]
  [./fmg]
    type = FileMeshGenerator
    file = exodus_pysq_dx1_AnoI.e
  []

  [central_interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = fmg
    master_block = '1'
    paired_block = '2'
    primary_block = '1'
    new_boundary = '70'
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
    prop_values = '0.01    0.013  1.2  221010.0 1256.64   1  0  0          1.0     0'
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
    block = '1 2'
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
    block = '1'
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
    block = '1'
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
    block = '1'
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi
      theta = polar_theta
      M0s = 1.0 #amplitude of the RandomConstrainedVectorFieldIC
      component  = 2
    [../]
  [../]

  [./potential_H_int]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
  [../]
[]

[AuxVariables]
  #--------------------------------------------#
  #                                            #
  #  field to seed IC that obeys constraint    #
  #                                            #
  #--------------------------------------------#
  [./azimuth_phi]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 0.001
      max = 0.002
      seed = 2
    [../]
  [../]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 1.5707
      max = 1.5708
      seed = 37
    [../]
  [../]

  [./mag_s]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]

  [./H_x]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./H_y]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./H_z]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
[]

[AuxKernels]
  [./mag_mag]
    type = VectorMag
    variable = mag_s
    vector_x = mag_x
    vector_y = mag_y
    vector_z = mag_z
    block = '1'
    execute_on = 'timestep_end final'
  [../]

  [./hxo]
    type = DemagFieldAux
    component = 0
    variable = H_x
    block = '1 2'
    execute_on = 'timestep_end final'
  [../]
  [./hyo]
    type = DemagFieldAux
    component = 1
    variable = H_y
    block = '1 2'
    execute_on = 'timestep_end final'
  [../]
  [./hzo]
    type = DemagFieldAux
    component = 2
    variable = H_z
    block = '1 2'
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
    block = '1'
  [../]
  [./mag_y_time]
    type = TimeDerivative
    variable = mag_y
    block = '1'
  [../]
  [./mag_z_time]
    type = TimeDerivative
    variable = mag_z
    block = '1'
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
  #    Magnetostatic Poisson equation     #
  #                                       #
  #---------------------------------------#

  [./int_pot_lap]
    type = Electrostatics
    variable = potential_H_int
    block = '1 2'
  [../]
  [./int_bc_pot_lap]
    type = MagHStrongCart
    variable = potential_H_int
    block = '1'
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
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

  [./bc_int_pot_boundary]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '1 2 3 4 5 6'
  [../]

  [./bc_surface_mag_x]
    type = NeumannBC
    variable = mag_x
    value = 0.0
    boundary = '70'
  [../]
  [./bc_surface_mag_y]
    type = NeumannBC
    variable = mag_y
    value = 0.0
    boundary = '70'
  [../]
  [./bc_surface_mag_z]
    type = NeumannBC
    variable = mag_z
    value = 0.0
    boundary = '70'
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

  [./<M>]
    type = ElementAverageValue
    variable = mag_s
    block = '1'
    execute_on = 'timestep_end final'
  [../]

  [./<mx>]
    type = ElementAverageValue
    variable = mag_x
    block = '1'
    execute_on = 'timestep_end final'
  [../]
  [./<my>]
    type = ElementAverageValue
    variable = mag_y
    block = '1'
    execute_on = 'timestep_end final'
  [../]
  [./<mz>]
    type = ElementAverageValue
    variable = mag_z
    block = '1'
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
    block = '1'
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
    block = '1'
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
    block = '1'
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

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  [./TimeIntegrator]
    type = ImplicitEuler
  [../]
  dtmin = 1e-9
  dtmax = 1.0e-5
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 18
    growth_factor = 1.3
    cutback_factor = 0.8
    dt = 1.0e-8
  [../]
  verbose = true
  num_steps = 2
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_magnetostatic_brick
    elemental_as_nodal = true
    interval = 1
  [../]
  [./outCSV]
    type = CSV
    file_base = out_magnetostatic_brick
  [../]
[]
