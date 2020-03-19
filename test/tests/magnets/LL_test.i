
[Mesh]
  file = exodus_cylS_flat4_brick.e
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
    prop_names = ' alpha    Ae    Ms g0 mu0 alpha_long'
    prop_values = '0.4  0.013  1.0 0.1761 1256.1 500'
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
    [./InitialCondition]
      type = RandomIC
      min = -5.0
      max = 5.0
    [../]
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
      min = 1.3
      max = 3.0
    [../]
  [../]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 0.
      max = 0.2
    [../]
  [../]

  [./mag_s]
    order = FIRST
    family = LAGRANGE
    block = '1'
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
    execute_on = 'initial linear nonlinear timestep_end'
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
    permittivity = 1.0 #a dummy variable at the moment since we use the "electrostatics" kernel
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
  #          LLB constraint term          #
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
  #  far from the ferromagnetic body      #
  #                                       #
  #---------------------------------------#
  
  [./bc_int_pot_boundary]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '1 2 3 4 5 6'
  [../]
[]

[Postprocessors]
  #---------------------------------------#
  #                                       #
  #             Average |M|               #
  #                                       #
  #---------------------------------------#

  [./<M>]
    type = ElementAverageValue
    variable = mag_s
    block = '1'
    execute_on = 'initial timestep_end'
  [../]

  #---------------------------------------#
  #                                       #
  #   Calculate exchange energy of        #
  #   the magnetic body                   #
  #                                       #
  #---------------------------------------#

  [./Fexch]
    type = MagneticExchangeEnergy
    execute_on = 'initial timestep_end'
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
    execute_on = 'initial timestep_end'
    block = '1'
  [../]

  [./Ftot]
    type = LinearCombinationPostprocessor 
    pp_names = 'Fexch Fdemag'
    pp_coefs = ' 1 1 ' 
    execute_on = 'initial timestep_end'
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
    petsc_options_value = '    121               1e-8      1e-8      1e-6       lu'
  [../]
[]

[Executioner]
  type = Transient            
  solve_type = 'NEWTON' #'NEWTON' #'PJFNK'
  [./TimeIntegrator]
   type = ImplicitEuler #Heun
  [../]
 # scheme = 'implicit-euler'   #, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-6
  dtmax = 5.0e-3
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 18
    growth_factor = 1.3
    cutback_factor = 0.8
    dt = 1.0e-7
  [../]
  verbose = true
  num_steps = 10
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_test
    interval = 5
    elemental_as_nodal = true
  [../]
[]
