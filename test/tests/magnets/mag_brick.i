
[Mesh]
  [./fmg]
    type = FileMeshGenerator
    file = exodus_dx2_x10_y10_z30_b100.e
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
  [./constants]
    type = GenericConstantMaterial
    prop_names = ' alpha    Ae    Ms          g0          mu0         nx ny nz '
    prop_values = '0.01    0.013   1.2      221010.0   1256.64        0  1  1  '
  [../]

  [./a_long]
    type = GenericFunctionMaterial
    prop_names = 'alpha_long'
    prop_values = 'bc_func_1'
  [../]

  [./aniso]
    type = GenericConstantMaterial
    prop_names = 'K1 K2'
    prop_values = '20.0 0'        #  positive is unaxial
  [../]
  [./permitivitty_1]
    type = GenericConstantMaterial
    prop_names = 'permittivity'  # dummy variable at the moment since we use the "electrostatics" kernel
    prop_values = '1.0'
    block = '1 2'
  [../]
[]

[Functions]

  ###############################
  ##                           ##
  ## Define the function for   ##
  ##       alpha_long          ##
  ##                           ##
  ## here is just a (large)    ##
  ## constant                  ##
  ##                           ##
  ###############################

  [./bc_func_1]
    type = ParsedFunction
    expression = 'st'
    symbol_names = 'st'
    symbol_values = '1e3'
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
      M0s = 1.0
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
      M0s = 1.0
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
    [./InitialCondition]
      type = RandomIC
      min = 0.3
      max = 0.31
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

  [./azimuth_phi_mag_min1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 0.00001
      seed = 2
    [../]
  [../]
  [./polar_theta_mag_min1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 0.00001
      seed = 37
    [../]
  [../]
[]

[AuxKernels]

  #---------------------------------------#
  #                                       #
  #       compute magnitude of m          #
  #                                       #
  #---------------------------------------#

  [./mag_mag]
    type = VectorMag
    variable = mag_s
    vector_x = mag_x
    vector_y = mag_y
    vector_z = mag_z
    block = '1'
    execute_on = 'timestep_end final'
  [../]

  #---------------------------------------#
  #                                       #
  #    compute demag field grad*Phi       #
  #                                       #
  #---------------------------------------#

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
  #       Exchange stiffness              #
  #                                       #
  #---------------------------------------#

  [./dllg_x_exch]
    type = MasterExchangeCartLLG
    variable = mag_x
    component = 0
  [../]
  [./dllg_y_exch]
    type = MasterExchangeCartLLG
    variable = mag_y
    component = 1
  [../]
  [./dllg_z_exch]
    type = MasterExchangeCartLLG
    variable = mag_z
    component = 2
  [../]

  #---------------------------------------#
  #                                       #
  #         anisotropy                    #
  #                                       #
  #---------------------------------------#

  [./d_aM_x]
    type = MasterAnisotropyCartLLG
    variable = mag_x
    component = 0
  [../]
  [./d_aM_y]
    type = MasterAnisotropyCartLLG
    variable = mag_y
    component = 1
  [../]
  [./d_aM_z]
    type = MasterAnisotropyCartLLG
    variable = mag_z
    component = 2
  [../]

  #---------------------------------------#
  #                                       #
  #    demagnetization field              #
  #                                       #
  #---------------------------------------#

  [./d_HM_x]
    type = MasterInteractionCartLLG
    variable = mag_x
    component = 0
  [../]
  [./d_HM_y]
    type = MasterInteractionCartLLG
    variable = mag_y
    component = 1
  [../]
  [./d_HM_z]
    type = MasterInteractionCartLLG
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
  [../]

  #---------------------------------------#
  #                                       #
  #     LLB constraint terms              #
  #                                       #
  #---------------------------------------#

  [./llb_x]
    type = MasterLongitudinalLLB
    variable = mag_x
    component = 0
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
  [../]
  [./llb_y]
    type = MasterLongitudinalLLB
    variable = mag_y
    component = 1
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
  [../]

  [./llb_z]
    type = MasterLongitudinalLLB
    variable = mag_z
    component = 2
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
  [../]
[]

[BCs]
  #---------------------------------------#
  #                                       #
  #  ground the magnetostatic potential   #
  #  at boundaries of the bounding box    #
  #                                       #
  #---------------------------------------#

  [./bc_int_pot_boundary]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '1 2 3 4 5 6'
  [../]

  #---------------------------------------#
  #                                       #
  #  enforce Neumann condition (m*n = 0)  #
  #  at the boundary of the brick         #
  #                                       #
  #         Note: I don't think this      #
  #               does anything but we    #
  #               leave it in regardless  #
  #                                       #
  #---------------------------------------#

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

  #---------------------------------------#
  #                                       #
  #       track dt step size              #
  #                                       #
  #---------------------------------------#

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
    execute_on = 'initial timestep_end final'
  [../]

  [./<mx>]
    type = ElementAverageValue
    variable = mag_x
    block = '1'
    execute_on = 'initial timestep_end final'
  [../]
  [./<my>]
    type = ElementAverageValue
    variable = mag_y
    block = '1'
    execute_on = 'initial timestep_end final'
  [../]
  [./<mz>]
    type = ElementAverageValue
    variable = mag_z
    block = '1'
    execute_on = 'initial timestep_end final'
  [../]

  #---------------------------------------#
  #                                       #
  #   Calculate anisotropy energy of      #
  #   the magnetic body                   #
  #                                       #
  #---------------------------------------#

  [./Fa]
    type = MasterMagneticAnisotropyEnergy
    execute_on = 'initial timestep_end final'
    block = '1'
    energy_scale = 6241.51  #converts results to eV
  [../]

  #---------------------------------------#
  #                                       #
  #   Calculate exchange energy of        #
  #   the magnetic body                   #
  #                                       #
  #---------------------------------------#

  [./Fexch]
    type = MasterMagneticExchangeEnergy
    execute_on = 'timestep_end final'
    block = '1'
    energy_scale = 6241.51  #converts results to eV
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
    energy_scale = 6241.51  #converts results to eV
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
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
  [../]

  #---------------------------------------#
  #                                       #
  #   add all the energy contributions    #
  #   and calculate their percent change  #
  #                                       #
  #---------------------------------------#

  [./Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'Fexch Fdemag Fa'
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
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121               1e-8      1e-8      1e-5      bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  [./TimeIntegrator]
    type = ImplicitEuler
  [../]
  dtmin = 1e-12
  dtmax = 1.0e-2  #10 ns

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 10
    linear_iteration_ratio = 100
    dt = 1.0e-8   #10 fs

  [../]

  verbose = true

  num_steps = 2
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_mag_brick
    elemental_as_nodal = true
    time_step_interval = 1
  [../]
  [./outCSV]
    type = CSV
    file_base = out_mag_brick
  [../]
[]
