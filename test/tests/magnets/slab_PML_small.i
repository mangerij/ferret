
[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z

  potential_H_int = potential_H_int

  g0 = 1.0
  Hscale = 0.004519239
  mu0 = 1.256637e-06
  y0pmlminus = -25
  deltasyminus = 1.0e04
  deltawyminus = 5
  deltapyminus = 0.1
[]

alphadef = 0.02

[Mesh]
  [./mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 25
    ny = 43
    nz = 5
    xmin = -50
    xmax = 50
    ymin = -35
    ymax = 50
    zmin = -10
    zmax = 10
  [../]
  [./infinite_domain]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    bottom_left = '-50 -35 -10'
    top_right = '50 50 10'
    block_id = 3
    block_name = infinite_domain
  [../]
  [./vacuum_box]
    type = SubdomainBoundingBoxGenerator
    input = infinite_domain
    bottom_left = '-50 -25 -10'
    top_right = '50 50 10'
    block_id = 2
    block_name = vacuum
  [../]
  [./boundary1]
    type = SideSetsBetweenSubdomainsGenerator
    input = vacuum_box
    new_boundary = id_boundary
    paired_block = vacuum
    primary_block = 3 #mesh #infinite_domain
  [../]
  [./brick]
    type = SubdomainBoundingBoxGenerator
    input = boundary1 #vacuum_box
    bottom_left = '-10 -10 -1.5'
    top_right = '10 10 1.5'
    block_id = 1
    block_name = brick
  [../]
[../]



[Materials]
  ############################################################################
  ##
  ##       material constants used.
  ##
  ##
  ############################################################################

  [./constants]
    type = GenericConstantMaterial
    prop_names = ' alpha           permittivity  Ae        Ms'
    prop_values = '${alphadef}     1.0           1.3e-05   1.2'
    block = '1'
  [../]

#NOTE: g0 is g*mu0*Ms/2 as defined by Hertel
#alpha is chosen to be 1.0 as in the muMag paper

 [./a_long]
    type = GenericFunctionMaterial
    prop_names = 'alpha_long'
    prop_values = 'bc_func_1'
    block = '1'
 [../]
 [./constantsv]
    type = GenericConstantMaterial
    prop_names = ' permittivity'
    prop_values = '1.0'
    block = '3 2'
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
    vals = '1.e1'  #3?
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
    block= '1'
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
  [./phi2]
   order = FIRST
   family = LAGRANGE
   block = 3
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
      min = 1.5708
      max = 1.5709
      seed = 2
    [../]
  [../]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 1.5708
      max = 1.5709
      seed = 37
    [../]
  [../]

  [./mag_s]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./H_y]
    order = FIRST
    family = MONOMIAL
    block = '1 2'
  [../]
  [./H_y_v]
    order = FIRST
    family = MONOMIAL
    block = '3'
  [../]
[]


[AuxKernels]
  [./mag_mag]
    type = VectorMag
    variable = mag_s
    vector_x = mag_x
    vector_y = mag_y
    vector_z = mag_z
    execute_on = 'initial timestep_end final'
    block = '1'
  [../]
  [./Hy_brick]
    type = DemagFieldAux
    variable = H_y
    component = 1
    block = '1 2'
  [../]
  [./Hy_id]
    type = DemagFieldAuxPML
    phi1 = phi2
    variable = H_y_v
    component = 1
    block = '3'
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
  #    Local magnetic exchange            #
  #                                       #
  #---------------------------------------#

  [./dllg_x_exch]
    type = MasterExchangeCartLLG
    variable = mag_x
    component = 0
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
    block = '1'
  [../]
  [./dllg_y_exch]
    type = MasterExchangeCartLLG
    variable = mag_y
    component = 1
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
    block = '1'
  [../]
  [./dllg_z_exch]
    type = MasterExchangeCartLLG
    variable = mag_z
    component = 2
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
    block = '1'
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
    block = '1'
  [../]
  [./d_HM_y]
    type = MasterInteractionCartLLG
    variable = mag_y
    component = 1
    block = '1'
  [../]
  [./d_HM_z]
    type = MasterInteractionCartLLG
    variable = mag_z
    component = 2
    block = '1'
  [../]

  #---------------------------------------#
  #                                       #
  #          LLB constraint terms         #
  #                                       #
  #---------------------------------------#

  [./llb1_x]
    type = MasterLongitudinalLLB
    variable = mag_x
    component = 0
    block = '1'
  [../]
  [./llb1_y]
    type = MasterLongitudinalLLB
    variable = mag_y
    component = 1
    block = '1'
  [../]

  [./llb1_z]
    type = MasterLongitudinalLLB
    variable = mag_z
    component = 2
    block = '1'
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
  [./infinite_domain]
    type = MagneticPMLCart
    variable = phi2
    component = 1
    block = '3'
  [../]

[]

[InterfaceKernels]
  [./PML_boundary]
    type = InterfaceEquality
    variable = potential_H_int
    boundary = id_boundary
    neighbor_var = potential_H_int
    permittivity_neighbor = 1.
  [../]
[]

[BCs]
  [./vacuum_box]
    type = DirichletBC
    variable = phi2
    value = 0.
    boundary = 'bottom'
  [../]
  [./other_boundaries]
    type = DirichletBC
    variable = potential_H_int
    value = 0.
    boundary = 'left right front back top '
  [../]
  [./boundarydirichlet]
    type = CoupledDirichletBC
    variable = phi2
    coupled_var = potential_H_int
    boundary = id_boundary
  []
[]

[Postprocessors]
   [./dt]
     type = TimestepSize
   [../]

  #---------------------------------------#
  #                                       #
  #     Average M = |m|                   #
  #                                       #
  #---------------------------------------#

  [./M1]
    type = ElementAverageValue
    variable = mag_s
    execute_on = 'initial timestep_end final'
    block = '1'
  [../]

  [./<mx>]
    type = ElementAverageValue
    variable = mag_x
    execute_on = 'initial timestep_end final'
    block = '1'
  [../]
  [./<my>]
    type = ElementAverageValue
    variable = mag_y
    execute_on = 'initial timestep_end final'
    block = '1'
  [../]
  [./<mz>]
    type = ElementAverageValue
    variable = mag_z
    execute_on = 'initial timestep_end final'
    block = '1'
  [../]

  #---------------------------------------#
  #                                       #
  #   Calculate exchange energy of        #
  #   the magnetic body                   #
  #                                       #
  #---------------------------------------#

  [./Fexch]
    type = MasterMagneticExchangeEnergy
    energy_scale = 0.001  #converts results to eV
    execute_on = 'initial timestep_end final'
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
    energy_scale = 0.001  #converts results to eV
    execute_on = 'initial timestep_begin timestep_end final'
    block = '1'
  [../]


  #---------------------------------------#
  #                                       #
  #   Calculate excess energy from missed #
  #   LLB targets                         #
  #                                       #
  #---------------------------------------#

  [./Fllb1]
    type = MagneticExcessLLBEnergy
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
    execute_on = 'initial timestep_end final'
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
    pp_names = 'Fexch Fdemag'
    pp_coefs = ' 1.0 1.0'
    execute_on = 'initial timestep_end final'
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
    expression = 'perc_change <= 1.0e-8'
  [../]
[]


[Preconditioning]

  #---------------------------------------#
  #                                       #
  #            Solver options             #
  #                                       #
  #---------------------------------------#

  [./smp]
    type = SMP #FDP
#     petsc_options_iname = ' -pc_type -mat_fd_coloring_err -mat_fd_type'
#     petsc_options_value = 'lu   1.e-06   ds'
    full = true
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type -sub_pc_type '
    petsc_options_value = '    100               1e-12      1e-9      1e-8     bjacobi   ilu'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  automatic_scaling = true

  [./TimeIntegrator]
    type = NewmarkBeta
  [../]

  dtmin = 1.e-5
  dtmax = 1.e-3
  end_time = 2

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 25  #usually 10
    linear_iteration_ratio = 100
    dt = 1e-4
    growth_factor = 1.1
    cutback_factor = 0.75
  [../]
  num_steps = 2
[../]



[Outputs]
  print_linear_residuals = false
  [pgraph]
  type = PerfGraphOutput
  execute_on = final
  level = 2
[]
  [./out]
    type = Exodus
    file_base = slab_PML_small_25_002
    elemental_as_nodal = true
    time_step_interval = 1
    execute_on = 'initial timestep_end'
  [../]
  [./outCSV]
    type = CSV
    file_base = slab_PML_small_25_002
    time_step_interval = 1
  [../]
[]
