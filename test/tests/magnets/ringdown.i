

  ##############################
  ##
  ## UNITS:
  ##    
  ##   gamma = (2.2101*10^5 ) m/C
  ##
  ## NOTE:
  ##   gamma*Hscale = 1/ns
  ##
  ##   coefficients given in (pg/nm*ns) 
  ##   which is equivalent to an energy/vol
  ##
  ##   Effective fields are 1/(mu0*Ms)*coeff
  ##   which gives units of aC/(nm*mus)
  ##
  ##   Energies are natively printed in units 
  ##    of 0.160218 pg nm^2 / mus^2 
  ##    or 1.60218*10^{-22} J
  ##    or 0.001 eV
  ##
  ##############################

################
#
#  LLG alpha:
#
################

alphadef = 0.02


[Mesh]
  [./mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 50
    ny = 50
    nz = 10
    xmin = -50
    xmax = 50
    ymin = -50
    ymax = 50
    zmin = -10
    zmax = 10
  [../]  
  [./vacuum_box]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    bottom_left = '-50 -50 -10'
    top_right = '50 50 10'
    block_id = 2
    block_name = vacuum 
  [../]
  [./brick]
    type = SubdomainBoundingBoxGenerator
    input = vacuum_box
    bottom_left = '-10 -10 -1.5'
    top_right = '10 10 1.5'
    block_id = 1
    block_name = brick
  [../]
[../]


[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z

  potential_H_int = potential_H_int
    
  Hscale = 0.004519239
  g0 = 1.0
  mu0 = 1.256637e-06

[]

[Materials]
  ############################################################################
  ##
  ##       material constants used.
  ##
  ############################################################################

  [./constants] 
    type = GenericConstantMaterial
    prop_names = ' alpha                 Ae      Ms   permittivity'
    prop_values = '${alphadef}          1.3e-05  1.2  1.'
    block = '1'
  [../]

  [./a_long]
    type = GenericFunctionMaterial
    prop_names = 'alpha_long'
    prop_values = 'bc_func_1'
    block = '1' 
  [../]
 [./constantsv]
    type = GenericConstantMaterial
    prop_names = ' alpha                Ae      Ms  permittivity'
    prop_values = '1                   1.e-05   0.  1. '
    block = '2'
  [../]
  [./a_longv]
    type = GenericFunctionMaterial
    prop_names = 'alpha_long'
    prop_values = 'bc_func_1'
    block = '2'
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
    value = 'st'
    vars = 'st'
    vals = '1.e3'  #3?
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
    block= '1'
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
    [./InitialCondition]
      type = RandomIC
      min = 0.
      max = 0.01
      seed = 2
    [../]
  [../]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
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
  [../]

  [./H_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  [./H_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
  [../]
  [./H_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 1.0
    [../]
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
    block = '1'
  [../]
  [./dllg_y_exch]
    type = MasterExchangeCartLLG
    variable = mag_y
    component = 1
    block = '1'
  [../]
  [./dllg_z_exch]
    type = MasterExchangeCartLLG
    variable = mag_z
    component = 2
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

#  [./llb1_x]
#    type = MasterLongitudinalLLB
#    variable = mag_x
#    component = 0
#    block = '1'
#  [../]
#  [./llb1_y]
#    type = MasterLongitudinalLLB
#    variable = mag_y
#    component = 1
#    block = '1'
#  [../]
#
#  [./llb1_z]
#    type = MasterLongitudinalLLB
#    variable = mag_z
#    component = 2
#    block = '1'
#  [../]

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

[]

[BCs]
  [./vacuum_box]
    type = DirichletBC
    value = 0.
    variable = potential_H_int
    boundary = '0 1 2 3 4 5'
  [../]
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
    energy_scale = 1.  #converts results to eV
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
    energy_scale = 1.  #converts results to eV
    execute_on = 'initial timestep_end final'
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
  [mag]
    type = RenormalizeVector
    v = 'mag_x mag_y mag_z'
    norm = 1
    execute_on = 'TIMESTEP_END'
#    force_praux = true
  []
[]

[Preconditioning]
    active = smp
 [./muPBP]
    type = PBP
    solve_order = 'mag_x mag_y mag_z potential_H_int'
    preconditioner = 'AMG ILU'
    off_diag_row = 'mag_x mag_y mag_z'
    off_diag_column = 'mag_x mag_y mag_z'
 [../]
 [./muFSP]
    type = FSP
    topsplit = 'magpot'
    [./magpot]
      splitting = 'mag pot'
      splitting_type = additive
    [../]
    [./mag]
      vars = 'mag_x mag_y mag_z'
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type -pc_sub_type '
    petsc_options_value = '    40               1e-20      1e-6      1e-6     bjacobi  ilu'
    [../]
    [./pot]
      vars = potential_H_int
      petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type -pc_sub_type '
      petsc_options_value = '    40               1e-12      1e-6      1e-6     bjacobi ilu'
    [../]
  [../]

  #---------------------------------------#
  #                                       #
  #            Solver options             #
  #                                       #
  #---------------------------------------#

  [./smp]
    type = SMP
    full = true
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    40               1e-8      1e-8      1e-4      bjacobi' 
  [../]
[]

[Executioner]
  type = Transient            
  solve_type = 'NEWTON'
#  num_grids = 8

  [./TimeIntegrator]
    type = NewmarkBeta #LStableDirk4 #NewmarkBeta
  [../]

  dtmin = 1.e-5
  dtmax = 2.e-3
  end_time = 5.
  automatic_scaling = true

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 15  #usually 10
    iteration_window = 2
    linear_iteration_ratio = 1000
    dt = 1.e-4
    growth_factor = 1.1
    cutback_factor = 0.75
  [../]
  num_steps = 20
[../]



[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = ringdown
    elemental_as_nodal = true
    interval = 10
  [../]
  [./outCSV]
    type = CSV
    file_base = ringdown
  [../]
[]
