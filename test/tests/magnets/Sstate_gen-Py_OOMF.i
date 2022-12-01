

  ##############################
  ##
  ## UNITS:
  ##    
  ##   gamma = (2.2101*10^5 / 2 pi) m/C
  ##
  ##   gamma/(mu*0Ms) = 48291.9 nm*mus/pg
  ##
  ##   coefficients given in (pg/nm*mus) 
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


Nx = 100
Ny = 25
Nz = 3

#

xMin = 0.0
yMin = 0.0
zMin = 1.5

xMax = 500.0
yMax = 125.0
zMax = 3.0

alphadef = 1.0

#
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${Nx}
    ny = ${Ny}
    nz = ${Nz}
    xmin = ${xMin}
    xmax = ${xMax}
    ymin = ${yMin}
    ymax = ${yMax}
    zmin = ${zMin}
    zmax = ${zMax}
    elem_type = HEX8
  []
[]

#

[UserObjects]
  [./kill]
    type = Terminator
    expression = 'perc_change <= 1.0e-8'
  [../]
  [./reader_nearest]
    type = PropertyReadFile
    prop_file_name = 'fulldata_state.csv'
    read_type = 'voronoi'
    nprop = 6  # number of columns in CSV
    nblock = 1
    nvoronoi = 7500
    #use_random_voronoi = true
  [../]
[]
#

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z

  potential_H_int = potential_H_int


  mu0 = 1.0
  Hscale = 1.0
  g0 = 17680.8
[]

[Materials]
  ############################################################################
  ##
  ##      material constants used.
  ##
  ##
  ############################################################################

  [./constants] 
    type = GenericConstantMaterial
    prop_names = ' alpha           g0mu0Ms           permittivity Ae      Ms   '
    prop_values = '${alphadef}     34989.1        1.0        13.0   1.0  '
  [../]

#NOTE: g0 is g*mu0*Ms/2 as defined by Hertel 
#alpha is chosen to be 1.0 as in the muMag paper

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
    value = 'st'
    vars = 'st'
    vals = '1e1'  #3?
  [../]

  [./node_mx]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = 'reader_nearest'
    read_type = 'voronoi'
    column_number = '3'
  [../]
  [./node_my]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = 'reader_nearest'
    read_type = 'voronoi'
    column_number = '4'
  [../]
  [./node_mz]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = 'reader_nearest'
    read_type = 'voronoi'
    column_number = '5'
  [../]
[]
#

[Variables]
  [./mag_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = node_mx
    [../]
  [../]
  [./mag_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = node_my
    [../]
  [../]
  [./mag_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = node_mz
    [../]
  [../]

  [./potential_H_int]
    order = FIRST
    family = LAGRANGE
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
      min = 0.234979
      max = 0.244979
      seed = 2
    [../]
  [../]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 1.46742
      max = 1.47742
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
  #    Local magnetic exchange            #
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
  #          LLB constraint terms         #
  #                                       #
  #---------------------------------------#

  [./llb1_x]
    type = LongitudinalLLB
    variable = mag_x
    component = 0
  [../]
  [./llb1_y]
    type = LongitudinalLLB
    variable = mag_y
    component = 1
  [../]

  [./llb1_z]
    type = LongitudinalLLB
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
  [../]
  [./int_bc_pot_lap]
    type = MagHStrongCart
    variable = potential_H_int
  [../]

[]

[BCs]

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
  [../]

  [./<mx>]
    type = ElementAverageValue
    variable = mag_x
    execute_on = 'initial timestep_end final'
  [../]
  [./<my>]
    type = ElementAverageValue
    variable = mag_y
    execute_on = 'initial timestep_end final'
  [../]
  [./<mz>]
    type = ElementAverageValue
    variable = mag_z
    execute_on = 'initial timestep_end final'
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
    execute_on = 'initial timestep_end final'
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

[Preconditioning]

  #---------------------------------------#
  #                                       #
  #            Solver options             #
  #                                       #
  #---------------------------------------#

  [./smp]
    type = SMP
    full = true
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type '
    petsc_options_value = '    526               1e-5      1e-3      1e-6     bjacobi'
  [../]
[]
#
[Executioner]
  type = Transient            
  solve_type = 'NEWTON'

  [./TimeIntegrator]
    type = NewmarkBeta
  [../]

  dtmin = 1e-18
  dtmax = 1.0

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 38  #usually 10
    linear_iteration_ratio = 1000
    dt = 1e-8
    growth_factor = 1.1
  [../]

  num_steps = 2
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_Sstate_gen-Py_111_OOMF
    elemental_as_nodal = true
    interval = 1
  [../]
  [./outCSV]
    type = CSV
    file_base = out_Sstate_gen-Py_111_OOMF
  [../]
[]
