Nx = 100
Ny = 25
Nz = 3

xMin = 0.0
yMin = 0.0
zMin = 1.5

xMax = 500.0
yMax = 125.0
zMax = 3.0

alphadef = 1.0

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

[UserObjects]
  [./reader_nearest]
    type = PropertyReadFile
    prop_file_name = 'fulldata_state.csv'
    read_type = 'voronoi'
    nprop = 6
    nblock = 1
    nvoronoi = 7500
  [../]
[]

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

  [./constants]
    type = GenericConstantMaterial
    prop_names = ' alpha           g0mu0Ms           permittivity Ae      Ms   '
    prop_values = '${alphadef}     34989.1        1.0        13.0   1.0  '
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
    symbol_values = '1e1'
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
    type = VectorMagnitudeAux
    variable = mag_s
    x = mag_x
    y = mag_y
    z = mag_z
    execute_on = 'initial timestep_end final'
  [../]
[]

[Kernels]

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

  [./Fexch]
    type = MasterMagneticExchangeEnergy
    energy_scale = 0.001
    execute_on = 'initial timestep_end final'
  [../]

  [./Fdemag]
    type = MagnetostaticEnergyCart
    energy_scale = 0.001
    execute_on = 'initial timestep_end final'
  [../]

  [./Fllb1]
    type = MagneticExcessLLBEnergy
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
    execute_on = 'initial timestep_end final'
  [../]

  [./Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'Fexch Fdemag'
    pp_coefs = ' 1.0 1.0'
    execute_on = 'initial timestep_end final'
  [../]

[]

[Preconditioning]

  [./smp]
    type = SMP
    full = true
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type '
    petsc_options_value = '    526               1e-5      1e-3      1e-6     bjacobi'
  [../]
[]
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
    optimal_iterations = 38
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
    time_step_interval = 1
  [../]
  [./outCSV]
    type = CSV
    file_base = out_Sstate_gen-Py_111_OOMF
  [../]
[]
