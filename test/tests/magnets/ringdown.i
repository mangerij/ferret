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

  [./bc_func_1]
    type = ParsedFunction
    expression = 'st'
    symbol_names = 'st'
    symbol_values = '1.e3'
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
      M0s = 1.0
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
[]

[AuxVariables]
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
    type = VectorMagnitudeAux
    variable = mag_s
    x = mag_x
    y = mag_y
    z = mag_z
    execute_on = 'initial timestep_end final'
    block = '1'
  [../]
[]

[Kernels]

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

  [./Fexch]
    type = MasterMagneticExchangeEnergy
    energy_scale = 1.
    execute_on = 'initial timestep_end final'
    block = '1'
  [../]

  [./Fdemag]
    type = MagnetostaticEnergyCart
    energy_scale = 1.
    execute_on = 'initial timestep_end final'
    block = '1'
  [../]

  [./Fllb1]
    type = MagneticExcessLLBEnergy
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
    execute_on = 'initial timestep_end final'
    block = '1'
  [../]

  [./Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'Fexch Fdemag'
    pp_coefs = ' 1.0 1.0'
    execute_on = 'initial timestep_end final'
  [../]

[]

[UserObjects]
  [mag]
    type = RenormalizeVector
    v = 'mag_x mag_y mag_z'
    norm = 1
    execute_on = 'TIMESTEP_END'
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

  [./TimeIntegrator]
    type = NewmarkBeta
  [../]

  dtmin = 1.e-5
  dtmax = 2.e-3
  end_time = 5.
  automatic_scaling = true

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 15
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
    time_step_interval = 10
  [../]
  [./outCSV]
    type = CSV
    file_base = ringdown
  [../]
[]
