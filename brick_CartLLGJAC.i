
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 7
  ny = 7
  nz = 7
  xmin = -0.1
  xmax = 0.1
  ymin = -0.1
  ymax = 0.1
  zmin = -0.1
  zmax = 0.1
  elem_type = HEX8
[]

[GlobalParams]
 
  potential_H_int = potential_H_int
  potential_H_ext = potential_H_ext

  alpha = -0.85

  Ae = 0.013
  Ms = 0.8
  g0 = 176.1

  permittivity = 1.0
  mu0 = 1256.64
[]

#theta = 0  not a problem!!!

[Variables]
  [./lagrange]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 0.99999999
      max = 1.00000001
      seed = 3 
    [../]
  [../]

  [./mag_x]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      component  = 0
      phi = phi
      theta = theta
    [../]
  [../]
  [./mag_y]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      component  = 1
      phi = phi
      theta = theta
    [../]
  [../]
  [./mag_z]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      component  = 2
      phi = phi
      theta = theta
    [../]
  [../]
  [./potential_H_int]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
  [../]
  [./potential_H_ext]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
  [../]
[]

[AuxVariables]
  [./theta]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 0.00001
      max = 3.14159
      seed = 3 
    [../]
  [../]
  [./phi]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 0.00001
      max = 6.28319
      seed = 3 
    [../]
  [../]

[]



[AuxKernels]

[]

[Kernels]
  ## Time dependence
  [./mag_x_time]
    type = TimeDerivativeScaled
    variable = mag_x
    time_scale = 1.0
    block = '1'
  [../]
  [./mag_y_time]
    type = TimeDerivativeScaled
    variable = mag_y
    time_scale = 1.0
    block = '1'
  [../]

  [./mag_z_time]
    type = TimeDerivativeScaled
    variable = mag_z
    time_scale = 1.0
    block = '1'
  [../]


   #LLG simple

  # Exchange term

  [./d_llg_exch_x]
    type = ExchangeCartLLG
    variable = mag_x
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
    component = 0
  [../]
  [./d_llg_exch_y]
    type = ExchangeCartLLG
    variable = mag_y
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
    component = 1
  [../]
  [./d_llg_exch_z]
    type = ExchangeCartLLG
    variable = mag_z
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
    component = 2
  [../]

  # Anisotropy term. TURNED OFF

  # Magnetic interaction term

  [./d_HM_x]
    type = InteractionCartLLG
    variable = mag_x
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
    component = 0
  [../]
  [./d_HM_y]
    type = InteractionCartLLG
    variable = mag_y
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
    component = 1
  [../]
  [./d_HM_z]
    type = InteractionCartLLG
    variable = mag_z
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
    component = 2
  [../]

  # Magnetostatic Poisson equation

  [./ext_pot_lap]
    type = Electrostatics
    variable = potential_H_ext
    block = '1 2'
  [../]
  [./int_pot_lap]
    type = Electrostatics
    variable = potential_H_int
    block = '1 2'
  [../]
  [./int_bc_pot_lap]
    type = MagHStrongCart
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
    variable = potential_H_int
    block = '1'
  [../]

  [./l_mult]
    type = LagrangeConstraintCartLLG
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
    variable = lagrange
    block = '1'
  [../]
[]

[BCs]
   [./bc_int_pot_R]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = '1'
  [../]

  [./bc_int_pot_L]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = '2'
  [../]

  [./bc_int_pot_T]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = '3'
  [../]
  [./bc_int_pot_Bo]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = '4'
  [../]

  [./bc_int_pot_F]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = '5'
  [../]
  [./bc_int_pot_Ba]
    type = PresetBC
    variable = potential_H_int
    value = 0.0
    boundary = '6'
  [../]

  [./bc_ext_pot_front]
    type = PresetBC
    variable = potential_H_ext
    value = 0.0
    boundary = '1'
  [../]
  [./bc_ext_pot_back]
    type = PresetBC
    variable = potential_H_ext
    value = 0.0
    boundary = '2'
  [../]
[]

[Postprocessors]
  [./Fexchange]
    type = MagneticExchangeEnergy
    execute_on = 'initial timestep_end'
    block = '1'
    magnetic_x = mag_x
    magnetic_y = mag_y
    magnetic_z = mag_z
  [../]
  [./Fdemag]
    type = MagnetostaticEnergyCart
    execute_on = 'initial timestep_end'
    magnetic_x = mag_x
    magnetic_y = mag_y
    magnetic_z = mag_z
    block = '1'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor 
    pp_names = 'Fexchange Fdemag'
    pp_coefs = ' 1 1 ' 
    execute_on = 'initial timestep_end'
  [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason'
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type '
    petsc_options_value = '    121               1e-10      1e-8      1e-8       lu'
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[Executioner]
  type = Transient            
  solve_type = 'NEWTON'
  scheme = 'implicit-euler'   #, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-16
  dtmax = 5.0e-5
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 12
    growth_factor = 2.0
    cutback_factor = 0.4
    dt = 1.0e-8
  [../]
  verbose = true
[]

[Outputs]
  print_linear_residuals = true
  [./out]
    type = Exodus
    file_base = outCartLLG_test_VB
    interval = 1
    elemental_as_nodal = true
  [../]
[]
