
[Mesh]
  file = exodus_cylH_flat4_brick.e
[]

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z
  magnetic_x = mag_x
  magnetic_y = mag_y
  magnetic_z = mag_z

  potential_H_int = potential_H_int

  alpha = 0.4
  Ae = 0.013

  Ms = 1.0 #0.8

  M0s = 1.0 #amplitude of the RandomConstrainedVectorFieldIC

  g0 = 176.1

  permittivity = 1.0 #a dummy variable at the moment since we use the "electrostatics" kernel

  mu0 = 1256.1 

  norm = mag_s  # variable used for Norm Kernels

  alpha_long = 500.0

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
      min = 1.1
      max = 1.8
    [../]
  [../]

  [./mag_norm_x]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./mag_norm_y]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./mag_norm_z]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]


  [./mag_s]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]

 
  [./placer]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 0.9999999
      max = 1.0000001
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
    block = '1'
    execute_on = 'initial linear nonlinear timestep_end'
  [../]
[]

[Kernels]
   ## Time dependence
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


   #LLG simple

  # Exchange term

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

  # Magnetic interaction term

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

  # Magnetostatic Poisson equation

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
  # Ground the magnetostatic potential far from the ferromagnetic body.
  
  [./bc_int_pot_boundary]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = '1 2 3 4 5 6'
  [../]
[]

[Postprocessors]

  [./aveMs]
    type = ElementAverageValue
    variable = mag_s
    block = '1'
    execute_on = 'initial timestep_end'
  [../]


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
    petsc_options = '-snes_ksp_ew -snes_converged_reason'
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121               1e-8      1e-8      1e-6       lu'
  [../]
[]

[Debug]
  show_var_residual_norms = false
[]

[Executioner]
  type = Transient            
  solve_type = 'PJFNK'
  [./TimeIntegrator]
   type = Heun
  [../]
 # scheme = 'implicit-euler'   #, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-16
  dtmax = 1.0e-2
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 18
    growth_factor = 2.0
    cutback_factor = 0.6
    dt = 1.0e-7
  [../]
  verbose = true
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = ouLLG_Norm_n26
    interval = 5
    elemental_as_nodal = true
  [../]
[]
