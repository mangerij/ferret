
[Mesh]
  file = exodus_cyl_flat4_brick.e
[]

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z

  lambda = lambda

  potential_H_int = potential_H_int

  time_scale = 1.0
  alpha = 0.85
  Ae = 0.02
  Ms = 0.8
  M0s = 1.0

  g0 = 176.0

  mu0 = 1254.5

  permittivity = 1.0 #a dummy variable at the moment since we use the "electrostatics" kernel

  norm = mag_s

  eps = 1e-8 #this should be a local quantity depending on the deviation at the specific mesh point (not a global quantity).
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
      min = 0.0
      max = 6.283185307178
      seed = 3
    [../]
  [../]

  [./lambda]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 10000.0
      max = 10001.0
      seed = 3
    [../]
    scaling = 1e7
  [../]
[]

[AuxVariables]
  [./azimuth_phi]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 6.28318
    [../]
  [../]
  [./polar_theta]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 3.14159
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
    type = TimeLLG
    variable = mag_x
    component  = 0
    block = '1'
  [../]
  [./mag_y_time]
    type = TimeLLG
    variable = mag_y
    component  = 1
    block = '1'
  [../]
  [./mag_z_time]
    type = TimeLLG
    variable = mag_z
    component  = 2
    block = '1'
  [../]


   #LLG simple

  # Exchange term

  [./dllg_x_exch]
    type = RHSExchangeCartLLG
    variable = mag_x
    component = 0
  [../]
  [./dllg_y_exch]
    type = RHSExchangeCartLLG
    variable = mag_y
    component = 1
  [../]
  [./dllg_z_exch]
    type = RHSExchangeCartLLG
    variable = mag_z
    component = 2
  [../]

  # Magnetic interaction term

  [./d_HM_x]
    type = RHSInteractionCartLLG
    variable = mag_x
    component = 0
  [../]
  [./d_HM_y]
    type = RHSInteractionCartLLG
    variable = mag_y
    component = 1
  [../]
  [./d_HM_z]
    type = RHSInteractionCartLLG
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


  # Constraints

  [./d_LM_x]
    type = MagConstraintCartLLG
    variable = mag_x
    component = 0
  [../]
  [./d_LM_y]
    type = MagConstraintCartLLG
    variable = mag_y
    component = 1
  [../]
  [./d_LM_z]
    type = MagConstraintCartLLG
    variable = mag_z
    component = 2
  [../]

  [./d_LM_L]
    type = LagrangeConstraintCartLLG
    variable = lambda
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
#  [./Fdemag]
#    type = MagnetostaticEnergy
#    execute_on = 'initial timestep_end'
#    block = '1'
#  [../]
#  [./Ftot]
#    type = LinearCombinationPostprocessor 
#    pp_names = 'Fexchange Fdemag'
#    pp_coefs = ' 1 1 ' 
#    execute_on = 'initial timestep_end'
#  [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason'
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121               1e-10      1e-8      1e-6       lu'
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
  dtmax = 1.0
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 12
    growth_factor = 2.0
    cutback_factor = 0.4
    dt = 1.0e-10
  [../]
  verbose = true
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = ouLLG_Constr_2
    interval = 1
    elemental_as_nodal = true
  [../]
[]
