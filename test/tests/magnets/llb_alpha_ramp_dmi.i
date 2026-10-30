[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 600
  xmin = -30
  xmax = 30
[]

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z
  g0 = 1
  mu0 = 1
  Hscale = 1
  llb = true
[]

[Functions]
  [ap]
    type = ParsedFunction
    expression = '0.25 + -0.002*x'
  []
  [al]
    type = ParsedFunction
    expression = '0.15 + 0.002*x'
  []
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'alpha Ms Ae K1 nx ny nz me chi_par D'
    prop_values = '0.1 1 1 1 0 0 1 0.8 0.05 0.3'
  []
  [alphas]
    type = GenericFunctionMaterial
    prop_names = 'alpha_perp alpha_par'
    prop_values = 'ap al'
  []
  [grad_alphas]
    type = GenericConstantVectorMaterial
    prop_names = 'grad_alpha_perp grad_alpha_par'
    prop_values = '-0.002 0 0 0.002 0 0'
  []
[]

[Variables]
  [mag_x]
    [InitialCondition]
      type = FunctionIC
      function = '0.8/cosh(x)'
    []
  []
  [mag_y]
    [InitialCondition]
      type = FunctionIC
      function = '0'
    []
  []
  [mag_z]
    [InitialCondition]
      type = FunctionIC
      function = '0.8*tanh(x)'
    []
  []
[]

[AuxVariables]
  [uz]
  []
[]

[AuxKernels]
  [uz]
    type = ParsedAux
    variable = uz
    coupled_variables = 'mag_x mag_y mag_z'
    expression = 'mag_z/sqrt(mag_x^2 + mag_y^2 + mag_z^2)'
  []
[]

[Kernels]
  [time_x]
    type = TimeDerivative
    variable = mag_x
  []
  [time_y]
    type = TimeDerivative
    variable = mag_y
  []
  [time_z]
    type = TimeDerivative
    variable = mag_z
  []
  [MasterExchangeCartLLG_x]
    type = MasterExchangeCartLLG
    variable = mag_x
    component = 0
  []
  [MasterExchangeCartLLG_y]
    type = MasterExchangeCartLLG
    variable = mag_y
    component = 1
  []
  [MasterExchangeCartLLG_z]
    type = MasterExchangeCartLLG
    variable = mag_z
    component = 2
  []
  [MasterAnisotropyCartLLG_x]
    type = MasterAnisotropyCartLLG
    variable = mag_x
    component = 0
  []
  [MasterAnisotropyCartLLG_y]
    type = MasterAnisotropyCartLLG
    variable = mag_y
    component = 1
  []
  [MasterAnisotropyCartLLG_z]
    type = MasterAnisotropyCartLLG
    variable = mag_z
    component = 2
  []
  [MasterDMICartLLG_x]
    type = MasterDMICartLLG
    variable = mag_x
    component = 0
    dmi_type = interfacial
  []
  [MasterDMICartLLG_y]
    type = MasterDMICartLLG
    variable = mag_y
    component = 1
    dmi_type = interfacial
  []
  [MasterDMICartLLG_z]
    type = MasterDMICartLLG
    variable = mag_z
    component = 2
    dmi_type = interfacial
  []
  [MasterLongitudinalLLB_x]
    type = MasterLongitudinalLLB
    variable = mag_x
    component = 0
    form = llb
    temperature = 0.5
    Tc = 1
  []
  [MasterLongitudinalLLB_y]
    type = MasterLongitudinalLLB
    variable = mag_y
    component = 1
    form = llb
    temperature = 0.5
    Tc = 1
  []
  [MasterLongitudinalLLB_z]
    type = MasterLongitudinalLLB
    variable = mag_z
    component = 2
    form = llb
    temperature = 0.5
    Tc = 1
  []
[]

[UserObjects]
  [assert]
    type = Terminator
    expression = 'sim_time > 10 & abs(q_rel) > 0.05'
    fail_mode = HARD
    error_level = ERROR
    message = 'LLB wall in an alpha ramp moves: the grad(alpha) terms are missing or wrong'
    execute_on = 'timestep_end'
  []
[]

[Postprocessors]
  [Iuz]
    type = ElementIntegralVariablePostprocessor
    variable = uz
  []
  [q]
    type = ParsedPostprocessor
    expression = '-0.5*Iuz'
    pp_names = 'Iuz'
  []
  [sim_time]
    type = TimePostprocessor
  []
  [q_rel]
    type = ParsedPostprocessor
    expression = 'q/(0.00396039604*sim_time)'
    pp_names = 'q sim_time'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = bdf2
  dt = 0.05
  end_time = 20
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_rel_step_tol = 1e-12
  line_search = none
  nl_forced_its = 1
[]

[Outputs]
  file_base = out_llb_alpha_ramp_dmi
  csv = true
[]
