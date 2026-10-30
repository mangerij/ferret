[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 300
  xmin = -30
  xmax = 30
[]

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z
  potential_H_int = zero
  Hext_x = hx
  Hext_y = zero
  Hext_z = zero
  g0 = 1
  mu0 = 1
  Hscale = 1
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'alpha Ms Ae K1 nx ny nz Kh nhx nhy nhz'
    prop_values = '0.1 1 1 1 1 0 0 0.5 0 0 1'
  []
[]

[Variables]
  [mag_x]
    [InitialCondition]
      type = FunctionIC
      function = '-tanh(x)'
    []
  []
  [mag_y]
    [InitialCondition]
      type = FunctionIC
      function = '1/cosh(x)'
    []
  []
  [mag_z]
    [InitialCondition]
      type = FunctionIC
      function = '0'
    []
  []
[]

[AuxVariables]
  [zero]
  []
  [hx]
    initial_condition = 0.1
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
  [MasterInteractionCartLLGHConst_x]
    type = MasterInteractionCartLLGHConst
    variable = mag_x
    component = 0
  []
  [MasterInteractionCartLLGHConst_y]
    type = MasterInteractionCartLLGHConst
    variable = mag_y
    component = 1
  []
  [MasterInteractionCartLLGHConst_z]
    type = MasterInteractionCartLLGHConst
    variable = mag_z
    component = 2
  []
[]

[UserObjects]
  [renormalize]
    type = PointwiseRenormalizeVector
    v = 'mag_x mag_y mag_z'
    execute_on = timestep_end
    force_preaux = true
  []
  [assert]
    type = Terminator
    expression = 'sim_time > 20 & abs(w/w_walker - 1) > 0.01'
    fail_mode = HARD
    error_level = ERROR
    message = 'Walker breakdown: the wall precession rate is off (H - H_W sin 2phi)/(1 + alpha^2)'
    execute_on = 'timestep_end'
  []
[]

[Postprocessors]
  [Ix]
    type = ElementIntegralVariablePostprocessor
    variable = mag_x
  []
  [dIx]
    type = ChangeOverTimePostprocessor
    postprocessor = Ix
  []
  [Iy]
    type = ElementIntegralVariablePostprocessor
    variable = mag_y
  []
  [dIy]
    type = ChangeOverTimePostprocessor
    postprocessor = Iy
  []
  [Iz]
    type = ElementIntegralVariablePostprocessor
    variable = mag_z
  []
  [dIz]
    type = ChangeOverTimePostprocessor
    postprocessor = Iz
  []
  [dt]
    type = TimestepSize
  []
  [sim_time]
    type = TimePostprocessor
  []
  [w]
    type = ParsedPostprocessor
    expression = '(Iy*dIz - Iz*dIy)/((Iy*Iy + Iz*Iz)*dt)'
    pp_names = 'Iy Iz dIy dIz dt'
  []
  [w_walker]
    type = ParsedPostprocessor
    expression = '(0.1 - 0.05*sin(2*atan2(Iz, Iy)))/1.01'
    pp_names = 'Iy Iz'
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
  end_time = 40
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_rel_step_tol = 1e-12
  line_search = none
  nl_forced_its = 1
[]

[Outputs]
  file_base = out_walker_llg_above
  csv = true
[]
