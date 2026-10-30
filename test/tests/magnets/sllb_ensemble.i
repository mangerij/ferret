[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 10
  ny = 10
  nz = 2
  xmax = 53.13292846
  ymax = 53.13292846
  zmax = 10.62658569
[]

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z
  g0 = 2.211e-4
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'alpha alpha_perp alpha_par Ms me chi_par sigma_perp'
    prop_values = '1.0 0.6868686869 0.6262626263 1000000 0.3123200632 1.029614007e-08 125820.8133'
  []
[]

[Variables]
  [mag_x]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  []
  [mag_y]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  []
  [mag_z]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.3123200632
  []
[]

[AuxVariables]
  [zero]
  []
  [one]
    initial_condition = 1
  []
  [mm]
    order = CONSTANT
    family = MONOMIAL
  []
  [mm2]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [mm]
    type = ParsedAux
    variable = mm
    coupled_variables = 'mag_x mag_y mag_z'
    expression = 'sqrt(mag_x^2 + mag_y^2 + mag_z^2)'
  []
  [mm2]
    type = ParsedAux
    variable = mm2
    coupled_variables = 'mag_x mag_y mag_z'
    expression = 'mag_x^2 + mag_y^2 + mag_z^2'
  []
[]

[UserObjects]
  [hx_noise]
    type = LocalConservedNormalNoise
    conserved = false
    white_noise = true
    seed = 1
  []
  [hy_noise]
    type = LocalConservedNormalNoise
    conserved = false
    white_noise = true
    seed = 10000001
  []
  [hz_noise]
    type = LocalConservedNormalNoise
    conserved = false
    white_noise = true
    seed = 20000001
  []
  [px_noise]
    type = LocalConservedNormalNoise
    conserved = false
    white_noise = true
    seed = 30000001
  []
  [py_noise]
    type = LocalConservedNormalNoise
    conserved = false
    white_noise = true
    seed = 40000001
  []
  [pz_noise]
    type = LocalConservedNormalNoise
    conserved = false
    white_noise = true
    seed = 50000001
  []
  [assert]
    type = Terminator
    expression = 'sim_time > 0.002230747728 & (abs(m_avg - 0.31308) > 0.008 | abs(m_std/0.02159 - 1) > 0.2)'
    fail_mode = HARD
    error_level = ERROR
    message = 'sLLB ensemble: <|m|> or its spread is outside the Boltzmann window'
    execute_on = 'timestep_end'
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
  [MasterLongitudinalLLB_x]
    type = MasterLongitudinalLLB
    variable = mag_x
    component = 0
    form = llb
    temperature = 620
    Tc = 660
  []
  [MasterLongitudinalLLB_y]
    type = MasterLongitudinalLLB
    variable = mag_y
    component = 1
    form = llb
    temperature = 620
    Tc = 660
  []
  [MasterLongitudinalLLB_z]
    type = MasterLongitudinalLLB
    variable = mag_z
    component = 2
    form = llb
    temperature = 620
    Tc = 660
  []
  [MasterInteractionCartLLGHConst_x]
    type = MasterInteractionCartLLGHConst
    variable = mag_x
    component = 0
    thermal_field = true
    noise_x = hx_noise
    noise_y = hy_noise
    noise_z = hz_noise
    potential_H_int = zero
    Hext_x = zero
    Hext_y = zero
    Hext_z = zero
  []
  [MasterInteractionCartLLGHConst_y]
    type = MasterInteractionCartLLGHConst
    variable = mag_y
    component = 1
    thermal_field = true
    noise_x = hx_noise
    noise_y = hy_noise
    noise_z = hz_noise
    potential_H_int = zero
    Hext_x = zero
    Hext_y = zero
    Hext_z = zero
  []
  [MasterInteractionCartLLGHConst_z]
    type = MasterInteractionCartLLGHConst
    variable = mag_z
    component = 2
    thermal_field = true
    noise_x = hx_noise
    noise_y = hy_noise
    noise_z = hz_noise
    potential_H_int = zero
    Hext_x = zero
    Hext_y = zero
    Hext_z = zero
  []
  [eta_par_x]
    type = LocalConservedLangevinNoise
    variable = mag_x
    noise = px_noise
    amplitude = 30.71179306
    T_bath = one
  []
  [eta_par_y]
    type = LocalConservedLangevinNoise
    variable = mag_y
    noise = py_noise
    amplitude = 30.71179306
    T_bath = one
  []
  [eta_par_z]
    type = LocalConservedLangevinNoise
    variable = mag_z
    noise = pz_noise
    amplitude = 30.71179306
    T_bath = one
  []
[]

[Postprocessors]
  [m_avg]
    type = ElementAverageValue
    variable = mm
  []
  [m2_avg]
    type = ElementAverageValue
    variable = mm2
  []
  [m_std]
    type = ParsedPostprocessor
    expression = 'sqrt(m2_avg - m_avg*m_avg)'
    pp_names = 'm2_avg m_avg'
  []
  [sim_time]
    type = TimePostprocessor
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
  scheme = implicit-euler
  dt = 1.487165152e-05
  num_steps = 300
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_rel_step_tol = 1e-12
  line_search = bt
  nl_forced_its = 1
[]

[Outputs]
  file_base = out_sllb_ensemble
  csv = true
[]
