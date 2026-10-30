[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 800
  xmax = 54.45427266
[]

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z
  g0 = 2.211e-4
  mu0 = 1.2566370614e-6
  dmi_type = bulk
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'alpha alpha_perp alpha_par Ms Ae D me chi_par'
    prop_values = '1.0 0.7142857143 0.5714285714 860000 13000000 3000000 0.4680904102 1.87707808e-09'
  []
[]

[Variables]
  [mag_x]
    initial_condition = 0.01
  []
  [mag_y]
    [InitialCondition]
      type = FunctionIC
      function = '0.4680904102*cos(0.1153846154*x)'
    []
  []
  [mag_z]
    [InitialCondition]
      type = FunctionIC
      function = '0.4680904102*sin(0.1153846154*x)'
    []
  []
[]

[BCs]
  [Periodic]
    [all]
      variable = 'mag_x mag_y mag_z'
      auto_direction = x
    []
  []
[]

[AuxVariables]
  [mm]
  []
[]

[AuxKernels]
  [mm]
    type = ParsedAux
    variable = mm
    coupled_variables = 'mag_x mag_y mag_z'
    expression = 'sqrt(mag_x^2 + mag_y^2 + mag_z^2)'
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
    llb = true
  []
  [MasterExchangeCartLLG_y]
    type = MasterExchangeCartLLG
    variable = mag_y
    component = 1
    llb = true
  []
  [MasterExchangeCartLLG_z]
    type = MasterExchangeCartLLG
    variable = mag_z
    component = 2
    llb = true
  []
  [MasterDMICartLLG_x]
    type = MasterDMICartLLG
    variable = mag_x
    component = 0
    llb = true
  []
  [MasterDMICartLLG_y]
    type = MasterDMICartLLG
    variable = mag_y
    component = 1
    llb = true
  []
  [MasterDMICartLLG_z]
    type = MasterDMICartLLG
    variable = mag_z
    component = 2
    llb = true
  []
  [MasterLongitudinalLLB_x]
    type = MasterLongitudinalLLB
    variable = mag_x
    component = 0
    form = llb
    temperature = 600
    Tc = 700
  []
  [MasterLongitudinalLLB_y]
    type = MasterLongitudinalLLB
    variable = mag_y
    component = 1
    form = llb
    temperature = 600
    Tc = 700
  []
  [MasterLongitudinalLLB_z]
    type = MasterLongitudinalLLB
    variable = mag_z
    component = 2
    form = llb
    temperature = 600
    Tc = 700
  []
[]

[Postprocessors]
  [sim_time]
    type = TimePostprocessor
  []
  [m_avg]
    type = ElementAverageValue
    variable = mm
  []
  [E_ex]
    type = MasterMagneticExchangeEnergy
  []
  [E_dmi]
    type = MasterMagneticDMIEnergy
  []
  [shift_rel]
    type = ParsedPostprocessor
    expression = '(m_avg/0.4680904102 - 1)/0.0006010528864'
    pp_names = 'm_avg'
  []
  [dens_rel]
    type = ParsedPostprocessor
    expression = '(E_ex + E_dmi)/54.45427266/(-173076.9231*m_avg*m_avg)'
    pp_names = 'E_ex E_dmi m_avg'
  []
[]

[UserObjects]
  [assert]
    type = Terminator
    expression = 'sim_time > 1.0 & (abs(shift_rel - 1) > 0.02 | abs(dens_rel - 1) > 1e-4)'
    fail_mode = HARD
    error_level = ERROR
    message = 'LLB bulk-DMI helix: |m| shift or energy density differs from the analytic LLB result'
    execute_on = 'timestep_end'
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
  end_time = 1.2
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-5
    growth_factor = 1.3
    optimal_iterations = 8
  []
  dtmax = 0.05
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_rel_step_tol = 1e-12
  line_search = none
  nl_forced_its = 1
[]

[Outputs]
  file_base = out_llb_helix
  csv = true
[]
