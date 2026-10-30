[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmax = 2
  ymax = 2
  zmax = 2
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
    prop_names = 'alpha alpha_par Ms me chi_par'
    prop_values = '0.1 0.06262626263 1e6 0.3123200632 1.029614007e-08'
  []
[]

[Variables]
  [mag_x]
    initial_condition = 0
  []
  [mag_y]
    initial_condition = 0
  []
  [mag_z]
    initial_condition = 0.3126323832
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
[]

[Postprocessors]
  [mx]
    type = ElementAverageValue
    variable = mag_x
  []
  [my]
    type = ElementAverageValue
    variable = mag_y
  []
  [mz]
    type = ElementAverageValue
    variable = mag_z
  []
  [sim_time]
    type = TimePostprocessor
  []
  [ratio]
    type = ParsedPostprocessor
    expression = '(sqrt(mx*mx + my*my + mz*mz) - 0.3123200632)/(0.0003123200632)'
    pp_names = 'mx my mz'
  []
[]

[UserObjects]
  [assert]
    type = Terminator
    expression = 'sim_time > 0.0001502036803 & abs(ratio - exp(-sim_time/0.0007510184017)) > 2e-3'
    fail_mode = HARD
    error_level = ERROR
    message = 'LLB longitudinal relaxation deviates from exp(-t/tau_par) at T = 620 K'
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
  dt = 7.510184017e-06
  end_time = 0.001126527603
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_rel_step_tol = 1e-12
  line_search = none
  nl_forced_its = 1
[]

[Outputs]
  file_base = out_llb_tau_par
  csv = true
[]
