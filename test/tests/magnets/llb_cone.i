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
  mu0 = 1.2566370614e-6
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'alpha alpha_perp alpha_par Ms K1 K2 nx ny nz me chi_par'
    prop_values = '1.0 0.8571428571 0.2857142857 860000 292180.3674 -426846.8353 0 0 1 0.8273504024 1.799363238e-10'
  []
[]

[Variables]
  [mag_x]
    initial_condition = 0.5318105875
  []
  [mag_y]
    initial_condition = 0
  []
  [mag_z]
    initial_condition = 0.6337871782
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
  [MasterAnisotropyCartLLG_x]
    type = MasterAnisotropyCartLLG
    variable = mag_x
    component = 0
    llb = true
  []
  [MasterAnisotropyCartLLG_y]
    type = MasterAnisotropyCartLLG
    variable = mag_y
    component = 1
    llb = true
  []
  [MasterAnisotropyCartLLG_z]
    type = MasterAnisotropyCartLLG
    variable = mag_z
    component = 2
    llb = true
  []
  [MasterLongitudinalLLB_x]
    type = MasterLongitudinalLLB
    variable = mag_x
    component = 0
    form = llb
    temperature = 300
    Tc = 700
  []
  [MasterLongitudinalLLB_y]
    type = MasterLongitudinalLLB
    variable = mag_y
    component = 1
    form = llb
    temperature = 300
    Tc = 700
  []
  [MasterLongitudinalLLB_z]
    type = MasterLongitudinalLLB
    variable = mag_z
    component = 2
    form = llb
    temperature = 300
    Tc = 700
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
  [theta]
    type = ParsedPostprocessor
    expression = 'acos(mz/sqrt(mx*mx + my*my + mz*mz))*180/3.141592653589793'
    pp_names = 'mx my mz'
  []
  [mrel]
    type = ParsedPostprocessor
    expression = 'sqrt(mx*mx + my*my + mz*mz)/0.8273504024'
    pp_names = 'mx my mz'
  []
[]

[UserObjects]
  [assert]
    type = Terminator
    expression = 'sim_time > 1.8 & (abs(theta - 45.005576) > 2e-5 | abs(mrel - 1.00009733) > 1e-7)'
    fail_mode = HARD
    error_level = ERROR
    message = 'LLB easy cone: angle or |m| differs from the free-energy minimiser'
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
  end_time = 2.0
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-4
    growth_factor = 1.5
    optimal_iterations = 8
  []
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_rel_step_tol = 1e-12
  line_search = none
  nl_forced_its = 1
[]

[Outputs]
  file_base = out_llb_cone
  csv = true
[]
