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
  potential_H_int = zero
  Hext_x = zero
  Hext_y = hy
  Hext_z = hz
  g0 = 1
  mu0 = 1
  Hscale = 1
  llb = true
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'alpha alpha_perp alpha_par Ms K1 nx ny nz Kh nhx nhy nhz me chi_par'
    prop_values = '0.1 0.08333333333 0.03333333333 1 1 1 0 0 0.5 0 0 1 0.8 0.05'
  []
[]

[Variables]
  [mag_x]
    initial_condition = 0.8
  []
  [mag_y]
    initial_condition = 0
  []
  [mag_z]
    initial_condition = 0
  []
[]

[AuxVariables]
  [zero]
  []
  [hy]
    initial_condition = 0.1
  []
  [hz]
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
    expression = 'sim_time > 90 & (abs(my/0.05 - 1) > 1e-10 | abs(mz/0.03333333333 - 1) > 1e-10 | abs(sqrt(mx*mx + my*my + mz*mz)/0.8 - 1) > 1e-10)'
    fail_mode = HARD
    error_level = ERROR
    message = 'biaxial LLB equilibrium: m_y, m_z or |m| off the exact values'
    execute_on = 'timestep_end'
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
  dt = 1
  end_time = 120
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_rel_step_tol = 1e-12
  line_search = none
  nl_forced_its = 1
[]

[Outputs]
  file_base = out_llb_biaxial_chi
  csv = true
[]
