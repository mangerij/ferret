[Mesh]
  [line]
    type = CartesianMeshGenerator
    dim = 1
    dx = '5 2 5'
    ix = '5 8 5'
    subdomain_id = '0 1 0'
  []
  [names]
    type = RenameBlockGenerator
    input = line
    old_block = '0 1'
    new_block = 'vacuum film'
  []
[]

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z
  g0 = 2.211e-4
  mu0 = 1.2566370614e-6
  potential_H_int = phi
  Hext_x = zero
  Hext_y = zero
  Hext_z = Hz
[]

[Materials]
  [film]
    type = GenericConstantMaterial
    prop_names = 'alpha alpha_perp alpha_par Ms Ae me chi_par permittivity'
    prop_values = '0.02 0.01428571429 0.01142857143 860000 13000000 0.4680904102 1.87707808e-09 1.0'
    block = film
  []
  [vacuum]
    type = GenericConstantMaterial
    prop_names = 'Ms permittivity'
    prop_values = '0 1.0'
    block = vacuum
  []
[]

[Variables]
  [mag_x]
    block = film
    initial_condition = 0.01633611973
  []
  [mag_y]
    block = film
    initial_condition = 0
  []
  [mag_z]
    block = film
    initial_condition = 0.4678052622
  []
  [phi]
  []
[]

[AuxVariables]
  [zero]
    block = film
  []
  [Hz]
    block = film
    initial_condition = 159154.9431
  []
  [hd_res]
    order = CONSTANT
    family = MONOMIAL
    block = film
  []
  [dphidx]
    order = CONSTANT
    family = MONOMIAL
    block = film
  []
[]

[AuxKernels]
  [dphidx]
    type = VariableGradientComponent
    variable = dphidx
    gradient_variable = phi
    component = x
    block = film
  []
  [hd_res]
    type = ParsedAux
    variable = hd_res
    coupled_variables = 'dphidx mag_x'
    expression = 'abs(dphidx - 860000*mag_x)/860000'
    block = film
  []
[]

[Kernels]
  [time_x]
    type = TimeDerivative
    variable = mag_x
    block = film
  []
  [time_y]
    type = TimeDerivative
    variable = mag_y
    block = film
  []
  [time_z]
    type = TimeDerivative
    variable = mag_z
    block = film
  []
  [MasterExchangeCartLLG_x]
    type = MasterExchangeCartLLG
    variable = mag_x
    component = 0
    block = film
    llb = true
  []
  [MasterExchangeCartLLG_y]
    type = MasterExchangeCartLLG
    variable = mag_y
    component = 1
    block = film
    llb = true
  []
  [MasterExchangeCartLLG_z]
    type = MasterExchangeCartLLG
    variable = mag_z
    component = 2
    block = film
    llb = true
  []
  [MasterInteractionCartLLGHConst_x]
    type = MasterInteractionCartLLGHConst
    variable = mag_x
    component = 0
    block = film
    llb = true
  []
  [MasterInteractionCartLLGHConst_y]
    type = MasterInteractionCartLLGHConst
    variable = mag_y
    component = 1
    block = film
    llb = true
  []
  [MasterInteractionCartLLGHConst_z]
    type = MasterInteractionCartLLGHConst
    variable = mag_z
    component = 2
    block = film
    llb = true
  []
  [MasterLongitudinalLLB_x]
    type = MasterLongitudinalLLB
    variable = mag_x
    component = 0
    block = film
    form = llb
    temperature = 600
    Tc = 700
  []
  [MasterLongitudinalLLB_y]
    type = MasterLongitudinalLLB
    variable = mag_y
    component = 1
    block = film
    form = llb
    temperature = 600
    Tc = 700
  []
  [MasterLongitudinalLLB_z]
    type = MasterLongitudinalLLB
    variable = mag_z
    component = 2
    block = film
    form = llb
    temperature = 600
    Tc = 700
  []
  [phi_laplace]
    type = Electrostatics
    variable = phi
  []
  [phi_source]
    type = MagHStrongCart
    variable = phi
    block = film
  []
[]

[BCs]
  [pin]
    type = DirichletBC
    variable = phi
    boundary = left
    value = 0
  []
[]

[Postprocessors]
  [mx]
    type = ElementAverageValue
    variable = mag_x
    block = film
  []
  [my]
    type = ElementAverageValue
    variable = mag_y
    block = film
  []
  [mz]
    type = ElementAverageValue
    variable = mag_z
    block = film
  []
  [sim_time]
    type = TimePostprocessor
  []
  [hd_err]
    type = ElementExtremeValue
    variable = hd_res
    block = film
  []
  [mrel]
    type = ParsedPostprocessor
    expression = 'sqrt(mx*mx + my*my + mz*mz)/0.4683891565'
    pp_names = 'mx my mz'
  []
[]

[UserObjects]
  [assert]
    type = Terminator
    expression = 'hd_err > 1e-6 | (sim_time > 0.25 & abs(mrel - 1) > 2e-6)'
    fail_mode = HARD
    error_level = ERROR
    message = 'LLB film: demag field is not -Ms0 m_x, or |m| misses the para-process value'
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
  dt = 0.0009506011268
  end_time = 0.3
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8
  nl_rel_step_tol = 1e-12
  line_search = none
  nl_forced_its = 1
[]

[Outputs]
  file_base = out_llb_demag_film
  csv = true
[]
