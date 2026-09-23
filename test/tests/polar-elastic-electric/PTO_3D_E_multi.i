L    = 4.0
dx   = 1.0

TC   = 25.0
um   = -0.01
l0   = 1.0
eps_r = 10.0

dt_polar = 0.5
dt_mech  = 1.0

args = 'TC=${TC};um=${um};l0=${l0};L=${L};dx=${dx}'
args_polar = '${args};eps_r=${eps_r};dt_polar=${dt_polar}'
args_mech  = '${args}'

[Problem]
  solve = false
[]

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]

[MultiApps]
  [polar]
    type = TransientMultiApp
    input_files = PTO_3D_E_multi_polar.i
    cli_args = '${args_polar}'
    execute_on = TIMESTEP_END
    execution_order_group = 0
    sub_cycling = true
  []
  [mech]
    type = TransientMultiApp
    input_files = PTO_3D_E_multi_mech.i
    cli_args = '${args_mech}'
    execute_on = TIMESTEP_END
    execution_order_group = 1
  []
[]

[Transfers]
  [u_to_polar]
    check_multiapp_execute_on = false
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = mech
    to_multi_app = polar
    source_variable = 'u_x u_y u_z'
    variable = 'u_x u_y u_z'
    execute_on = TIMESTEP_END
  []
  [P_to_mech]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = polar
    to_multi_app = mech
    source_variable = 'polar_x polar_y polar_z'
    variable = 'polar_x polar_y polar_z'
    execute_on = TIMESTEP_END
    execute_after_from_multiapp = true
  []
  [pp_Fbulk]
    type = MultiAppPostprocessorTransfer
    from_multi_app = polar
    from_postprocessor = Fbulk
    to_postprocessor = Fbulk
    reduction_type = average
    execute_on = TIMESTEP_END
  []
  [pp_Fwall]
    type = MultiAppPostprocessorTransfer
    from_multi_app = polar
    from_postprocessor = Fwall
    to_postprocessor = Fwall
    reduction_type = average
    execute_on = TIMESTEP_END
  []
  [pp_Felec]
    type = MultiAppPostprocessorTransfer
    from_multi_app = polar
    from_postprocessor = Felec
    to_postprocessor = Felec
    reduction_type = average
    execute_on = TIMESTEP_END
  []
  [pp_Felastic]
    type = MultiAppPostprocessorTransfer
    from_multi_app = mech
    from_postprocessor = Felastic
    to_postprocessor = Felastic
    reduction_type = average
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  [Fbulk]
    type = Receiver
  []
  [Fwall]
    type = Receiver
  []
  [Felec]
    type = Receiver
  []
  [Felastic]
    type = Receiver
  []
  [Ftotal]
    type = LinearCombinationPostprocessor
    pp_names = 'Fbulk Fwall Felastic Felec'
    pp_coefs = '1 1 1 1'
    execute_on = 'timestep_end'
  []
[]

[Executioner]
  type = Transient
  dt = ${dt_mech}
  num_steps = 3
[]
