um = 0.01
dt_mech = 1.0
dt_polar = 0.1
T = 100.0
alpha1 = ${fparse 0.0404439830432*(1.0/tanh(54.0/T) - 1.0/tanh(54.0/30.0))}
alpha11 = 1.69629782974
alpha12 = 4.44268479218
beta1 = ${fparse 0.000132*(1.0/tanh(145.0/T) - 1.0/tanh(145.0/105.0))}
beta11 = 1.69e-07
beta12 = 3.88e-07
t1111 = 0.000299792012938
t1122 = -2.28639273072e-05
t1212 = -0.000629824117823
C11 = 336
C12 = 107
C44 = 127
Q11 = 0.0457466385977
Q12 = -0.0134813276811
Q44 = 0.00957174265355
R11 = 8.7e-06
R12 = -7.8e-06
R44 = -9.2e-06
p0x = 0.15
p0y = 0
p0z = 0
a0x = 0
a0y = 7
a0z = 0
ts_A = 1e-4
end_time = 3000.0

args = 'um=${um};C11=${C11};C12=${C12};C44=${C44};Q11=${Q11};Q12=${Q12};Q44=${Q44};R11=${R11};R12=${R12};R44=${R44}'
pargs = 'alpha1=${alpha1};alpha11=${alpha11};alpha12=${alpha12};beta1=${beta1};beta11=${beta11};beta12=${beta12};t1111=${t1111};t1122=${t1122};t1212=${t1212};p0x=${p0x};p0y=${p0y};p0z=${p0z};a0x=${a0x};a0y=${a0y};a0z=${a0z};dt_polar=${dt_polar};ts_A=${ts_A}'

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
  []
[]

[Problem]
  solve = false
[]

[Variables]
  [dummy] []
[]

[MultiApps]
  [polar]
    type = TransientMultiApp
    input_files = STO_pertsev_polar.i
    cli_args = '${args};${pargs}'
    execute_on = TIMESTEP_END
    execution_order_group = 0
    sub_cycling = true
  []
  [mech]
    type = TransientMultiApp
    input_files = STO_pertsev_mech.i
    cli_args = '${args}'
    execute_on = TIMESTEP_END
    execution_order_group = 1
  []
[]

[Transfers]
  [u_to_polar]
    type = MultiAppCopyTransfer
    from_multi_app = mech
    to_multi_app = polar
    source_variable = 'u_x u_y u_z'
    variable = 'u_x u_y u_z'
    execute_on = TIMESTEP_END
  []
  [OP_to_mech]
    type = MultiAppCopyTransfer
    from_multi_app = polar
    to_multi_app = mech
    source_variable = 'polar_x polar_y polar_z antiphase_A_x antiphase_A_y antiphase_A_z'
    variable = 'polar_x polar_y polar_z antiphase_A_x antiphase_A_y antiphase_A_z'
    execute_on = TIMESTEP_END
    execute_after_from_multiapp = true
  []
  [pp_Px]
    type = MultiAppPostprocessorTransfer
    from_multi_app = polar
    from_postprocessor = Px
    to_postprocessor = Px
    reduction_type = average
    execute_on = TIMESTEP_END
  []
  [pp_Py]
    type = MultiAppPostprocessorTransfer
    from_multi_app = polar
    from_postprocessor = Py
    to_postprocessor = Py
    reduction_type = average
    execute_on = TIMESTEP_END
  []
  [pp_Pz]
    type = MultiAppPostprocessorTransfer
    from_multi_app = polar
    from_postprocessor = Pz
    to_postprocessor = Pz
    reduction_type = average
    execute_on = TIMESTEP_END
  []
  [pp_Ax]
    type = MultiAppPostprocessorTransfer
    from_multi_app = polar
    from_postprocessor = Ax
    to_postprocessor = Ax
    reduction_type = average
    execute_on = TIMESTEP_END
  []
  [pp_Ay]
    type = MultiAppPostprocessorTransfer
    from_multi_app = polar
    from_postprocessor = Ay
    to_postprocessor = Ay
    reduction_type = average
    execute_on = TIMESTEP_END
  []
  [pp_Az]
    type = MultiAppPostprocessorTransfer
    from_multi_app = polar
    from_postprocessor = Az
    to_postprocessor = Az
    reduction_type = average
    execute_on = TIMESTEP_END
  []
  [pp_drift]
    type = MultiAppPostprocessorTransfer
    from_multi_app = polar
    from_postprocessor = drift
    to_postprocessor = drift
    reduction_type = maximum
    execute_on = TIMESTEP_END
  []
  [pp_Fbulk]
    type = MultiAppPostprocessorTransfer
    from_multi_app = polar
    from_postprocessor = Fbulk
    to_postprocessor = Fbulk
    reduction_type = average
    execute_on = TIMESTEP_END
  []
  [pp_Froto]
    type = MultiAppPostprocessorTransfer
    from_multi_app = polar
    from_postprocessor = Froto
    to_postprocessor = Froto
    reduction_type = average
    execute_on = TIMESTEP_END
  []
  [pp_Fcouple]
    type = MultiAppPostprocessorTransfer
    from_multi_app = polar
    from_postprocessor = Fcouple
    to_postprocessor = Fcouple
    reduction_type = average
    execute_on = TIMESTEP_END
  []
  [pp_Felastic_true]
    type = MultiAppPostprocessorTransfer
    from_multi_app = mech
    from_postprocessor = Felastic_true
    to_postprocessor = Felastic_true
    reduction_type = average
    execute_on = TIMESTEP_END
  []
  [pp_exx]
    type = MultiAppPostprocessorTransfer
    from_multi_app = mech
    from_postprocessor = exx
    to_postprocessor = exx
    reduction_type = average
    execute_on = TIMESTEP_END
  []
  [pp_ezz]
    type = MultiAppPostprocessorTransfer
    from_multi_app = mech
    from_postprocessor = ezz
    to_postprocessor = ezz
    reduction_type = average
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  [Px]
    type = Receiver
  []
  [Py]
    type = Receiver
  []
  [Pz]
    type = Receiver
  []
  [Ax]
    type = Receiver
  []
  [Ay]
    type = Receiver
  []
  [Az]
    type = Receiver
  []
  [drift]
    type = Receiver
  []
  [Fbulk]
    type = Receiver
  []
  [Froto]
    type = Receiver
  []
  [Fcouple]
    type = Receiver
  []
  [Felastic_true]
    type = Receiver
  []
  [exx]
    type = Receiver
  []
  [ezz]
    type = Receiver
  []
  [Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'Fbulk Froto Fcouple'
    pp_coefs = '1 1 1'
  []
  [Ftotal]
    type = LinearCombinationPostprocessor
    pp_names = 'Fbulk Froto Fcouple Felastic_true'
    pp_coefs = '1 1 1 1'
  []
  [dt_pp]
    type = TimestepSize
  []
  [rate]
    type = EnergyRatePostprocessor
    postprocessor = Ftot
    dt = dt_pp
    execute_on = TIMESTEP_END
  []
[]

[Executioner]
  type = Transient
  dt = ${dt_mech}
  end_time = ${end_time}

  num_steps = 4
[]

[Outputs]
  print_linear_residuals = false
  [console]
    type = Console
    print_mesh_changed_info = false
    outlier_variable_norms = false
  []
  [csv]
    type = CSV
    execute_on = 'FINAL'
  []
[]
