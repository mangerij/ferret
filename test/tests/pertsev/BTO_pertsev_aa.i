T    = 298.15
TC   = ${fparse T - 273.15}
um = 0.005
Lx = 2.0
Ly = 2.0
tf = 2.0
nx = 2
ny = 2
nz = 8
dt_polar = 0.05
dt_mech  = 0.25
end_time = 400.0
p0x = 0.22
p0y = 0.22
p0z = 0.0
noise = 0.01

C11 = 175.549
C12 = 84.639
umzz = ${fparse -2.0*C12/C11*um}

args = 'TC=${TC};um=${um};umzz=${umzz};Lx=${Lx};Ly=${Ly};tf=${tf};nx=${nx};ny=${ny};nz=${nz}'

[Problem]
  solve = false
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${ny}
    nz = ${nz}
    xmin = 0.0
    xmax = ${Lx}
    ymin = 0.0
    ymax = ${Ly}
    zmin = 0.0
    zmax = ${tf}
    elem_type = HEX8
  []
[]

[MultiApps]
  [polar]
    type = TransientMultiApp
    input_files = BTO_pertsev_polar.i
    cli_args = '${args};dt_polar=${dt_polar};p0x=${p0x};p0y=${p0y};p0z=${p0z};noise=${noise}'
    execute_on = TIMESTEP_END
    execution_order_group = 0
    sub_cycling = true
  []
  [mech]
    type = TransientMultiApp
    input_files = BTO_pertsev_mech.i
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
  [P_to_mech]
    type = MultiAppCopyTransfer
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
  [pp_Pmag_max]
    type = MultiAppPostprocessorTransfer
    from_multi_app = polar
    from_postprocessor = Pmag_max
    to_postprocessor = Pmag_max
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
  [pp_Felastic]
    type = MultiAppPostprocessorTransfer
    from_multi_app = mech
    from_postprocessor = Felastic
    to_postprocessor = Felastic
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
  [pp_exz]
    type = MultiAppPostprocessorTransfer
    from_multi_app = mech
    from_postprocessor = exz
    to_postprocessor = exz
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
  [Px]
    type = Receiver
  []
  [Py]
    type = Receiver
  []
  [Pz]
    type = Receiver
  []
  [Pmag_max]
    type = Receiver
  []
  [Felastic]
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
  [exz]
    type = Receiver
  []
  [Ftotal]
    type = LinearCombinationPostprocessor
    pp_names = 'Fbulk Fwall Felastic'
    pp_coefs = '1 1 1'
    execute_on = 'timestep_end'
  []
  [t]
    type = TimePostprocessor
    execute_on = 'timestep_end'
  []
[]

[Executioner]
  type = Transient
  dt = ${dt_mech}
  end_time = ${end_time}

  num_steps = 4
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = 'FINAL'
  []
[]
