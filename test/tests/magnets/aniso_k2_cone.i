# Easy cone from K2: K1 = 1, K2 = -1 relaxes m to cos^2(theta) = K1/(2|K2|) = 0.5 (mz = 0.70711); mesh scaled down from 40 x 20 x 10.
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 4
  ny = 2
  nz = 1
  xmax = 4
  ymax = 2
  zmax = 1
[]

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z
  g0 = 1.0
  mu0 = 1.0
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'alpha Ae Ms K1 K2 nx ny nz'
    prop_values = '1.0 0.1 1.0 1.0 -1.0 0 0 1'
  []
[]

[Variables]
  [mag_x]
    [InitialCondition]
      type = FunctionIC
      function = 'sin(0.1 + 0.1*x)'
    []
  []
  [mag_y]
  []
  [mag_z]
    [InitialCondition]
      type = FunctionIC
      function = 'cos(0.1 + 0.1*x)'
    []
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
  [exchange_x]
    type = MasterExchangeCartLLG
    variable = mag_x
    component = 0
  []
  [exchange_y]
    type = MasterExchangeCartLLG
    variable = mag_y
    component = 1
  []
  [exchange_z]
    type = MasterExchangeCartLLG
    variable = mag_z
    component = 2
  []
  [anisotropy_x]
    type = MasterAnisotropyCartLLG
    variable = mag_x
    component = 0
  []
  [anisotropy_y]
    type = MasterAnisotropyCartLLG
    variable = mag_y
    component = 1
  []
  [anisotropy_z]
    type = MasterAnisotropyCartLLG
    variable = mag_z
    component = 2
  []
[]

[Correctors]
  [unit_magnetization]
    type = PointwiseRenormalizeVector
    v = 'mag_x mag_y mag_z'
    norm = 1
    execute_on = 'timestep_end'
  []
[]

[Postprocessors]
  [mz]
    type = ElementAverageValue
    variable = mag_z
    execute_on = 'initial timestep_end'
  []
  [Fexch]
    type = MasterMagneticExchangeEnergy
    execute_on = 'initial timestep_end'
  []
  [Faniso]
    type = MasterMagneticAnisotropyEnergy
    execute_on = 'initial timestep_end'
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
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  dt = 1.0
  end_time = 30
[]

[Outputs]
  file_base = out_aniso_k2_cone
  csv = true
[]
