[Mesh]
  # SP6 2D benchmark: num_sectors = 80, rings = 50 (1 nm rings), relaxed (IterationAdaptiveDT, dtmax 0.01) until rms dm/dt < 1e-6 /ns.
  [disk]
    type = ConcentricCircleMeshGenerator
    num_sectors = 6
    radii = '50'
    rings = '10'
    has_outer_square = false
    preserve_volumes = false
  []
[]

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z
  g0 = 0.0002211
  mu0 = 1.256637061e-06
  Hscale = 1
  dmi_type = interfacial
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'alpha Ms Ae D K1 nx ny nz'
    prop_values = '1 860000 13000000 -3000000 400000 0 0 1'
  []
[]

[Variables]
  [mag_x]
    [InitialCondition]
      type = FunctionIC
      function = '(x/sqrt(x*x+y*y+1e-24)*(1-abs(2*exp(-(x*x+y*y)/625)-1)))/sqrt((1-abs(2*exp(-(x*x+y*y)/625)-1))^2*(x*x+y*y)/(x*x+y*y+1e-24)+(2*exp(-(x*x+y*y)/625)-1)^2)'
    []
  []
  [mag_y]
    [InitialCondition]
      type = FunctionIC
      function = '(y/sqrt(x*x+y*y+1e-24)*(1-abs(2*exp(-(x*x+y*y)/625)-1)))/sqrt((1-abs(2*exp(-(x*x+y*y)/625)-1))^2*(x*x+y*y)/(x*x+y*y+1e-24)+(2*exp(-(x*x+y*y)/625)-1)^2)'
    []
  []
  [mag_z]
    [InitialCondition]
      type = FunctionIC
      function = '(2*exp(-(x*x+y*y)/625)-1)/sqrt((1-abs(2*exp(-(x*x+y*y)/625)-1))^2*(x*x+y*y)/(x*x+y*y+1e-24)+(2*exp(-(x*x+y*y)/625)-1)^2)'
    []
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
  [dmi_x]
    type = MasterDMICartLLG
    variable = mag_x
    component = 0
  []
  [dmi_y]
    type = MasterDMICartLLG
    variable = mag_y
    component = 1
  []
  [dmi_z]
    type = MasterDMICartLLG
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

[Postprocessors]
  [mx]
    type = ElementAverageValue
    variable = mag_x
  []
  [mz]
    type = ElementAverageValue
    variable = mag_z
  []
  [E_exchange]
    type = MasterMagneticExchangeEnergy
    energy_scale = 2e-27
  []
  [E_dmi]
    type = MasterMagneticDMIEnergy
    energy_scale = 2e-27
  []
  [E_anisotropy]
    type = MasterMagneticAnisotropyEnergy
    energy_scale = 2e-27
  []
  [E_total]
    type = LinearCombinationPostprocessor
    pp_names = 'E_exchange E_dmi E_anisotropy'
    pp_coefs = '1 1 1'
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
  dt = 0.005
  num_steps = 10
[]

[Outputs]
  file_base = out_sp6_2d
  exodus = true
  csv = true
[]
