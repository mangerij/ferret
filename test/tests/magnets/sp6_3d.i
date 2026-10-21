[Mesh]
  # SP6 3D benchmark: num_sectors = 72, rings = 46, num_layers = 21 (~1-2 nm), relaxed (IterationAdaptiveDT, dtmax 0.003) until rms dm/dt < 1e-6 /ns.
  [disk]
    type = ConcentricCircleMeshGenerator
    num_sectors = 4
    radii = '91.5'
    rings = '8'
    has_outer_square = false
    preserve_volumes = false
  []
  [cylinder]
    type = AdvancedExtruderGenerator
    input = disk
    heights = '21'
    num_layers = '3'
    direction = '0 0 1'
  []
  [center]
    type = TransformGenerator
    input = cylinder
    transform = TRANSLATE
    vector_value = '0 0 -10.5'
  []
[]

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z
  g0 = 0.0002211
  mu0 = 1.256637061e-06
  Hscale = 1
  dmi_type = bulk
  potential_H_int = potential_H_int
  Hext_x = Hext_x
  Hext_y = Hext_y
  Hext_z = Hext_z
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'alpha Ms Ae D'
    prop_values = '1 384000 8780000 1580000'
  []
[]

[Variables]
  [mag_x]
    [InitialCondition]
      type = FunctionIC
      function = '(-y/sqrt(x*x+y*y+1e-24)*(1-abs((1-2*exp(-(x*x+y*y)/1296)))))/sqrt((1-abs((1-2*exp(-(x*x+y*y)/1296))))^2*(x*x+y*y)/(x*x+y*y+1e-24)+(1-2*exp(-(x*x+y*y)/1296))^2)'
    []
  []
  [mag_y]
    [InitialCondition]
      type = FunctionIC
      function = '(x/sqrt(x*x+y*y+1e-24)*(1-abs((1-2*exp(-(x*x+y*y)/1296)))))/sqrt((1-abs((1-2*exp(-(x*x+y*y)/1296))))^2*(x*x+y*y)/(x*x+y*y+1e-24)+(1-2*exp(-(x*x+y*y)/1296))^2)'
    []
  []
  [mag_z]
    [InitialCondition]
      type = FunctionIC
      function = '(1-2*exp(-(x*x+y*y)/1296))/sqrt((1-abs((1-2*exp(-(x*x+y*y)/1296))))^2*(x*x+y*y)/(x*x+y*y+1e-24)+(1-2*exp(-(x*x+y*y)/1296))^2)'
    []
  []
[]

[AuxVariables]
  [potential_H_int]
  []
  [Hext_x]
  []
  [Hext_y]
  []
  [Hext_z]
    initial_condition = 318309.8862
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
  [zeeman_x]
    type = MasterInteractionCartLLGHConst
    variable = mag_x
    component = 0
  []
  [zeeman_y]
    type = MasterInteractionCartLLGHConst
    variable = mag_y
    component = 1
  []
  [zeeman_z]
    type = MasterInteractionCartLLGHConst
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
    energy_scale = 1e-27
  []
  [E_dmi]
    type = MasterMagneticDMIEnergy
    energy_scale = 1e-27
  []
  [E_zeeman]
    type = MasterMagneticZeemanEnergyCart
    energy_scale = 1e-27
  []
  [E_total]
    type = LinearCombinationPostprocessor
    pp_names = 'E_exchange E_dmi E_zeeman'
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
  num_steps = 5
[]

[Outputs]
  file_base = out_sp6_3d
  exodus = true
  csv = true
[]
