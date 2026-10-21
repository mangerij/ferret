[Mesh]
  # muMAG SP5 benchmark mesh: 2.5 nm film cells (40 x 40 in-plane, 4 through the 10 nm), 4 film-size vacuum shells then x1.7
  # grading out to an 800 x 800 x 350 nm box (graded tensor mesh g4a); relax the vortex first (alpha 0.5, u = 0), then
  # restart with the current for 14 ns (NewmarkBeta, IterationAdaptiveDT, dtmax 0.005).
  [box]
    type = CartesianMeshGenerator
    dim = 3
    dx = '5.862556 2.931278 17.587667 2.931278 5.862556'
    dy = '5.862556 2.931278 17.587667 2.931278 5.862556'
    dz = '1.758767 0.879383 1.758767 0.879383 1.758767'
    ix = '1 1 6 1 1'
    iy = '1 1 6 1 1'
    iz = '1 1 2 1 1'
  []
  [center]
    type = TransformGenerator
    input = box
    transform = TRANSLATE
    vector_value = '-17.5876675 -17.5876675 -3.517534'
  []
  [film_block]
    type = SubdomainBoundingBoxGenerator
    input = center
    block_id = 1
    bottom_left = '-8.795 -8.795 -0.8795'
    top_right = '8.795 8.795 0.8795'
  []
  [block_names]
    type = RenameBlockGenerator
    input = film_block
    old_block = '0 1'
    new_block = 'vacuum film'
  []
[]

[GlobalParams]
  mag_x = mag_x
  mag_y = mag_y
  mag_z = mag_z
  potential_H_int = potential_H_int
  Hext_x = Hext_x
  Hext_y = Hext_y
  Hext_z = Hext_z
  mu0 = 1.256637
  g0 = 140.206695
  Hscale = 1.0
[]

[Materials]
  [film_properties]
    type = GenericConstantMaterial
    prop_names = 'alpha Ae Ms permittivity'
    prop_values = '0.1 1.0 1.261566 1.0'
    block = film
  []
  [vacuum_properties]
    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '1.0'
    block = vacuum
  []
[]

[Variables]
  [mag_x]
    block = film
    [InitialCondition]
      type = FunctionIC
      function = '-y/sqrt(x*x + y*y + 1.7587667438*1.7587667438)'
    []
  []
  [mag_y]
    block = film
    [InitialCondition]
      type = FunctionIC
      function = 'x/sqrt(x*x + y*y + 1.7587667438*1.7587667438)'
    []
  []
  [mag_z]
    block = film
    [InitialCondition]
      type = FunctionIC
      function = '1.7587667438/sqrt(x*x + y*y + 1.7587667438*1.7587667438)'
    []
  []
  [potential_H_int]
    block = 'film vacuum'
  []
[]

[AuxVariables]
  [mag_s]
    block = film
  []
  [Hext_x]
    block = film
  []
  [Hext_y]
    block = film
  []
  [Hext_z]
    block = film
  []
[]

[AuxKernels]
  [mag_norm]
    type = VectorMagnitudeAux
    variable = mag_s
    x = mag_x
    y = mag_y
    z = mag_z
    execute_on = 'initial timestep_end'
    block = film
  []
[]

[Correctors]
  [unit_magnetization]
    type = PointwiseRenormalizeVector
    v = 'mag_x mag_y mag_z'
    norm = 1
    execute_on = 'initial timestep_end final'
    force_preaux = true
  []
[]

[Kernels]
  [time_x]
    type = TimeDerivative
    variable = mag_x
    block = film
  []
  [exchange_x]
    type = MasterExchangeCartLLG
    variable = mag_x
    component = 0
    block = film
  []
  [demag_x]
    type = MasterInteractionCartLLGHConst
    variable = mag_x
    component = 0
    block = film
  []
  [stt_x]
    type = MasterZLSpinTransferTorqueCartLLG
    variable = mag_x
    component = 0
    current_direction = 0
    u = -12.6930
    xi = 0.05
    block = film
  []
  [time_y]
    type = TimeDerivative
    variable = mag_y
    block = film
  []
  [exchange_y]
    type = MasterExchangeCartLLG
    variable = mag_y
    component = 1
    block = film
  []
  [demag_y]
    type = MasterInteractionCartLLGHConst
    variable = mag_y
    component = 1
    block = film
  []
  [stt_y]
    type = MasterZLSpinTransferTorqueCartLLG
    variable = mag_y
    component = 1
    current_direction = 0
    u = -12.6930
    xi = 0.05
    block = film
  []
  [time_z]
    type = TimeDerivative
    variable = mag_z
    block = film
  []
  [exchange_z]
    type = MasterExchangeCartLLG
    variable = mag_z
    component = 2
    block = film
  []
  [demag_z]
    type = MasterInteractionCartLLGHConst
    variable = mag_z
    component = 2
    block = film
  []
  [stt_z]
    type = MasterZLSpinTransferTorqueCartLLG
    variable = mag_z
    component = 2
    current_direction = 0
    u = -12.6930
    xi = 0.05
    block = film
  []
  [magnetostatic_laplace]
    type = Electrostatics
    variable = potential_H_int
    block = 'film vacuum'
  []
  [magnetostatic_source]
    type = MagHStrongCart
    variable = potential_H_int
    block = film
  []
[]

[BCs]
  [outer_vacuum]
    type = DirichletBC
    variable = potential_H_int
    boundary = '0 1 2 3 4 5'
    value = 0.0
  []
[]

[Postprocessors]
  [mx]
    type = ElementAverageValue
    variable = mag_x
    block = film
    execute_on = 'initial timestep_end'
  []
  [my]
    type = ElementAverageValue
    variable = mag_y
    block = film
    execute_on = 'initial timestep_end'
  []
  [mz]
    type = ElementAverageValue
    variable = mag_z
    block = film
    execute_on = 'initial timestep_end'
  []
  [mag_norm_min]
    type = NodalExtremeValue
    variable = mag_s
    value_type = min
    block = film
    execute_on = 'initial timestep_end'
  []
  [mag_norm_max]
    type = NodalExtremeValue
    variable = mag_s
    value_type = max
    block = film
    execute_on = 'initial timestep_end'
  []
  [Fexch]
    type = MasterMagneticExchangeEnergy
    energy_scale = 1.0
    block = film
    execute_on = 'initial timestep_end'
  []
  [Fdemag]
    type = MagnetostaticEnergyCart
    energy_scale = 1.0
    block = film
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
  [TimeIntegrator]
    type = NewmarkBeta
  []
  petsc_options_iname = '-pc_type -snes_atol -snes_rtol'
  petsc_options_value = 'lu       1e-8       1e-8'
  dtmin = 1e-6
  dtmax = 5e-3
  num_steps = 10
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 15
    iteration_window = 2
    linear_iteration_ratio = 1000
    dt = 1e-4
    growth_factor = 1.1
    cutback_factor = 0.7
  []
[]

[Outputs]
  file_base = out_sp5
  exodus = true
  csv = true
[]
