[Mesh]
  # SP7 (MnF2 spin flop, Rezende et al., JAP 126, 151101 (2019)): this one-element macrospin is the benchmark mesh. Run to
  # end_time = 1.6 ns for the full flop (m1z = m2z = 0.0958), or set alpha = 0 and sweep H0 = 0-150 kOe for the AFMR branches.
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
    elem_type = HEX8
  []
[]

[GlobalParams]
  mag1_x = mag1_x
  mag1_y = mag1_y
  mag1_z = mag1_z
  mag2_x = mag2_x
  mag2_y = mag2_y
  mag2_z = mag2_z
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'H0 Ms g0 He Ha alpha'
    prop_values = '100.0 1.0 17.59291886 526.0 8.2 0.01'
  []
[]

[Variables]
  [mag1_x]
    initial_condition = 0.00099999983333334
  []
  [mag1_y]
    initial_condition = 0.0
  []
  [mag1_z]
    initial_condition = 0.99999950000004
  []
  [mag2_x]
    initial_condition = -0.00099999983333335
  []
  [mag2_y]
    initial_condition = 0.0
  []
  [mag2_z]
    initial_condition = -0.99999950000004
  []
[]

[Kernels]
  [mag1_x_time]
    type = TimeDerivative
    variable = mag1_x
  []
  [mag1_y_time]
    type = TimeDerivative
    variable = mag1_y
  []
  [mag1_z_time]
    type = TimeDerivative
    variable = mag1_z
  []
  [mag2_x_time]
    type = TimeDerivative
    variable = mag2_x
  []
  [mag2_y_time]
    type = TimeDerivative
    variable = mag2_y
  []
  [mag2_z_time]
    type = TimeDerivative
    variable = mag2_z
  []
  [afmr1_x]
    type = UniaxialAFMSublattice
    variable = mag1_x
    mag_sub = 0
    component = 0
  []
  [afmr1_y]
    type = UniaxialAFMSublattice
    variable = mag1_y
    mag_sub = 0
    component = 1
  []
  [afmr1_z]
    type = UniaxialAFMSublattice
    variable = mag1_z
    mag_sub = 0
    component = 2
  []
  [afmr2_x]
    type = UniaxialAFMSublattice
    variable = mag2_x
    mag_sub = 1
    component = 0
  []
  [afmr2_y]
    type = UniaxialAFMSublattice
    variable = mag2_y
    mag_sub = 1
    component = 1
  []
  [afmr2_z]
    type = UniaxialAFMSublattice
    variable = mag2_z
    mag_sub = 1
    component = 2
  []
[]

[Correctors]
  [unit_m1]
    type = PointwiseRenormalizeVector
    v = 'mag1_x mag1_y mag1_z'
    norm = 1
    execute_on = 'initial timestep_end final'
    force_preaux = true
  []
  [unit_m2]
    type = PointwiseRenormalizeVector
    v = 'mag2_x mag2_y mag2_z'
    norm = 1
    execute_on = 'initial timestep_end final'
    force_preaux = true
  []
[]

[Postprocessors]
  [m1x]
    type = PointValue
    variable = mag1_x
    point = '0 0 0'
    execute_on = 'initial timestep_end'
  []
  [m1y]
    type = PointValue
    variable = mag1_y
    point = '0 0 0'
    execute_on = 'initial timestep_end'
  []
  [m1z]
    type = PointValue
    variable = mag1_z
    point = '0 0 0'
    execute_on = 'initial timestep_end'
  []
  [m2x]
    type = PointValue
    variable = mag2_x
    point = '0 0 0'
    execute_on = 'initial timestep_end'
  []
  [m2y]
    type = PointValue
    variable = mag2_y
    point = '0 0 0'
    execute_on = 'initial timestep_end'
  []
  [m2z]
    type = PointValue
    variable = mag2_z
    point = '0 0 0'
    execute_on = 'initial timestep_end'
  []
  [energy]
    type = ParsedPostprocessor
    expression = '-100.0*(m1z + m2z) + 526.0*(m1x*m2x + m1y*m2y + m1z*m2z) - 0.5*8.2*(m1z*m1z + m2z*m2z)'
    pp_names = 'm1x m1y m1z m2x m2y m2z'
    enable_jit = false
    execute_on = 'initial timestep_end'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -snes_rtol -snes_atol'
    petsc_options_value = 'lu 1e-10 1e-10'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  [TimeIntegrator]
    type = NewmarkBeta
  []
  dt = 5e-5
  num_steps = 20
[]

[Outputs]
  file_base = out_sp7
  csv = true
[]
