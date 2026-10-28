# AFM skyrmion eigenmodes (SP8; Kravchuk et al., PRB 99, 184429 (2019)): two sublattices, interfacial DMI, field along z,
# undamped NewmarkBeta ring-down of a deformed Neel skyrmion (mu = 0 + mu = 2 kick), reduced units ell = Hsf = gamma = 1,
# xi = 0.05, delta = 0.9, h = 0.1 (Ae = xi/2, D = 2 delta xi/pi, Ha = xi, He = 1/(2 xi), H0 = h).
# To scale to the benchmark: radii = '9 25', rings = '180 80', num_sectors = 40; start from the state relaxed by the same
# kernels (alpha = 0.02, implicit Euler, PointwiseRenormalizeVector) and kick it; run dt = 0.25 to t = 300 and fit the
# frequencies of int n_z cos(mu phi) dA (LFB mu = 0: 0.1508, mu = +-2: 0.4411 / 0.4301 at h = 0.1).
[Mesh]
  [disk]
    type = ConcentricCircleMeshGenerator
    num_sectors = 4
    radii = '4 8'
    rings = '8 4'
    has_outer_square = false
    preserve_volumes = false
  []
[]

[GlobalParams]
  mag1_x = mag1_x
  mag1_y = mag1_y
  mag1_z = mag1_z
  mag2_x = mag2_x
  mag2_y = mag2_y
  mag2_z = mag2_z
  dmi_type = interfacial
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'alpha Ms g0 Ae D He Ha H0'
    prop_values = '0.0 1.0 1.0 0.025 0.0286478897565 10 0.05 0.1'
  []
[]

[Functions]
  [nx]
    type = ParsedFunction
    expression = 'r := sqrt(x*x + y*y + 1e-24); R := 2.5*(1.01 + 0.02*(x*x - y*y)/(r*r)); x/r/cosh(r - R)'
  []
  [ny]
    type = ParsedFunction
    expression = 'r := sqrt(x*x + y*y + 1e-24); R := 2.5*(1.01 + 0.02*(x*x - y*y)/(r*r)); y/r/cosh(r - R)'
  []
  [nz]
    type = ParsedFunction
    expression = 'r := sqrt(x*x + y*y + 1e-24); R := 2.5*(1.01 + 0.02*(x*x - y*y)/(r*r)); tanh(r - R)'
  []
  [mnx]
    type = ParsedFunction
    symbol_names = 'f'
    symbol_values = 'nx'
    expression = '-f'
  []
  [mny]
    type = ParsedFunction
    symbol_names = 'f'
    symbol_values = 'ny'
    expression = '-f'
  []
  [mnz]
    type = ParsedFunction
    symbol_names = 'f'
    symbol_values = 'nz'
    expression = '-f'
  []
[]

[Variables]
  [mag1_x]
    [InitialCondition]
      type = FunctionIC
      function = nx
    []
  []
  [mag1_y]
    [InitialCondition]
      type = FunctionIC
      function = ny
    []
  []
  [mag1_z]
    [InitialCondition]
      type = FunctionIC
      function = nz
    []
  []
  [mag2_x]
    [InitialCondition]
      type = FunctionIC
      function = mnx
    []
  []
  [mag2_y]
    [InitialCondition]
      type = FunctionIC
      function = mny
    []
  []
  [mag2_z]
    [InitialCondition]
      type = FunctionIC
      function = mnz
    []
  []
[]

[Kernels]
  [t1x]
    type = TimeDerivative
    variable = mag1_x
  []
  [ex1x]
    type = MasterExchangeCartLLG
    variable = mag1_x
    component = 0
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
  []
  [dmi1x]
    type = MasterDMICartLLG
    variable = mag1_x
    component = 0
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
  []
  [afm1x]
    type = UniaxialAFMSublattice
    variable = mag1_x
    mag_sub = 0
    component = 0
  []
  [t1y]
    type = TimeDerivative
    variable = mag1_y
  []
  [ex1y]
    type = MasterExchangeCartLLG
    variable = mag1_y
    component = 1
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
  []
  [dmi1y]
    type = MasterDMICartLLG
    variable = mag1_y
    component = 1
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
  []
  [afm1y]
    type = UniaxialAFMSublattice
    variable = mag1_y
    mag_sub = 0
    component = 1
  []
  [t1z]
    type = TimeDerivative
    variable = mag1_z
  []
  [ex1z]
    type = MasterExchangeCartLLG
    variable = mag1_z
    component = 2
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
  []
  [dmi1z]
    type = MasterDMICartLLG
    variable = mag1_z
    component = 2
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
  []
  [afm1z]
    type = UniaxialAFMSublattice
    variable = mag1_z
    mag_sub = 0
    component = 2
  []
  [t2x]
    type = TimeDerivative
    variable = mag2_x
  []
  [ex2x]
    type = MasterExchangeCartLLG
    variable = mag2_x
    component = 0
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
  []
  [dmi2x]
    type = MasterDMICartLLG
    variable = mag2_x
    component = 0
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
  []
  [afm2x]
    type = UniaxialAFMSublattice
    variable = mag2_x
    mag_sub = 1
    component = 0
  []
  [t2y]
    type = TimeDerivative
    variable = mag2_y
  []
  [ex2y]
    type = MasterExchangeCartLLG
    variable = mag2_y
    component = 1
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
  []
  [dmi2y]
    type = MasterDMICartLLG
    variable = mag2_y
    component = 1
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
  []
  [afm2y]
    type = UniaxialAFMSublattice
    variable = mag2_y
    mag_sub = 1
    component = 1
  []
  [t2z]
    type = TimeDerivative
    variable = mag2_z
  []
  [ex2z]
    type = MasterExchangeCartLLG
    variable = mag2_z
    component = 2
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
  []
  [dmi2z]
    type = MasterDMICartLLG
    variable = mag2_z
    component = 2
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
  []
  [afm2z]
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
    execute_on = 'initial timestep_end'
    force_preaux = true
  []
  [unit_m2]
    type = PointwiseRenormalizeVector
    v = 'mag2_x mag2_y mag2_z'
    norm = 1
    execute_on = 'initial timestep_end'
    force_preaux = true
  []
[]

[AuxVariables]
  [n_z]
  []
  [nzc2]
  []
  [m_z]
  []
[]

[AuxKernels]
  [n_z]
    type = ParsedAux
    enable_jit = false
    variable = n_z
    coupled_variables = 'mag1_z mag2_z'
    expression = '(mag1_z - mag2_z)/2'
    execute_on = 'initial timestep_end'
  []
  [nzc2]
    type = ParsedAux
    enable_jit = false
    variable = nzc2
    coupled_variables = 'mag1_z mag2_z'
    use_xyzt = true
    expression = '(mag1_z - mag2_z)/2*(x*x - y*y)/(x*x + y*y + 1e-24)'
    execute_on = 'initial timestep_end'
  []
  [m_z]
    type = ParsedAux
    enable_jit = false
    variable = m_z
    coupled_variables = 'mag1_z mag2_z'
    expression = '(mag1_z + mag2_z)/2'
    execute_on = 'initial timestep_end'
  []
[]

[Postprocessors]
  [nz_int]
    type = ElementIntegralVariablePostprocessor
    variable = n_z
    execute_on = 'initial timestep_end'
  []
  [nzc2_int]
    type = ElementIntegralVariablePostprocessor
    variable = nzc2
    execute_on = 'initial timestep_end'
  []
  [mz_int]
    type = ElementIntegralVariablePostprocessor
    variable = m_z
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
  [TimeIntegrator]
    type = NewmarkBeta
  []
  dt = 0.25
  num_steps = 5
[]

[Outputs]
  file_base = out_afm_skyr_modes
  csv = true
[]
