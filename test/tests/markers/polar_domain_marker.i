###############################################################
##
##  Unit test for Ferret's PolarDomainMarker (src/markers/PolarDomainMarker.C)
##
##  PolarDomainMarker is MOOSE's ValueRangeMarker with two deliberate
##  differences, and this test pins down both:
##
##    1. it compares abs(u) instead of u, so the two 180-degree domains
##       of a ferroelectric are treated identically;
##    2. its "inside" state is DO_NOTHING rather than REFINE, so the
##       marker can only ever COARSEN or hold -- REFINE is reachable
##       ONLY through third_state.  (This is why third_state = REFINE
##       below is not optional decoration: without it this marker can
##       never refine anything.)
##
##  No Landau physics is involved.  The driving field is a static,
##  analytic ramp
##
##      P_x(x) = x        on x in [-1, 1]
##
##  chosen over a tanh wall on purpose: a linear ramp puts the marker's
##  decision boundaries at exactly known x, and makes every band several
##  elements wide.  With a tanh profile of realistic width the buffer band
##  is narrower than one element, no quadrature point ever lands inside
##  it, and the marker silently never refines -- a trap worth avoiding in
##  a regression test.
##
##  With lower_bound/upper_bound = 0.4/1.0 and buffer_size = 0.2:
##
##      |x| in [0.4, 1.0]  -> DO_NOTHING   (bulk of both domains)
##      |x| in [0.2, 0.4)  -> REFINE       (third_state, the "wall flanks")
##      |x| <  0.2         -> COARSEN      (the "wall core")
##
##  The mesh is split into two blocks at x = 0 so that the marker field
##  can be integrated over each half separately.  Because P_x is ODD in x
##  and the marker takes abs(), the marked pattern must be MIRROR
##  SYMMETRIC: marker_int_left == marker_int_right exactly.  The
##  Terminator in [UserObjects] asserts that, which makes the abs()
##  behaviour a gold-independent assertion -- if abs() ever regressed to
##  a signed comparison (or to an integer overload) the x < 0 half would
##  be marked COARSEN instead of DO_NOTHING and the test would fail with
##  a clear message rather than a mystery exodiff.
##
###############################################################

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 2
    xmin = -1.0
    xmax = 1.0
    ymin = -0.1
    ymax = 0.1
    elem_type = QUAD4
  []

  ##  block 1 = x > 0, block 0 = x < 0, purely so the two halves can be
  ##  integrated independently for the symmetry assertion.
  [right_half]
    type = SubdomainBoundingBoxGenerator
    input = gen
    bottom_left = '0.0 -0.1 0.0'
    top_right = '1.0 0.1 0.0'
    block_id = 1
    location = INSIDE
  []
[]

[Functions]
  [./ramp]
    type = ParsedFunction
    expression = 'x'
  [../]
[]

##  A trivial diffusion problem exists only to give the Transient
##  executioner something to solve; the marker is driven entirely by the
##  AuxVariable below.
[Variables]
  [./dummy]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxKernels]
  [./polar_x]
    type = FunctionAux
    variable = polar_x
    function = ramp

    ##  NOTE: 'initial' alone is not enough.  The marker is evaluated
    ##  inside the adaptivity cycle; if polar_x has not been recomputed
    ##  on the newly adapted mesh it reads as zero, |0| falls outside
    ##  every band, and the marker coarsens everywhere.
    execute_on = 'initial linear timestep_begin timestep_end'
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = dummy
  [../]
  [./td]
    type = TimeDerivative
    variable = dummy
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = dummy
    boundary = left
    value = 0.0
  [../]
  [./right]
    type = DirichletBC
    variable = dummy
    boundary = right
    value = 1.0
  [../]
[]

[Adaptivity]
  marker = polar_domain
  steps = 1
  max_h_level = 3
  [./Markers]
    [./polar_domain]
      type = PolarDomainMarker
      variable = polar_x
      lower_bound = 0.4
      upper_bound = 1.0
      buffer_size = 0.2
      third_state = REFINE
    [../]
  [../]
[]

[Postprocessors]
  [./n_elem]
    type = NumElements
    execute_on = 'initial timestep_end'
  [../]

  ##  Volume integrals of the marker field over each half.  The marker
  ##  values are the MarkerValue enum: COARSEN = 0, DO_NOTHING = 1,
  ##  REFINE = 2 (DONT_MARK = -1).
  [./marker_int_left]
    type = ElementIntegralVariablePostprocessor
    variable = polar_domain
    block = 0
    execute_on = 'initial timestep_end'
  [../]
  [./marker_int_right]
    type = ElementIntegralVariablePostprocessor
    variable = polar_domain
    block = 1
    execute_on = 'initial timestep_end'
  [../]
  [./marker_asymmetry]
    type = LinearCombinationPostprocessor
    pp_names = 'marker_int_left marker_int_right'
    pp_coefs = '1 -1'
    execute_on = 'initial timestep_end'
  [../]
[]

[UserObjects]

  ###############################################
  ##
  ##  Gold-independent assertion: abs() in PolarDomainMarker must make
  ##  the marking symmetric about x = 0 for this odd input field.
  ##  error_level = ERROR is required -- fail_mode = HARD on its own
  ##  stops the run but still exits 0, which the test harness would
  ##  record as a pass.
  ##
  ###############################################

  [./abs_symmetry_assert]
    type = Terminator
    expression = 'abs(marker_asymmetry) > 1.0e-12'
    fail_mode = HARD
    error_level = ERROR
    execute_on = 'TIMESTEP_END'
    message = 'PolarDomainMarker marked the two 180-degree domains differently: abs(u) is not being applied symmetrically.'
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 5
  dt = 0.1
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_abs_tol = 1e-10
[]

[Outputs]
  ##  A single file_base so that a variant test can rename every output
  ##  with one override (Outputs/file_base=...).
  ##  CSV only: adaptivity writes a separate Exodus file per adapted step
  ##  (.e-s002, .e-s003, ...), which is awkward to gold and unnecessary here
  ##  -- every quantity this test asserts is a postprocessor.
  file_base = out_polar_domain_marker
  csv = true
[]
