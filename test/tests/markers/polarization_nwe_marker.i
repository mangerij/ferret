###############################################################
##
##  Unit test for Ferret's PolarizationNWEMarker
##  (src/markers/PolarizationNWEMarker.C)
##
##  !!  THIS INPUT CURRENTLY CANNOT RUN.  PolarizationNWEMarker fails  !!
##  !!  during construction for EVERY input, with                      !!
##  !!                                                                 !!
##  !!    ERROR: no MaterialPropertyName parameter named               !!
##  !!           "prop_getter_suffix" found.                           !!
##  !!                                                                 !!
##  !!  Cause: the class derives from QuadraturePointMarker, whose     !!
##  !!  constructor builds a MaterialPropertyInterface and therefore   !!
##  !!  reads 'prop_getter_suffix'.  But its validParams() starts from !!
##  !!  Marker::validParams() instead of QuadraturePointMarker::       !!
##  !!  validParams(), so MaterialPropertyInterface::validParams() is  !!
##  !!  never folded in and that parameter does not exist.             !!
##  !!  (PolarDomainMarker, the working sibling in this directory,     !!
##  !!  correctly starts from QuadraturePointMarker::validParams().)   !!
##  !!                                                                 !!
##  !!  A second, latent defect is waiting behind the first: the       !!
##  !!  member 'bool _invert' is declared in the header, read in       !!
##  !!  computeQpMarker(), but never added as a parameter and never    !!
##  !!  initialized -- so once construction is fixed the refine/coarsen!!
##  !!  sense is chosen by an uninitialized byte.                      !!
##  !!                                                                 !!
##  !!  This input is therefore registered in 'tests' as skipped.      !!
##  !!  Both defects are one-line fixes in the marker source, which is !!
##  !!  out of scope for a test-only change.  Un-skip it once they are !!
##  !!  fixed -- the input below is otherwise complete and exercises   !!
##  !!  the marker's postprocessor gate.                               !!
##
##  What it is intended to cover.  PolarizationNWEMarker gates ordinary
##  threshold marking on a GLOBAL scalar: with
##
##      ten_per = 0.20 * Bulk_Polar          (despite the variable name)
##
##  it marks REFINE everywhere unless ExtremeValue > ten_per, and only
##  above that gate does it compare the local value to refine/coarsen.
##  Note also that AMRoff = true returns REFINE everywhere, which is the
##  OPPOSITE of what the parameter's own description claims.
##
##  Here maxP = 1.0 and Bulk_Polar = 1.0, so ten_per = 0.2 < maxP and the
##  gate is OPEN, selecting the threshold branch:
##
##      P_x > 0.5   -> REFINE
##      P_x < 0.1   -> COARSEN
##      otherwise   -> third_state (DO_NOTHING)
##
##  Unlike PolarDomainMarker this marker does NOT take abs(), so the
##  marking is deliberately asymmetric about x = 0 -- the whole x < 0
##  half falls below 'coarsen'.  The two half-domain integrals below
##  record that asymmetry rather than asserting it away.
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
  marker = polar_nwe
  steps = 1
  max_h_level = 3
  [./Markers]
    [./polar_nwe]
      type = PolarizationNWEMarker
      variable = polar_x
      refine = 0.5
      coarsen = 0.1
      ExtremeValue = maxP
      Bulk_Polar = 1.0
      third_state = DO_NOTHING
    [../]
  [../]
[]

[Postprocessors]

  ##  The global scalar the marker gates on.  It must be available before
  ##  the marker runs, hence the broad execute_on.
  [./maxP]
    type = NodalExtremeValue
    variable = polar_x
    value_type = max
    execute_on = 'initial linear timestep_begin timestep_end'
  [../]

  [./n_elem]
    type = NumElements
    execute_on = 'initial timestep_end'
  [../]
  [./marker_int_left]
    type = ElementIntegralVariablePostprocessor
    variable = polar_nwe
    block = 0
    execute_on = 'initial timestep_end'
  [../]
  [./marker_int_right]
    type = ElementIntegralVariablePostprocessor
    variable = polar_nwe
    block = 1
    execute_on = 'initial timestep_end'
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
  file_base = out_polarization_nwe_marker
  csv = true
[]
