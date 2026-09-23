TC   = 25.0
um   = -0.01
l0   = 1.0
L    = 4.0
nx   = 4

alpha1 = ${fparse 3.8e-4*(TC - 479.0)}
G110   = ${fparse l0*l0*abs(alpha1)}

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${nx}
    nz = ${nx}
    xmin = 0.0
    xmax = ${L}
    ymin = 0.0
    ymax = ${L}
    zmin = 0.0
    zmax = ${L}
    elem_type = HEX8
  []
  [pin]
    type = ExtraNodesetGenerator
    input = gen
    new_boundary = 'pin_node'
    coord = '2.0 2.0 2.0'
    use_closest_node = true
  []
[]

[GlobalParams]
  displacements = 'u_x u_y u_z'
  global_strain = global_strain
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = SMALL
    add_variables = true
    eigenstrain_names = 'ferro'
    generate_output = 'strain_xx strain_zz'
  []
[]

[Ferret/CubicParentFEPhaseField]
  electrostatics = false
  elastic = true
  alpha_ijkl = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
  alpha_ijkl_val = '${alpha1} -0.073 0.75 0.26 0.61 -3.7 0.0 0.0 0.0 0.0'
  G_ij = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
  G_ij_val = '${G110} 0.6 0.0 0.3 0.3'
  Q_ij = 'Q11 Q12 Q44'
  Q_ij_val = '0.089 -0.026 0.03375'
  C_ij = 'C11 C12 C44'
  C_ij_val = '174.269 79.029 47.62'
[]

[ICs]
  [Px]
    type = RandomIC
    variable = polar_x
    min = -0.5
    max = 0.5
    seed = 1
  []
  [Py]
    type = RandomIC
    variable = polar_y
    min = -0.5
    max = 0.5
    seed = 2
  []
  [Pz]
    type = RandomIC
    variable = polar_z
    min = -0.5
    max = 0.5
    seed = 3
  []
[]

[AuxVariables]
  [Pmag]
    order = FIRST
    family = LAGRANGE
  []
  [wE]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [Pmag]
    type = ParsedAux
    variable = Pmag
    coupled_variables = 'polar_x polar_y polar_z'
    expression = 'sqrt(polar_x^2 + polar_y^2 + polar_z^2)'
    execute_on = 'initial timestep_end'
  []
  [wE]
    type = WallEnergyDensity
    variable = wE
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
    execute_on = 'initial timestep_end'
  []
[]

[Materials]
  [misfit]
    type = GenericConstantRankTwoTensor
    tensor_name = global_strain
    tensor_values = '${um} 0 0   0 ${um} 0   0 0 0'
  []
[]

[BCs]
  [Periodic]
    [xyz]
      auto_direction = 'x y z'
      variable = 'polar_x polar_y polar_z u_x u_y u_z'
    []
  []
  [pin_ux]
    type = DirichletBC
    variable = u_x
    boundary = 'pin_node'
    value = 0
  []
  [pin_uy]
    type = DirichletBC
    variable = u_y
    boundary = 'pin_node'
    value = 0
  []
  [pin_uz]
    type = DirichletBC
    variable = u_z
    boundary = 'pin_node'
    value = 0
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -ksp_type -ksp_rtol -ksp_gmres_restart'
    petsc_options_value = 'bjacobi  ilu          gmres     1e-6      100'
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  l_max_its = 200
  end_time = 5.0
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5
    growth_factor = 1.3
    cutback_factor = 0.8
    optimal_iterations = 8
    linear_iteration_ratio = 1000
  []
  dtmax = 5.0

  num_steps = 3
[]

[Outputs]
  print_linear_residuals = false
  [exo]
    type = Exodus
  []
[]
