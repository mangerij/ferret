TC   = 25.0
um   = -0.01
l0   = 1.0
Lx   = 8.0
tf   = 2.0
nx   = 40
ny   = 10
lam  = 8.0
P_seed = 0.3
Px_seed = 0.01

eps_r = 10.0
eps0  = 0.0088542

alpha1 = ${fparse 3.8e-4*(TC - 479.0)}
G110   = ${fparse l0*l0*abs(alpha1)}

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${nx}
    ny = ${ny}
    xmin = 0.0
    xmax = ${Lx}
    ymin = 0.0
    ymax = ${tf}
    elem_type = QUAD4
  []
  [pin]
    type = ExtraNodesetGenerator
    input = gen
    new_boundary = 'pin_node'
    coord = '8.0 0.0 0.0'
    use_closest_node = true
  []
[]

[GlobalParams]
  displacements = 'u_x u_y'
  global_strain = global_strain
[]

[Functions]
  [ic_Py]
    type = ParsedFunction
    expression = '${P_seed}*cos(2*pi*x/${lam})'
  []
  [ic_Px]
    type = ParsedFunction
    expression = '${Px_seed}*sin(2*pi*x/${lam})'
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = SMALL
    add_variables = true
    eigenstrain_names = 'ferro'
    generate_output = 'strain_xx strain_yy'
  []
[]

[Ferret/CubicParentFEPhaseField]
  electrostatics = true
  elastic = true
  alpha_ijkl = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
  alpha_ijkl_val = '${alpha1} -0.073 0.75 0.26 0.61 -3.7 0.0 0.0 0.0 0.0'
  G_ij = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
  G_ij_val = '${G110} 0.6 0.0 0.3 0.3'
  Q_ij = 'Q11 Q12 Q44'
  Q_ij_val = '0.089 -0.026 0.03375'
  C_ij = 'C11 C12 C44'
  C_ij_val = '174.269 79.029 47.62'
  permittivity_val = '${fparse eps_r*eps0}'
[]

[ICs]
  [Px]
    type = FunctionIC
    variable = polar_x
    function = ic_Px
  []
  [Py]
    type = FunctionIC
    variable = polar_y
    function = ic_Py
  []
  [Pz]
    type = ConstantIC
    variable = polar_z
    value = 0.0
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
    tensor_values = '${um} 0 0   0 0 0   0 0 ${um}'
  []
[]

[BCs]
  [Periodic]
    [x]
      auto_direction = 'x'
      variable = 'polar_x polar_y polar_z potential_E_int u_x u_y'
    []
  []
  [phi_pin]
    type = DirichletBC
    variable = potential_E_int
    boundary = 'pin_node'
    value = 0
  []
  [clamp_x]
    type = DirichletBC
    variable = u_x
    boundary = 'bottom'
    value = 0
  []
  [clamp_y]
    type = DirichletBC
    variable = u_y
    boundary = 'bottom'
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
