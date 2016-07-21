#try elastic block with MortarPeriodicMesh

[Mesh]
  type = MortarPeriodicMesh
  dim = 3
  nx = 8
  ny = 8
  nz = 8
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  zmin = -0.5
  zmax = 0.5

  periodic_directions = 'x y z'

  [./MortarInterfaces] #define where the lm lives
    [./left_right]
      master = 1
      slave = 3
      subdomain = 10
    [../]
    [./up_down]
      master = 0
      slave = 2
      subdomain = 11
    [../]
    [./front_back]
      master = 4
      slave = 5
      subdomain = 12
    [../]
  [../]
[]

[MeshModifiers]
  [./cnode] #center node
    type = AddExtraNodeset
    coord = '0.0 0.0 0.0'
    new_boundary = 100
  [../]
  [./anode] #angular fixed node to remove the rigid body nullspace mode
    type = AddExtraNodeset
    coord = '0.5 0.5 0.5' 
    new_boundary = 101
  [../]
  #[./dnode]
  #  type = AddExtraNodeset
  #  coord = '0.0 -0.5 0.0'
  #  new_boundary = 102
  #[../]
[]

[GlobalParams]
  derivative_order = 2
  enable_jit = true
  displacements = 'disp_x disp_y disp_z'
 # prefactor = 0.01 #negative = tension, positive = compression
[]

[Variables]
  # Solute concentration variable
  [./c]
    [./InitialCondition]
      type = RandomIC
      min = 0.49
      max = 0.51
    [../]
    block = 0
  [../]
  [./w]
    block = 0
  [../]

  #elastic fields
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]

  ## Lagrange multipliers for gradient component periodicity
  #lm left right interface
  [./lm_left_right_xx]
    order = FIRST
    family = LAGRANGE
    block = 10
  [../]
  [./lm_left_right_xy]
    order = FIRST
    family = LAGRANGE
    block = 10
  [../]
  [./lm_left_right_xz]
    order = FIRST
    family = LAGRANGE
    block = 10
  [../]

  [./lm_left_right_yx]
    order = FIRST
    family = LAGRANGE
    block = 10
  [../]
  [./lm_left_right_yy]
    order = FIRST
    family = LAGRANGE
    block = 10
  [../]
  [./lm_left_right_yz]
    order = FIRST
    family = LAGRANGE
    block = 10
  [../]

  [./lm_left_right_zx]
    order = FIRST
    family = LAGRANGE
    block = 10
  [../]
  [./lm_left_right_zy]
    order = FIRST
    family = LAGRANGE
    block = 10
  [../]
  [./lm_left_right_zz]
    order = FIRST
    family = LAGRANGE
    block = 10
  [../]

  #lm up down interface
  [./lm_up_down_xx]
    order = FIRST
    family = LAGRANGE
    block = 11
  [../]
  [./lm_up_down_xy]
    order = FIRST
    family = LAGRANGE
    block = 11
  [../]
  [./lm_up_down_xz]
    order = FIRST
    family = LAGRANGE
    block = 11
  [../]

  [./lm_up_down_yx]
    order = FIRST
    family = LAGRANGE
    block = 11
  [../]
  [./lm_up_down_yy]
    order = FIRST
    family = LAGRANGE
    block = 11
  [../]
  [./lm_up_down_yz]
    order = FIRST
    family = LAGRANGE
    block = 11
  [../]

  [./lm_up_down_zx]
    order = FIRST
    family = LAGRANGE
    block = 11
  [../]
  [./lm_up_down_zy]
    order = FIRST
    family = LAGRANGE
    block = 11
  [../]
  [./lm_up_down_zz]
    order = FIRST
    family = LAGRANGE
    block = 11
  [../]

  #lm front back interface
  [./lm_front_back_xx]
    order = FIRST
    family = LAGRANGE
    block = 12
  [../]
  [./lm_front_back_xy]
    order = FIRST
    family = LAGRANGE
    block = 12
  [../]
  [./lm_front_back_xz]
    order = FIRST
    family = LAGRANGE
    block = 12
  [../]

  [./lm_front_back_yx]
    order = FIRST
    family = LAGRANGE
    block = 12
  [../]
  [./lm_front_back_yy]
    order = FIRST
    family = LAGRANGE
    block = 12
  [../]
  [./lm_front_back_yz]
    order = FIRST
    family = LAGRANGE
    block = 12
  [../]

  [./lm_front_back_zx]
    order = FIRST
    family = LAGRANGE
    block = 12
  [../]
  [./lm_front_back_zy]
    order = FIRST
    family = LAGRANGE
    block = 12
  [../]
  [./lm_front_back_zz]
    order = FIRST
    family = LAGRANGE
    block = 12
  [../]
[]

[Constraints]
  #left right LM constraints

 [./lr_disp_x_grad_x]
    type = EqualGradientConstraint
    variable = lm_left_right_xx
    interface = left_right
    component = 0
    master_variable = disp_x
  [../]
  [./lr_disp_x_grad_y]
    type = EqualGradientConstraint
    variable = lm_left_right_xy
    interface = left_right
    component = 1
    master_variable = disp_x
  [../]
  [./lr_disp_x_grad_z]
    type = EqualGradientConstraint
    variable = lm_left_right_xz
    interface = left_right
    component = 2
    master_variable = disp_x
  [../]

  [./lr_disp_y_grad_x]
    type = EqualGradientConstraint
    variable = lm_left_right_yx
    interface = left_right
    component = 0
    master_variable = disp_y
  [../]
  [./lr_disp_y_grad_y]
    type = EqualGradientConstraint
    variable = lm_left_right_yy
    interface = left_right
    component = 1
    master_variable = disp_y
  [../]
  [./lr_disp_y_grad_z]
    type = EqualGradientConstraint
    variable = lm_left_right_yz
    interface = left_right
    component = 2
    master_variable = disp_y
  [../]

  [./lr_disp_z_grad_x]
    type = EqualGradientConstraint
    variable = lm_left_right_zx
    interface = left_right
    component = 0
    master_variable = disp_z
  [../]
  [./lr_disp_z_grad_y]
    type = EqualGradientConstraint
    variable = lm_left_right_zy
    interface = left_right
    component = 1
    master_variable = disp_z
  [../]
  [./lr_disp_z_grad_z]
    type = EqualGradientConstraint
    variable = lm_left_right_zz
    interface = left_right
    component = 2
    master_variable = disp_z
  [../]

  #up down LM constraints

 [./ud_disp_x_grad_x]
    type = EqualGradientConstraint
    variable = lm_up_down_xx
    interface = up_down
    component = 0
    master_variable = disp_x
  [../]
  [./ud_disp_x_grad_y]
    type = EqualGradientConstraint
    variable = lm_up_down_xy
    interface = up_down
    component = 1
    master_variable = disp_x
  [../]
  [./ud_disp_x_grad_z]
    type = EqualGradientConstraint
    variable = lm_up_down_xz
    interface = up_down
    component = 2
    master_variable = disp_x
  [../]

  [./ud_disp_y_grad_x]
    type = EqualGradientConstraint
    variable = lm_up_down_yx
    interface = up_down
    component = 0
    master_variable = disp_y
  [../]
  [./ud_disp_y_grad_y]
    type = EqualGradientConstraint
    variable = lm_up_down_yy
    interface = up_down
    component = 1
    master_variable = disp_y
  [../]
  [./ud_disp_y_grad_z]
    type = EqualGradientConstraint
    variable = lm_up_down_yz
    interface = up_down
    component = 2
    master_variable = disp_y
  [../]

  [./ud_disp_z_grad_x]
    type = EqualGradientConstraint
    variable = lm_up_down_zx
    interface = up_down
    component = 0
    master_variable = disp_z
  [../]
  [./ud_disp_z_grad_y]
    type = EqualGradientConstraint
    variable = lm_up_down_zy
    interface = up_down
    component = 1
    master_variable = disp_z
  [../]
  [./ud_disp_z_grad_z]
    type = EqualGradientConstraint
    variable = lm_up_down_zz
    interface = up_down
    component = 2
    master_variable = disp_z
  [../]

  #front back LM constraints

  [./fb_disp_x_grad_x]
    type = EqualGradientConstraint
    variable = lm_front_back_xx
    interface = front_back
    component = 0
    master_variable = disp_x
  [../]
  [./fb_disp_x_grad_y]
    type = EqualGradientConstraint
    variable = lm_front_back_xy
    interface = front_back
    component = 1
    master_variable = disp_x
  [../]
  [./fb_disp_x_grad_z]
    type = EqualGradientConstraint
    variable = lm_front_back_xz
    interface = front_back
    component = 2
    master_variable = disp_x
  [../]

  [./fb_disp_y_grad_x]
    type = EqualGradientConstraint
    variable = lm_front_back_yx
    interface = front_back
    component = 0
    master_variable = disp_y
  [../]
  [./fb_disp_y_grad_y]
    type = EqualGradientConstraint
    variable = lm_front_back_yy
    interface = front_back
    component = 1
    master_variable = disp_y
  [../]
  [./fb_disp_y_grad_z]
    type = EqualGradientConstraint
    variable = lm_front_back_yz
    interface = front_back
    component = 2
    master_variable = disp_y
  [../]

  [./fb_disp_z_grad_x]
    type = EqualGradientConstraint
    variable = lm_front_back_zx
    interface = front_back
    component = 0
    master_variable = disp_z
  [../]
  [./fb_disp_z_grad_y]
    type = EqualGradientConstraint
    variable = lm_front_back_zy
    interface = front_back
    component = 1
    master_variable = disp_z
  [../]
  [./fb_disp_z_grad_z]
    type = EqualGradientConstraint
    variable = lm_front_back_zz
    interface = front_back
    component = 2
    master_variable = disp_z
  [../]

[]

[AuxVariables]
  [./local_energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx_elastic]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./stress_yy_elastic]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./stress_xy_elastic]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./stress_xz_elastic]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./stress_zz_elastic]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./stress_yz_elastic]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./strain_xx_elastic]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./strain_yy_elastic]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./strain_xy_elastic]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./strain_xz_elastic]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./strain_zz_elastic]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./strain_yz_elastic]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
[]


[AuxKernels]
  [./local_free_energy]
    type = TotalFreeEnergy
    block = 0
    execute_on = 'initial LINEAR'
    variable = local_energy
    interfacial_vars = 'c'
    kappa_names = 'kappa_c'
  [../]
  [./matl_e11]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    variable = strain_xx_elastic
    execute_on = 'timestep_end'
    block = 0
  [../]
  [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    variable = strain_xy_elastic
    execute_on = 'timestep_end'
    block = 0
  [../]
  [./matl_e13]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 2
    variable = strain_xz_elastic
    execute_on = 'timestep_end'
    block = 0
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    variable = strain_yy_elastic
    execute_on = 'timestep_end'
    block = 0
  [../]
  [./matl_e23]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 2
    variable = strain_yz_elastic
    execute_on = 'timestep_end'
    block = 0
  [../]
  [./matl_e33]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 2
    index_j = 2
    variable = strain_zz_elastic
    execute_on = 'timestep_end'
    block = 0
  [../]
  [./matl_s11]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = stress_xx_elastic
    execute_on = 'timestep_end'
    block = 0
  [../]
  [./matl_s12]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = stress_xy_elastic
    execute_on = 'timestep_end'
    block = 0
  [../]
  [./matl_s13]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
    variable = stress_xz_elastic
    execute_on = 'timestep_end'
    block = 0
  [../]
  [./matl_s22]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = stress_yy_elastic
    execute_on = 'timestep_end'
    block = 0
  [../]
  [./matl_s23]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    variable = stress_yz_elastic
    execute_on = 'timestep_end'
    block = 0
  [../]
  [./matl_s33]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    variable = stress_zz_elastic
    execute_on = 'timestep_end'
    block = 0
  [../]
[]

[Kernels]
  # Set up stress divergence kernels
  [./TensorMechanics]
  [../]

  # Cahn-Hilliard kernels
  [./c_dot]
    type = CoupledTimeDerivative
    variable = w
    v = c
  [../]
  [./c_res]
    type = SplitCHParsed
    variable = c
    f_name = F
    kappa_name = kappa_c
    w = w
  [../]
  [./w_res]
    type = SplitCHWRes
    variable = w
    mob_name = M
  [../]
[]

[Materials]
  # declare a few constants, such as mobilities (L,M) and interface gradient prefactors (kappa*)
  [./consts]
    type = GenericConstantMaterial
    block = '0 10 11 12'
    prop_names  = 'M   kappa_c dd'
    prop_values = '0.2 0.01   -0.3'
  [../]

  [./shear1]
    type = GenericConstantRankTwoTensor
    block = 0
    tensor_values = '0 0 0 0 0 0.5'
    tensor_name = shear1
  [../]
  [./shear2]
    type = GenericConstantRankTwoTensor
    block = 0
    tensor_values = '0 0 0 0 0 -0.5'
    tensor_name = shear2
  [../]
  [./expand3]
    type = GenericConstantRankTwoTensor
    block = 0
    tensor_values = '1 1 0 0 0 0'
    tensor_name = expand3
  [../]

  [./weight1]
    type = DerivativeParsedMaterial
    block = 0
    function = '0.3*c^2'
    f_name = weight1
    args = c
  [../]
  [./weight2]
    type = DerivativeParsedMaterial
    block = 0
    function = '0.3*(1-c)^2'
    f_name = weight2
    args = c
  [../]
  [./weight3]
    type = DerivativeParsedMaterial
    block = 0
    function = '4*(0.5-c)^2'
    f_name = weight3
    args = c
  [../]

  # matrix phase
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = 0
    C_ijkl = '1 1'
    fill_method = symmetric_isotropic
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = 0
    displacements = 'disp_x disp_y disp_z'
  [../]

  [./eigenstrain]
    type = CompositeEigenstrain
    block = 0
    tensors = 'shear1  shear2  expand3'
    weights = 'weight1 weight2 weight3'
    args = c
  [../]

  [./stress]
    type = ComputeLinearElasticStress
    block = 0
  [../]

  # chemical free energies
  [./chemical_free_energy]
    type = DerivativeParsedMaterial
    block = 0
    f_name = Fc
    function = '4*c^2*(1-c)^2'
    args = 'c'
    #outputs = exodus
    output_properties = Fc
  [../]

  # elastic free energies
  [./elastic_free_energy]
    type = ElasticEnergyMaterial
    f_name = Fe
    block = 0
    args = 'c'
    #outputs = exodus
    output_properties = Fe
  [../]

  # free energy (chemical + elastic)
  [./free_energy]
    type = DerivativeSumMaterial
    block = 0
    f_name = F
    sum_materials = 'Fc Fe'
    args = 'c'
  [../]
[]


[BCs]
  [./Periodic]
    [./up_down]
      primary = top
      secondary = bottom
      translation = '0 -1 0'
      variable = 'c w'
    [../]
    [./left_right]
      primary = left
      secondary = right
      translation = '1 0 0'
      variable = 'c w'
    [../]
    [./front_back]
      primary = front
      secondary = back
      translation = '0 0 -1'
      variable = 'c w'
    [../]
  [../]

  # fix center point location
  [./centerfix_x]
    type = PresetBC
    boundary = 100
    variable = disp_x
    value = 0
  [../]
  [./centerfix_y]
    type = PresetBC
    boundary = 100
    variable = disp_y
    value = 0
  [../]
  [./centerfix_z]
    type = PresetBC
    boundary = 100
    variable = disp_z
    value = 0
  [../]

  # fix side point x coordinate to inhibit rotation
  [./angularfix]
    type = PresetBC
    boundary = 101
    variable = disp_x
    value = 0
  [../]
[]

# We monitor the total free energy and the total solute concentration (should be constant)
[Postprocessors]
  [./total_free_energy]
    type = ElementIntegralVariablePostprocessor
    block = 0
    execute_on = 'initial TIMESTEP_END'
    variable = local_energy
  [../]
  [./total_solute]
    type = ElementIntegralVariablePostprocessor
    block = 0
    execute_on = 'initial TIMESTEP_END'
    variable = c
  [../]
  [./min]
    type = ElementExtremeValue
    block = 0
    execute_on = 'initial TIMESTEP_END'
    value_type = min
    variable = c
  [../]
  [./max]
    type = ElementExtremeValue
    block = 0
    execute_on = 'initial TIMESTEP_END'
    value_type = max
    variable = c
  [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason'
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = 'NEWTON'

  line_search = basic

  # mortar currently does not support MPI parallelization
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = ' lu       NONZERO               1e-10'

  l_max_its = 60
  nl_max_its = 12

  l_tol = 1.0e-4

  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10

  start_time = 0.0
  num_steps = 100

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 0.01
  [../]
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_elastic_mortar
    elemental_as_nodal = true
    interval = 1
  [../]
[]
