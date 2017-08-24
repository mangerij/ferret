#try elastic block with MortarPeriodicMesh

[Mesh]
  type = MortarPeriodicMesh
  dim = 3
  nx = 9
  ny = 9
  nz = 9
  xmin = -3.0
  xmax = 3.0
  ymin = -3.0
  ymax = 3.0
  zmin = -3.0
  zmax = 3.0

  periodic_directions = 'x y'

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
 #   [./front_back]
 #     master = 4
 #     slave = 5
 #     subdomain = 12
 #   [../]
  [../]
[]

#[MeshModifiers]
#  [./cnode] #center node
#    type = AddExtraNodeset
#    coord = '0.0 0.0 0.0'
#    new_boundary = 100
#  [../]
#  [./anode] #angular fixed node to remove the rigid body nullspace mode
#    type = AddExtraNodeset
#    coord = '-1.0 -1.0 -1.0'
#    new_boundary = 101
#  [../]
#[]

[GlobalParams]
  len_scale = 1.0
  alpha1 = -0.1722883 # (3.766(T-765.1)*10^5) C^{-2} nm^2 (T = 293 K)
  alpha11 = -0.07253
  alpha111 = 0.26
  alpha12 = 0.75
  alpha112 = 0.61
  alpha123 = -3.67
  G110 = 0.173
  G11_G110 = 0.6
  G12_G110 = 0
  G44_G110 = 0.3
  G44P_G110 = 0.3
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  displacements = 'disp_x disp_y disp_z'
  enable_jit = true
 # prefactor = 0.01 #negative = tension, positive = compression
[]

[Variables]
  [./polar_x]
    block = '0'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
  [../]
  [./polar_y]
    block = '0'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
  [../]
  [./polar_z]
    block = '0'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
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

  ##lm front back interface
  #[./lm_front_back_xx]
  #  order = FIRST
  #  family = LAGRANGE
  #  block = 12
  #[../]
  #[./lm_front_back_xy]
  #  order = FIRST
  #  family = LAGRANGE
  #  block = 12
  #[../]
  #[./lm_front_back_xz]
  #  order = FIRST
  #  family = LAGRANGE
  #  block = 12
  #[../]
  #
  #[./lm_front_back_yx]
  #  order = FIRST
  #  family = LAGRANGE
  #  block = 12
  #[../]
  #[./lm_front_back_yy]
  #  order = FIRST
  #  family = LAGRANGE
  #  block = 12
  #[../]
  #[./lm_front_back_yz]
  #  order = FIRST
  #  family = LAGRANGE
  #  block = 12
  #[../]
  #
  #[./lm_front_back_zx]
  #  order = FIRST
  #  family = LAGRANGE
  #  block = 12
  #[../]
  #[./lm_front_back_zy]
  #  order = FIRST
  #  family = LAGRANGE
  #  block = 12
  #[../]
  #[./lm_front_back_zz]
  #  order = FIRST
  #  family = LAGRANGE
  #  block = 12
  #[../]
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

  #[./fb_disp_x_grad_x]
  #  type = EqualGradientConstraint
  #  variable = lm_front_back_xx
  #  interface = front_back
  #  component = 0
  #  master_variable = disp_x
  #[../]
  #[./fb_disp_x_grad_y]
  #  type = EqualGradientConstraint
  #  variable = lm_front_back_xy
  #  interface = front_back
  #  component = 1
  #  master_variable = disp_x
  #[../]
  #[./fb_disp_x_grad_z]
  #  type = EqualGradientConstraint
  #  variable = lm_front_back_xz
  #  interface = front_back
  #  component = 2
  #  master_variable = disp_x
  #[../]
  #
  #[./fb_disp_y_grad_x]
  #  type = EqualGradientConstraint
  #  variable = lm_front_back_yx
  #  interface = front_back
  #  component = 0
  #  master_variable = disp_y
  #[../]
  #[./fb_disp_y_grad_y]
  #  type = EqualGradientConstraint
  #  variable = lm_front_back_yy
  #  interface = front_back
  #  component = 1
  #  master_variable = disp_y
  #[../]
  #[./fb_disp_y_grad_z]
  #  type = EqualGradientConstraint
  #  variable = lm_front_back_yz
  #  interface = front_back
  #  component = 2
  #  master_variable = disp_y
  #[../]
  #
  #[./fb_disp_z_grad_x]
  #  type = EqualGradientConstraint
  #  variable = lm_front_back_zx
  #  interface = front_back
  #  component = 0
  #  master_variable = disp_z
  #[../]
  #[./fb_disp_z_grad_y]
  #  type = EqualGradientConstraint
  #  variable = lm_front_back_zy
  #  interface = front_back
  #  component = 1
  #  master_variable = disp_z
  #[../]
  #[./fb_disp_z_grad_z]
  #  type = EqualGradientConstraint
  #  variable = lm_front_back_zz
  #  interface = front_back
  #  component = 2
  #  master_variable = disp_z
  #[../]

[]

[AuxVariables]
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
  #Bulk energy density
  [./bed_x]
    type = BulkEnergyDerivativeSixth
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = BulkEnergyDerivativeSixth
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = BulkEnergyDerivativeSixth
    variable = polar_z
    component = 2
  [../]
  ##Wall energy penalty
[./walled_x]
    type = WallEnergyDerivative
    variable = polar_x
    component = 0
 [../]
 [./walled_y]
    type = WallEnergyDerivative
    variable = polar_y
    component = 1
  [../]
  [./walled_z]
     type = WallEnergyDerivative
     variable = polar_z
     component = 2
  [../]
##Polarization-strain coupling

  [./ferroelectriccouplingp_xx]
    type = FerroelectricCouplingP
    variable = polar_x
    component = 0
  [../]
  [./ferroelectriccouplingp_yy]
    type = FerroelectricCouplingP
    variable = polar_y
    component = 1
  [../]
  [./ferroelectriccouplingp_zz]
    type = FerroelectricCouplingP
    variable = polar_z
    component = 2
  [../]


  [./ferroelectriccouplingX_xx]
    type = FerroelectricCouplingX
    block = '0'
    variable = disp_x
    component = 0
  [../]
  [./ferroelectriccouplingX_yy]
    type = FerroelectricCouplingX
    block = '0'
    variable = disp_y
    component = 1
  [../]
  [./ferroelectriccouplingX_zz]
    type = FerroelectricCouplingX
    block = '0'
    variable = disp_z
    component = 2
  [../]
  ##Electrostatics
  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
     permittivity = 0.08854187
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     block = '0'
     permittivity = 0.08854187
  [../]


  [./polar_electric_px]
     type = PolarElectricPStrong
     variable = polar_x
     component = 0
  [../]
  [./polar_electric_py]
     type = PolarElectricPStrong
     variable = polar_y
     component = 1
  [../]
  [./polar_electric_pz]
     type = PolarElectricPStrong
     variable = polar_z
     component = 2
  [../]
  ##Time dependence
  [./polar_x_time]
     type = TimeDerivativeScaled
     variable=polar_x
    time_scale = 1.0
  [../]
  [./polar_y_time]
     type = TimeDerivativeScaled
     variable=polar_y
    time_scale = 1.0
  [../]
  [./polar_z_time]
     type = TimeDerivativeScaled
     variable = polar_z
    time_scale = 1.0
  [../]
[]

[Materials]
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    #from MaterialsProject
    C_ijkl = '281 115.74 115.74 281 115.74 281 97.18 97.18 97.18'
    block = '0'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    block = '0'
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
    block = '0'
  [../]

  [./slab_ferroelectric]
    block = '0'
    type = ComputeElectrostrictiveTensor
    Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.03375 0.03375 0.03375'
    C_ijkl = '281 115.74 115.74 281 115.74 281 97.18 97.18 97.18'
  [../]
[]


[BCs]

  [./potential_grounded_bot]
    type = DirichletBC
    variable = potential_int
    boundary = 'front back'
    value = 0
  [../]

  [./disp_x_bot]
    type = DirichletBC
    variable = disp_x
    boundary = 'front back'
    value = 0
  [../]
  [./disp_y_bot]
    type = DirichletBC
    variable = disp_y
    boundary = 'front back'
    value = 0
  [../]
  [./disp_z_bot]
    type = DirichletBC
    variable = disp_z
    boundary = 'front back'
    value = 0
  [../]

  [./Periodic]
    [./up_down_pbc]
      primary = top
      secondary = bottom
      translation = '0 -6 0'
      variable = 'polar_x polar_y polar_z potential_int'
    [../]
    [./left_right_pbc]
      primary = left
      secondary = right
      translation = '6 0 0'
      variable = 'polar_x polar_y polar_z potential_int'
    [../]
    #[./front_back_pbc]
    #  primary = front
    #  secondary = back
    #  translation = '0 0 -6'
    #  variable = 'polar_x polar_y polar_z potential_int'
    #[../]
  [../]

#  # fix center point location
#  [./centerfix_x]
#    type = PresetBC
#    boundary = 100
#    variable = disp_x
#    value = 0
#  [../]
#  [./centerfix_y]
#    type = PresetBC
#    boundary = 100
#    variable = disp_y
#    value = 0
#  [../]
#  [./centerfix_z]
#    type = PresetBC
#    boundary = 100
#    variable = disp_z
#    value = 0
#  [../]
#
#  # fix side point x coordinate to inhibit rotation
#  [./angularfix]
#    type = PresetBC
#    boundary = 101
#    variable = disp_x
#    value = 0
#  [../]
[]

# We monitor the total free energy and the total solute concentration (should be constant)
[Postprocessors]
 [./Fbulk]
      type = BulkEnergy
      block = '0'
      execute_on = 'timestep_end'
    [../]
    [./Fwall]
      type = WallEnergy
      block = '0'
      execute_on = 'timestep_end'
    [../]
    [./Felastic]
      type = ElasticEnergy
      block = '0'
      execute_on = 'timestep_end'
    [../]
    [./Fcoupled]
      block = '0'
      type = CoupledEnergy
      execute_on = 'timestep_end'
    [../]
    [./Felec]
      block = '0'
      type = ElectrostaticEnergy
      permittivity = 0.08854187
      execute_on = 'timestep_end'
    [../]
    [./Ftotal]
      type = TotalEnergyFlow
      Fbulk = Fbulk
      Fwall = Fwall
      Fcoupled = Fcoupled
      Felec = Felec
      execute_on = 'timestep_end'
    [../]
    [./perc_change]
     type = PercentChangePostprocessor
     postprocessor = Ftotal
   [../]
[]

[UserObjects]
 [./kill]
  type = Terminator
  expression = 'perc_change <= 1.0e-3'
 [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type -pc_factor_shift_type -pc_factor_shift_amount'
    petsc_options_value = '    121                1e-6      1e-7    lu      NONZERO               1e-10'
  [../]
[]

[Executioner]
  type = Transient
    [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.7
    #iteration_window = 3
    optimal_iterations = 6 #should be 5 probably
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.8
[../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 0.7
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = outPTO_mortar_gradients_short
    elemental_as_nodal = true
    interval = 5
  [../]
  [./outcsv]
    type = CSV
    file_base = outPTO_mortar_gradients_short
  [../]
[]
