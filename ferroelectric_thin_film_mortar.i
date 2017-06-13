#
# Eigenstrain with Mortar gradient periodicity
#

[Mesh]
  type = MortarPeriodicMesh
  dim = 3
  nx = 10
  ny = 10
  nz = 10
  xmin = -3.5
  xmax = 3.5
  ymin = -3.5
  ymax = 3.5
  zmin = -3
  zmax = 3
  periodic_directions = 'x y'

  [./MortarInterfaces]
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
  [../]
[]

[MeshModifiers]
  [./cnode]
    type = AddExtraNodeset
    coord = '0.0 0.0 0.0'
    new_boundary = 100
  [../]
  [./anode]
    type = AddExtraNodeset
    coord = '-3.5 -3.5 -3.0'
    new_boundary = 101
  [../]
[]


[GlobalParams]

  len_scale = 1.0
 
  # (3.766(T-765.1)*10^5) C^{-2} nm^2 (T = 293 K) = -0.288 to -0.1722883 to 0.05
  alpha1 = -0.287
  alpha11 = -0.07253
  alpha111 = 0.26
  alpha12 = 0.75
  alpha112 = 0.61
  alpha123 = -3.67

  G110 = 0.173
  G11/G110 = 0.6
  G12/G110 = 0
  G44/G110 = 0.3
  G44P/G110 = 0.3

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z

  #negative = tension, positive = compression

  prefactor = 0.02

  displacements = 'disp_x disp_y disp_z'
[]

[AuxVariables]
  [./stress_xx_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_xx
  [../]
  [./stress_yy_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_yy
  [../]
  [./stress_xy_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_xy
  [../]
  [./stress_xz_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_xz
  [../]
  [./stress_zz_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_zz
  [../]
  [./stress_yz_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_yz
  [../]
  [./strain_xx_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_xx
  [../]
  [./strain_yy_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_yy
  [../]
  [./strain_xy_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_xy
  [../]
  [./strain_xz_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_xz
  [../]
  [./strain_zz_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_zz
  [../]
  [./strain_yz_elastic]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_yz
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
  [../]
  [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    variable = strain_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e13]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 2
    variable = strain_xz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    variable = strain_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e23]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 2
    variable = strain_yz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e33]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 2
    index_j = 2
    variable = strain_zz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s11]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = stress_xx_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s12]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = stress_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s13]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
    variable = stress_xz_elastic
    execute_on = 'timestep_end'
  [../]
 [./matl_s22]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = stress_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s23]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    variable = stress_yz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s33]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    variable = stress_zz_elastic
    execute_on = 'timestep_end'
  [../]
[]


[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]

  # Lagrange multipliers for gradient component periodicity
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
[]

[Constraints]
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
[]

[Kernels]
  #Elastic problem
  [./TensorMechanics]
  #This is an action block
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
    variable = disp_x
    component = 0
  [../]
  [./ferroelectriccouplingX_yy]
    type = FerroelectricCouplingX
    variable = disp_y
    component = 1
  [../]
  [./ferroelectriccouplingX_zz]
    type = FerroelectricCouplingX
    variable = disp_z
    component = 2
  [../]
  ##Electrostatics
  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
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
 [./eigen_strain_zz] #Use for stress-free strain (ie epitaxial)
   type = ComputeEigenstrain
   # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
   eigen_base = '1 0 0 0 1 0 0 0 0'
   eigenstrain_name = eigenstrain
 [../]


  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '281 115.74 115.74 281 115.74 281 97.18 97.18 97.18'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    eigenstrain_names = eigenstrain
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
  [../]

  [./slab_ferroelectric]
    type = ComputeElectrostrictiveTensor
    Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.03375 0.03375 0.03375'
    C_ijkl = '281 115.74 115.74 281 115.74 281 97.18 97.18 97.18' #PTO from MaterialsProject
  [../]
[]

[BCs]
  [./potential_int_1]
    type = DirichletBC
    variable = potential_int
    boundary = 'front'
    value = -0.0001
  [../]

  [./potential_int_2]
    type = DirichletBC
    variable = potential_int
    boundary = 'back'
    value = -0.0001
  [../]

#  [./bot_disp_x]
#    variable = disp_x
#    type = DirichletBC
#    value = 0.0
#    boundary = 'back'
#  [../]
#  [./bot_disp_y]
#    variable = disp_y
#    type = DirichletBC
#    value = 0.0
#    boundary = 'back'
#  [../]
#  [./bot_disp_z]
#    variable = disp_z
#    type = DirichletBC
#    value = 0.0
#    boundary = 'back'
#  [../]


  # fix center point location
  [./centerfix_x]
    type = PresetBC
    boundary = 100
    variable = disp_x
   value = 0.0
  [../]
  [./centerfix_y]
    type = PresetBC
    boundary = 100
    variable = disp_y
    value = 0.0
  [../]
  [./centerfix_z]
    type = PresetBC
    boundary = 100
    variable = disp_z
    value = 0.0
  [../]

  # fix side point x coordinate to inhibit rotation
  [./angularfix_x]
    type = PresetBC
    boundary = 101
    variable = disp_x
    value = 0.0
  [../]
  [./angularfix_y]
    type = PresetBC
    boundary = 101
    variable = disp_y
    value = 0.0
  [../]
  [./angularfix_z]
    type = PresetBC
    boundary = 101
    variable = disp_z
    value = 0.0
  [../]

  [./Periodic]

    [./TB_polar_x_pbc]
      variable = polar_x
      primary = 'bottom'
      secondary = 'top'
      translation = '0 7 0'
    [../]
    [./TB_polar_y_pbc]
      variable = polar_y
      primary = 'bottom'
      secondary = 'top'
      translation = '0 7 0'
    [../]
    [./TB_polar_z_pbc]
      variable = polar_z
      primary = 'bottom'
      secondary = 'top'
      translation = '0 7 0'
    [../]
    [./TB_potential_int_pbc]
      variable = potential_int
      primary = 'bottom'
      secondary = 'top'
      translation = '0 7 0'
    [../]

    [./RL_polar_x_pbc]
      variable = polar_x
      primary = 'right'
      secondary = 'left'
      translation = '-7 0 0'
    [../]
    [./RL_polar_y_pbc]
      variable = polar_y
      primary = 'right'
      secondary = 'left'
      translation = '-7 0 0'
    [../]
    [./RL_polar_z_pbc]
      variable = polar_z
      primary = 'right'
      secondary = 'left'
      translation = '-7 0 0'
    [../]
    [./RL_potential_pbc]
      variable = potential_int
      primary = 'right'
      secondary = 'left'
      translation = '-7 0 0'
    [../]
  [../]
[]

[Postprocessors]
   [./Fbulk]
      type = BulkEnergy
      execute_on = 'timestep_end'
    [../]
    [./Fwall]
      type = WallEnergy
      execute_on = 'timestep_end'
    [../]
    [./Felastic]
      type = ElasticEnergy
      execute_on = 'timestep_end'
    [../]
    [./Fcoupled]
      type = CoupledEnergy
      execute_on = 'timestep_end'
    [../]
    [./Felec]
      type = ElectrostaticEnergy
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
  execute_on = 'timestep_end'
 [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason -ksp_converged_reason'
    # mortar currently does not support MPI parallelization
    petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount -ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol'
    petsc_options_value = '    lu        NONZERO                     1e-10                 100             1e-10      1e-8      1e-2 '
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = 'PJFNK'


  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.05
    #iteration_window = 3
    optimal_iterations = 6 #should be 5 probably
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.8
  [../]
  dtmax = 0.4
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_film_tension_pbc
    elemental_as_nodal = true
    interval = 1
    execute_on = 'timestep_end'
  [../]
[]
