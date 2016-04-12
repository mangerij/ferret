
[Mesh]
  file = exodus_thinfilm.e
[]

[GlobalParams]
  len_scale = 1.0
  alpha1 = -0.1722883 # (3.766(T-765.1)*10^5) C^{-2} nm^2 (T = 673 K)
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
  Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.03375 0.03375 0.03375'
  #permittivity = 5.8
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
  #potential_ext = potential_ext
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  displacements = 'disp_x disp_y disp_z'
  #use_displaced_mesh = false
  #prefactor = 0.005 #negative = tension, positive = compression
[]



[Variables]
  [./polar_x]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./polar_y]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./polar_z]
    block = '1'
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
    block = '1 2'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
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
  #[./ferroelectriccouplingu_x]
  #  type = FerroelectricCouplingX
  #  variable = disp_x
  #  block = '1'
  #  component = 0
  #[../]
  #[./ferroelectriccouplingu_y]
  #  type = FerroelectricCouplingX
  #  variable = disp_y
  #  block = '1'
  #  component = 1
  #[../]
  #[./ferroelectriccouplingu_z]
  #  type = FerroelectricCouplingX
  #  variable = disp_z
  #  block = '1'
  #  component = 2
  #[../]
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
  ##Electrostatics
  [./polar_x_electric_E]
     type=PolarElectricEStrong
     variable = potential_int
     permittivity = 0.008812
  [../]
  [./FE_E_int]
     type=Electrostatics
     variable = potential_int
     permittivity = 0.008812
  [../]

#  [./vac_E_int]
#     type=Electrostatics
#     variable = potential_int
#     block = '2'
#     permittivity = 2.6
#  [../]

  [./polar_electric_px]
     type=PolarElectricPStrong
     variable = polar_x
     component = 0
  [../]
  [./polar_electric_py]
     type=PolarElectricPStrong
     variable = polar_y
     component = 1
  [../]
  [./polar_electric_pz]
     type=PolarElectricPStrong
     variable = polar_z
     component = 2
  [../]
  ##Time dependence
  [./polar_x_time]
     type=TimeDerivativeScaled
     variable=polar_x
    time_scale = 1.0
  [../]
  [./polar_y_time]
     type=TimeDerivativeScaled
     variable=polar_y
    time_scale = 1.0
  [../]
  [./polar_z_time]
     type=TimeDerivativeScaled
     variable = polar_z
    time_scale = 1.0
  [../]
[]


[BCs]

  #[./disp_x_1]
  #  type = DirichletBC
  #  variable = disp_x
  #  boundary = '1'
  #  value = 0.0
  #[../]
  #[./disp_y_1]
  #  type = DirichletBC
  #  variable = disp_y
  #  boundary = '1'
  #  value = 0.0
  #[../]
  #[./disp_z_1]
  #  type = DirichletBC
  #  variable = disp_z
  #  boundary = '1'
  #  value = 0.0
  #[../]

#[./potential_int_1]
#  type = DirichletBC
#  variable = potential_int
#  boundary = '1'
#  value = -0.0001
#[../]
#
#[./potential_int_2]
#  type = DirichletBC
#  variable = potential_int
#  boundary = '2'
#  value = -0.0001
#[../]

  [./bot_disp_x]
    variable = disp_x
    type = DirichletBC
    value = 0
    boundary = '7'
  [../]
  [./bot_disp_y]
    variable = disp_y
    type = DirichletBC
    value = 0
    boundary = '7'
  [../]
  [./bot_disp_z]
    variable = disp_z
    type = DirichletBC
    value = 0
    boundary = '7'
  [../]

 [./Periodic]
    [./TB_disp_x_pbc]
      variable = disp_x
      primary = '3'
      secondary = '5'
      translation = '0 40 0'
    [../]
    [./TB_disp_y_pbc]
      variable = disp_y
      primary = '3'
      secondary = '5'
      translation = '0 40 0'
    [../]
    [./TB_disp_z_pbc]
      variable = disp_z
      primary = '3'
      secondary = '5'
      translation = '0 40 0'
    [../]

    [./TB_polar_x_pbc]
      variable = polar_x
      primary = '3'
      secondary = '5'
      translation = '0 40 0'
    [../]
    [./TB_polar_y_pbc]
      variable = polar_y
      primary = '3'
      secondary = '5'
      translation = '0 40 0'
    [../]
    [./TB_polar_z_pbc]
      variable = polar_z
      primary = '3'
      secondary = '5'
      translation = '0 40 0'
    [../]
    [./TB_potential_int_pbc]
      variable = potential_int
      primary = '3'
      secondary = '5'
      translation = '0 40 0'
    [../]
  #
    [./TBsub_disp_x_pbc]
      variable = disp_x
      primary = '8'
      secondary = '10'
      translation = '0 40 0'
    [../]
    [./TBsub_disp_y_pbc]
      variable = disp_y
      primary = '8'
      secondary = '10'
      translation = '0 40 0'
    [../]
    [./TBsub_disp_z_pbc]
      variable = disp_z
      primary = '8'
      secondary = '10'
      translation = '0 40 0'
    [../]

    [./RL_disp_x_pbc]
      variable = disp_x
      primary = '4'
      secondary = '6'
      translation = '40 0 0'
    [../]
    [./RL_disp_y_pbc]
      variable = disp_y
      primary = '4'
      secondary = '6'
      translation = '40 0 0'
    [../]
    [./RL_disp_z_pbc]
      variable = disp_z
      primary = '4'
      secondary = '6'
      translation = '40 0 0'
    [../]

    [./RL_polar_x_pbc]
      variable = polar_x
      primary = '4'
      secondary = '6'
      translation = '45 0 0'
    [../]
    [./RL_polar_y_pbc]
      variable = polar_y
      primary = '4'
      secondary = '6'
      translation = '40 0 0'
    [../]
    [./RL_polar_z_pbc]
      variable = polar_z
      primary = '4'
      secondary = '6'
      translation = '40 0 0'
    [../]
    [./RL_potential_int_pbc]
      variable = potential_int
      primary = '4'
      secondary = '6'
      translation = '40 0 0'
    [../]

    [./RLsub_disp_x_pbc]
      variable = disp_x
      primary = '9'
      secondary = '11'
      translation = '40 0 0'
    [../]
    [./RLsub_disp_y_pbc]
      variable = disp_y
      primary = '9'
      secondary = '11'
      translation = '40 0 0'
    [../]
    [./RLsub_disp_z_pbc]
      variable = disp_z
      primary = '9'
      secondary = '11'
      translation = '40 0 0'
    [../]
  [../]
[]

[Materials]
 #  [./eigen_strain_xx_yy] #Use for stress-free strain (ie epitaxial)
 #   type = ComputeEigenstrain
 #   block = '1'
 #  # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
 #   eigen_base = '0 0 0 0 0 0 0 0 1'
 # [../]

  [./slab_ferroelectric]
    block = '1'
    type=LinearFerroelectricMaterial
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
  [../]

  [./PTO]
    type=LinearElasticMaterial
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
    block = '1'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
  [./STO]
    type=LinearElasticMaterial
    C_ijkl = '319 99.6 99.6 319 99.6 319 109.53 109.53 109.53'
    block = '2'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]

  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
    block = '1'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    block = '1'
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
    block = '1'
  [../]

  [./elasticity_tensor_2]
    type = ComputeElasticityTensor
    C_ijkl = '319 99.6 99.6 319 99.6 319 109.53 109.53 109.53'
    fill_method = symmetric9
    block = '2'
  [../]
  [./strain_2]
    type = ComputeSmallStrain
    block = '2'
  [../]
  [./stress_2]
    type = ComputeLinearElasticStress
    block = '2'
  [../]
[]
[Postprocessors]
#  [./volume]
#    type = VolumePostprocessor
#    block = '1'
#    use_displaced_mesh = true
#  [../]
  [./bulk_energy]
     type = BulkEnergy
     block = '1'
     execute_on = 'timestep_end'
    [../]
   [./wall_energy]
    type = WallEnergy
    block = '1'
    execute_on = 'timestep_end'
   [../]
    [./elastic_energy]
    type = ElasticEnergy
    block = '1'
    execute_on = 'timestep_end'
    [../]
    [./coupled_energy]
     block = '1'
    type = CoupledEnergy
    execute_on = 'timestep_end'
    [../]
    [./electrostatic_energy]
     block = '1'
     type = ElectrostaticEnergy
     execute_on = 'timestep_end'
     permittivity = 0.008812
    [../]
    [./total_energy_noelastic]
    type = TotalEnergyFlow
    bulk_energy = bulk_energy
    wall_energy = wall_energy
    bulk_energy_fourth = bulk_energy_fourth
    coupled_energy = coupled_energy
    electrostatic_energy = electrostatic_energy
    execute_on = 'timestep_end'
    [../]
    [./perc_change]
     type = PercentChangePostprocessor
     postprocessor = total_energy_noelastic
   [../]
  []
[]


[UserObjects]
 [./kill]
  type = Terminator
  expression = 'perc_change <= 5.0e-4'
 [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -options_left'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type '
    petsc_options_value = '    121                1e-8      1e-8     bjacobi'
  [../]
[]

#Limits exist on -snes_rtol =< 1e-10.

[Executioner]
  type = Transient
    [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.12
    #iteration_window = 3
    optimal_iterations = 5
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.95
[../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
#dt = 0.15
  dtmin = 1e-13
  dtmax = 0.3
#num_steps= 5000
  #num_steps = 2
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_PTO_thinfilm_unstr_nat
    elemental_as_nodal = true
    interval = 15
  [../]
[]
