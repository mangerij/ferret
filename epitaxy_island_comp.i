
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 30
  ny = 30
  nz = 12
  ymax = 10
  ymin = -10
  xmin = -10
  xmax = 10
  zmin = -4
  zmax = 4
[]


[GlobalParams]
  len_scale = 1.0
  alpha1 = -0.027721
  alpha11 = -0.64755
  alpha111 = 8.004
  alpha12 = 0.323
  alpha112 = 4.47
  alpha123 = 4.919
  G110 = 1.0
  G11/G110 = 0.51
  G12/G110 = 0.0
  G44/G110 = 0.01
  G44P/G110 = 0.01
  Q_mnkl = '0.11 -0.045 -0.045 0.11 -0.045 0.11 0.059 0.059 0.059'
  permittivity = 4.06
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
  C_ijkl = '178 96.4 96.4 178 96.4 178 122 122 122'
  prefactor = 0.005 #negative = tension, positive = compression
[]



[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-4
      max = 0.5e-4
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-4
      max = 0.5e-4
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-4
      max = 0.5e-4
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
    block = '0'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./disp_y]
    block = '0'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./disp_z]
    block = '0'
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
    block = '0'
    variable = stress_xx_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s12]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    block = '0'
    variable = stress_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s13]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
    block = '0'
    variable = stress_xz_elastic
    execute_on = 'timestep_end'
  [../]
 [./matl_s22]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    block = '0'
    variable = stress_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s23]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    block = '0'
    variable = stress_yz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s33]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    block = '0'
    variable = stress_zz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e11]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    block = '0'
    variable = strain_xx_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    block = '0'
    variable = strain_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e13]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 2
    block = '0'
    variable = strain_xz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    block = '0'
    variable = strain_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e23]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 2
    block = '0'
    variable = strain_yz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e33]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 2
    index_j = 2
    block = '0'
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
  [./ferroelectriccouplingu_x]
    type = FerroelectricCouplingX
    variable = disp_x
    component = 0
    block = '0'
  [../]
  [./ferroelectriccouplingu_y]
    type = FerroelectricCouplingX
    variable = disp_y
    component = 1
    block = '0'
  [../]
  [./ferroelectriccouplingu_z]
    type = FerroelectricCouplingX
    variable = disp_z
    component = 2
    block = '0'
  [../]
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
     block = '0'
  [../]
  [./FE_E_int]
     type=Electrostatics
     variable = potential_int
     block = '0'
  [../]

#  [./vac_E_int]
#     type=Electrostatics
#     variable = potential_int
#     block = '2'
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
#  [./top_disp_x]
#    variable = disp_x
#    type = DirichletBC
#    value = 0
#    boundary = 'front'
#  [../]
#  [./top_disp_y]
#    variable = disp_y
#    type = DirichletBC
#    value = 0
#    boundary = 'front'
#  [../]
#  [./top_disp_z]
#    variable = disp_z
#    type = DirichletBC
#    value = 0
#    boundary = 'front'
#  [../]
[./potential_int_1]
  type = DirichletBC
  variable = potential_int
  boundary = 'front'
  value = -0.00001
[../]

[./potential_int_2]
  type = DirichletBC
  variable = potential_int
  boundary = 'back'
  value = -0.00001
[../]

  [./bot_disp_x]
    variable = disp_x
    type = DirichletBC
    value = 0
    boundary = 'back'
  [../]
  [./bot_disp_y]
    variable = disp_y
    type = DirichletBC
    value = 0
    boundary = 'back'
  [../]
  [./bot_disp_z]
    variable = disp_z
    type = DirichletBC
    value = 0
    boundary = 'back'
  [../]
  #
  #[./top_disp_x]
  #  variable = disp_x
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'top'
  #[../]
  #[./top_disp_y]
  #  variable = disp_y
  #  type = DirichletBC
  #  value = -0.5
  #  boundary = 'top'
  #[../]
  #[./top_disp_z]
  #  variable = disp_z
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'top'
  #[../]
  #
  #[./bottom_disp_x]
  #  variable = disp_x
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'bottom'
  #[../]
  #[./bottom_disp_y]
  #  variable = disp_y
  #  type = DirichletBC
  #  value = 0.5
  #  boundary = 'bottom'
  #[../]
  #[./bottom_disp_z]
  #  variable = disp_z
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'bottom'
  #[../]
  #
  #
  #
  #
  #[./right_disp_x]
  #  variable = disp_x
  #  type = DirichletBC
  #  value = -0.5
  #  boundary = 'right'
  #[../]
  #[./right_disp_y]
  #  variable = disp_y
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'right'
  #[../]
  #[./right_disp_z]
  #  variable = disp_z
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'right'
  #[../]
  #
  #[./left_disp_x]
  #  variable = disp_x
  #  type = DirichletBC
  #  value = 0.5
  #  boundary = 'left'
  #[../]
  #[./left_disp_y]
  #  variable = disp_y
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'left'
  #[../]
  #[./left_disp_z]
  #  variable = disp_z
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'left'
  #[../]
  [./Periodic]
    [./TB_disp_x_pbc]
      variable = disp_x
      primary = 'top'
      secondary = 'bottom'
      translation = '0 -20 0'
    [../]
    [./TB_disp_y_pbc]
      variable = disp_y
      primary = 'top'
      secondary = 'bottom'
      translation = '0 -20 0'
    [../]
    [./TB_disp_z_pbc]
      variable = disp_z
      primary = 'top'
      secondary = 'bottom'
      translation = '0 -20 0'
    [../]

    [./TB_polar_x_pbc]
      variable = polar_x
      primary = 'top'
      secondary = 'bottom'
      translation = '0 -20 0'
    [../]
    [./TB_polar_y_pbc]
      variable = polar_y
      primary = 'top'
      secondary = 'bottom'
      translation = '0 -20 0'
    [../]
    [./TB_polar_z_pbc]
      variable = polar_z
      primary = 'top'
      secondary = 'bottom'
      translation = '0 -20 0'
    [../]
    [./TB_potential_int_pbc]
      variable = potential_int
      primary = 'top'
      secondary = 'bottom'
      translation = '0 -20 0'
    [../]

    [./RL_disp_x_pbc]
      variable = disp_x
      primary = 'right'
      secondary = 'left'
      translation = '-20 0 0'
    [../]
    [./RL_disp_y_pbc]
      variable = disp_y
      primary = 'right'
      secondary = 'left'
      translation = '-20 0 0'
    [../]
    [./RL_disp_z_pbc]
      variable = disp_z
      primary = 'right'
      secondary = 'left'
      translation = '-20 0 0'
    [../]

    [./RL_polar_x_pbc]
      variable = polar_x
      primary = 'right'
      secondary = 'left'
      translation = '-20 0 0'
    [../]
    [./RL_polar_y_pbc]
      variable = polar_y
      primary = 'right'
      secondary = 'left'
      translation = '-20 0 0'
    [../]
    [./RL_polar_z_pbc]
      variable = polar_z
      primary = 'right'
      secondary = 'left'
      translation = '-20 0 0'
    [../]
    [./RL_potential_int_pbc]
      variable = potential_int
      primary = 'right'
      secondary = 'left'
      translation = '-20 0 0'
    [../]
  [../]



[]

[Materials]
   [./eigen_strain_xx_yy] #Use for stress-free strain (ie epitaxial)
    type = ComputeEigenstrain
    block = '0'
   # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
    eigen_base = '1 0 0 0 1 0 0 0 0'
  [../]

  [./slab_ferroelectric]
    type=LinearFerroelectricMaterial
    block = '0'
  [../]

  [./BTO]
    type=LinearElasticMaterial
    block = '0'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = '0'
    fill_method = symmetric9
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = '0'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '0'
  [../]
[]

[Postprocessors]
  [./volume]
    type = VolumePostprocessor
    use_displaced_mesh = true
  [../]
  [./bulk_energy]
     type = BulkEnergy
     execute_on = 'timestep_end'
     block = '0'
    [../]
   [./wall_energy]
    type = WallEnergy
    execute_on = 'timestep_end'
    block = '0'
   [../]
    [./elastic_energy]
    type = ElasticEnergy
    execute_on = 'timestep_end'
    block = '0'
    [../]
    [./coupled_energy]
    type = CoupledEnergy
    execute_on = 'timestep_end'
    block = '0'
    [../]
    [./electrostatic_energy]
     type = ElectrostaticEnergy
     execute_on = 'timestep_end'
     block = '0'
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
 expression = 'perc_change <= 5.0e-5'
[../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -options_left'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121                1e-8      1e-8    gamg'
  [../]
[]

#Limits exist on -snes_rtol =< 1e-10.

[Executioner]
  type = Transient
    [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.8
    iteration_window = 10
    optimal_iterations = 5
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.85
[../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
#dt = 0.15
  dtmin = 1e-13
  dtmax = 1
#num_steps= 5000
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_BTO_coupling_slab_tens
    output_initial = true
    elemental_as_nodal = true
    interval = 2
  [../]
[]