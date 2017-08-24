
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 50
  ny = 50
  nz = 40
  ymax = 12
  ymin = -12
  xmin = -12
  xmax = 12
  zmin = -10
  zmax = 10
[]


[GlobalParams]
  #-------------------------------------------------------#
  #---------------------BiFeO3----------------------------#
  #---We zero the 6th order coefficients for now.---------#
  #-------------------------------------------------------#
  len_scale = 1.0
  alpha1 = -0.39347
  alpha11 = 0.542
  alpha111 = 0.0
  alpha12 = 0.154
  alpha112 = 0
  alpha123 = 0
  G110 = 0.98
  G11_G110 = 1.2
  G12/G110 = 0.0
  G44/G110 = 0.0
  G44P/G110 = 0.0
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
  #disp_x = disp_x
  #disp_y = disp_y
  #disp_z = disp_z
  #displacements = 'disp_x disp_y disp_z'
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
  #[./disp_x]
  #  order = FIRST
  #  family = LAGRANGE
  #[../]
  #[./disp_y]
  #  order = FIRST
  #  family = LAGRANGE
  #[../]
  #[./disp_z]
  #  order = FIRST
  #  family = LAGRANGE
  #[../]
[]

[AuxVariables]
  [./angle]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

#[AuxVariables]
#  [./stress_xx_elastic]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./stress_yy_elastic]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./stress_xy_elastic]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./stress_xz_elastic]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./stress_zz_elastic]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./stress_yz_elastic]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./strain_xx_elastic]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./strain_yy_elastic]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./strain_xy_elastic]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./strain_xz_elastic]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./strain_zz_elastic]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./strain_yz_elastic]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#[]
#

[AuxKernels]
  [./angle_phase]
    type = AngleAux
    block = '0'
    variable = angle
    execute_on = 'timestep_end'
  [../]
[]

#[AuxKernels]
#  [./matl_s11]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    index_i = 0
#    index_j = 0
#    block = '0'
#    variable = stress_xx_elastic
#    execute_on = 'timestep_end'
#  [../]
#  [./matl_s12]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    index_i = 0
#    index_j = 1
#    block = '0'
#    variable = stress_xy_elastic
#    execute_on = 'timestep_end'
#  [../]
#  [./matl_s13]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    index_i = 0
#    index_j = 2
#    block = '0'
#    variable = stress_xz_elastic
#    execute_on = 'timestep_end'
#  [../]
# [./matl_s22]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    index_i = 1
#    index_j = 1
#    block = '0'
#    variable = stress_yy_elastic
#    execute_on = 'timestep_end'
#  [../]
#  [./matl_s23]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    index_i = 1
#    index_j = 2
#    block = '0'
#    variable = stress_yz_elastic
#    execute_on = 'timestep_end'
#  [../]
#  [./matl_s33]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    index_i = 2
#    index_j = 2
#    block = '0'
#    variable = stress_zz_elastic
#    execute_on = 'timestep_end'
#  [../]
#  [./matl_e11]
#    type = RankTwoAux
#    rank_two_tensor = elastic_strain
#    index_i = 0
#    index_j = 0
#    block = '0'
#    variable = strain_xx_elastic
#    execute_on = 'timestep_end'
#  [../]
#  [./matl_e12]
#    type = RankTwoAux
#    rank_two_tensor = elastic_strain
#    index_i = 0
#    index_j = 1
#    block = '0'
#    variable = strain_xy_elastic
#    execute_on = 'timestep_end'
#  [../]
#  [./matl_e13]
#    type = RankTwoAux
#    rank_two_tensor = elastic_strain
#    index_i = 0
#    index_j = 2
#    block = '0'
#    variable = strain_xz_elastic
#    execute_on = 'timestep_end'
#  [../]
#  [./matl_e22]
#    type = RankTwoAux
#    rank_two_tensor = elastic_strain
#    index_i = 1
#    index_j = 1
#    block = '0'
#    variable = strain_yy_elastic
#    execute_on = 'timestep_end'
#  [../]
#  [./matl_e23]
#    type = RankTwoAux
#    rank_two_tensor = elastic_strain
#    index_i = 1
#    index_j = 2
#    block = '0'
#    variable = strain_yz_elastic
#    execute_on = 'timestep_end'
#  [../]
#  [./matl_e33]
#    type = RankTwoAux
#    rank_two_tensor = elastic_strain
#    index_i = 2
#    index_j = 2
#    block = '0'
#    variable = strain_zz_elastic
#    execute_on = 'timestep_end'
#  [../]
#[]

[Kernels]
  ##-------------------------#
  ##-----Elastic problem-----#
  ##-------------------------#
  #
  ##Stress-divergence eq.
  #[./TensorMechanics]
  ##This is an action block
  #[../]

  #------------------------#
  #--Thin film subproblem--#
  #------------------------#
  # NOTE: this is omitted  #
  #       for now          #


  #------------------------#
  #-Ferroelectric problem -#
  #------------------------#

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
###Polarization-strain coupling
#
#  [./ferroelectriccouplingp_xx]
#    type = FerroelectricCouplingP
#    variable = polar_x
#    component = 0
#  [../]
#  [./ferroelectriccouplingp_yy]
#    type = FerroelectricCouplingP
#    variable = polar_y
#    component = 1
#  [../]
#  [./ferroelectriccouplingp_zz]
#    type = FerroelectricCouplingP
#    variable = polar_z
#    component = 2
#  [../]
#  [./ferroelectriccouplingX_xx]
#    type = FerroelectricCouplingX
#    block = '0'
#    variable = disp_x
#    component = 0
#  [../]
#  [./ferroelectriccouplingX_yy]
#    type = FerroelectricCouplingX
#    block = '0'
#    variable = disp_y
#    component = 1
#  [../]
#  [./ferroelectriccouplingX_zz]
#    type = FerroelectricCouplingX
#    block = '0'
#    variable = disp_z
#    component = 2
#  [../]
  ##Electrostatics
  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     block = '0'
     permittivity = 0.1679
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

[BCs]
[./potential_int_short]
  type = DirichletBC
  variable = potential_int
  boundary = 'front back'
  value = -0.00001
[../]
  #
  #[./bot_disp_zero]
  #  variable = 'disp_x disp_y disp_z'
  #  type = DirichletBC
  #  value = 0
  #  boundary = 'back'
  #[../]
[]

#[Materials]
#  [./elasticity_tensor_1]
#    block = '0'
#    type = ComputeElasticityTensor
#    fill_method = symmetric21
#    # From Phys. Rev. B. 80, 052102 (2009)
#    # Note this is rhombohedral
#    # C11 C12 C13 C14 C15 C16 C22 C23 C24 C25 C26 C33 C34 C35 C36 C44 C45 C46 C55 C56 C66
#    C_ijkl = '213 111 49 19 0 0 213 49 -19 0 0 139 0 0 0 31 0 0 31 0 31'
#  [../]
#  [./strain_1]
#    block = '0'
#    type = ComputeSmallStrain
#  [../]
#  [./stress_1]
#    block = '0'
#    type = ComputeLinearElasticStress
#  [../]
#  [./slab_ferroelectric]
#    block = '0'
#    #from a Chen et al J. Appl. Phys. 103, 094111 (2008)
#    type = ComputeElectrostrictiveTensor
#    fill_methodQ = symmetric9
#    fill_methodC = symmetric21
#    # C11 C12 C13 C14 C15 C16 C22 C23 C24 C25 C26 C33 C34 C35 C36 C44 C45 C46 C55 C56 C66
#    Q_mnkl = '0.032 -0.016 -0.016 0.032 -0.016 0.032 0.01 0.01 0.01'
#    C_ijkl = '213 111 49 19 0 0 213 49 -19 0 0 139 0 0 0 31 0 0 31 0 31'
#  [../]
#[]

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
    #[./Felastic]
    #  type = ElasticEnergy
    #  block = '0'
    #  execute_on = 'timestep_end'
    #[../]
    #[./Fcoupled]
    #  block = '0'
    #  type = CoupledEnergy
    #  execute_on = 'timestep_end'
    #[../]
    [./Felec]
      block = '0'
      type = ElectrostaticEnergy
      permittivity = 0.167
      execute_on = 'timestep_end'
    [../]
    [./Ftotal]
      type = TotalEnergyFlowNoElast
      Fbulk = Fbulk
      Fwall = Fwall
      #Fcoupled = Fcoupled
      Felec = Felec
      execute_on = 'timestep_end'
    [../]
    [./perc_change]
     type = PercentChangePostprocessor
     postprocessor = Ftotal
   [../]
   [./num_NLin]
    type = NumNonlinearIterations
   [../]
   [./num_Lin]
    type = NumLinearIterations
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
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -options_left'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121                1e-8      1e-8    bjacobi'
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
  dtmax = 0.8
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_BFO_test_noElast
    output_initial = true
    elemental_as_nodal = true
    interval = 1
  [../]
[]
