[Mesh]
 file = embedded_single_sphere_4.e
[]

[GlobalParams]
  len_scale = 1.0
  alpha1 = -0.1722883 # (3.766(T-765.1)*10^5) C^{-2} nm^2 (T = 673 K)
  alpha11 = -0.07253
  alpha111 = 0.26
  alpha12 = 0.75
  alpha112 = 0.61
  alpha123 = -3.67
  G110 = 0.13
  G11/G110 = 0.6
  G12/G110 = 0
  G44/G110 = 0.3
  G44P/G110 = 0.3
  Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.03375 0.03375 0.03375'
#  permittivity = 0.5843763
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
  #potential_ext = potential_ext
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  displacements = 'disp_x disp_y disp_z'
  #initial_from_file_timestep = 120
[]


[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block = '1'
   # initial_from_file_var = polar_x
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block = '1'
  # initial_from_file_var = polar_y
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
    #initial_from_file_var = polar_y
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
    #initial_from_file_var = polar_z
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
    #initial_from_file_var = potential_int
  [../]
  #[./potential_ext]
  #  order=FIRST
  #  family = LAGRANGE
  #  #initial_from_file_var = potential_ext
  #[../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
    block = '1 2'
    #initial_from_file_var = disp_x
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
    block = '1 2'
    #scaling = 1e6
    #initial_from_file_var = disp_y
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
    #initial_from_file_var = disp_z
  [../]
[]

[AuxVariables]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_xx
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_yy
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_xy
  [../]
  [./stress_xz]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_xz
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_zz
  [../]
  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = stress_yz
  [../]
  [./elastic_strain_xx]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_xx
  [../]
  [./elastic_strain_yy]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_yy
  [../]
  [./elastic_strain_xy]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_xy
  [../]
  [./elastic_strain_xz]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_xz
  [../]
  [./elastic_strain_zz]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_zz
  [../]
  [./elastic_strain_yz]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_yz
  [../]
  [./curlmag]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./curlP]
    type = CurlP
    variable = curlmag
    execute_on = 'timestep_end'
  [../]
  [./matl_s11]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    block = '1 2'
    variable = stress_xx
    execute_on = 'timestep_end'
  [../]
  [./matl_s12]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    block = '1 2'
    variable = stress_xy
    execute_on = 'timestep_end'
  [../]
  [./matl_s13]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
    block = '1 2'
    variable = stress_xz
    execute_on = 'timestep_end'
  [../]
 [./matl_s22]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    block = '1 2'
    variable = stress_yy
    execute_on = 'timestep_end'
  [../]
  [./matl_s23]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    block = '1 2'
    variable = stress_yz
    execute_on = 'timestep_end'
  [../]
  [./matl_s33]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    block = '1 2'
    variable = stress_zz
    execute_on = 'timestep_end'
  [../]
  [./matl_e11]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    block = '1 2'
    variable = elastic_strain_xx
    execute_on = 'timestep_end'
  [../]
  [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    block = '1 2'
    variable = elastic_strain_xy
    execute_on = 'timestep_end'
  [../]
  [./matl_e13]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 2
    block = '1 2'
    variable = elastic_strain_xz
    execute_on = 'timestep_end'
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    block = '1 2'
    variable = elastic_strain_yy
    execute_on = 'timestep_end'
  [../]
  [./matl_e23]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 2
    block = '1 2'
    variable = elastic_strain_yz
    execute_on = 'timestep_end'
  [../]
  [./matl_e33]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 2
    index_j = 2
    block = '1 2'
    variable = elastic_strain_zz
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
     type=WallEnergyDerivative
     variable = polar_x
     component = 0
  [../]
 [./walled_y]
     type=WallEnergyDerivative
     variable = polar_y
     component = 1
  [../]
  [./walled_z]
     type=WallEnergyDerivative
     variable = polar_z
     component = 2
  [../]
  ##Polarization-strain coupling
  [./ferroelectriccouplingu_x]
     type = FerroelectricCouplingX
     variable = disp_x
     component = 0
     block = '1'
  [../]
  [./ferroelectriccouplingu_y]
     type = FerroelectricCouplingX
     variable = disp_y
     component = 1
     block = '1'
  [../]
  [./ferroelectriccouplingu_z]
     type = FerroelectricCouplingX
     variable = disp_z
     component = 2
     block = '1'
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
     permittivity = 0.5843763
     block = '1'
  [../]
  [./FE_E_int]
     type=Electrostatics
     variable = potential_int
     block = '1'
     permittivity = 0.5843763
  [../]
  [./DIE_E_int]
     type=Electrostatics
     variable = potential_int
     block = '2'
     permittivity = 175.31289
  [../]
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
    # time_scale = 1.8e-7
    time_scale = 1.0
  [../]
  [./polar_y_time]
     type=TimeDerivativeScaled
     variable=polar_y
    # time_scale = 1.8e-7
    time_scale = 1.0
  [../]
  [./polar_z_time]
     type=TimeDerivativeScaled
     variable = polar_z
    # time_scale = 1.8e-7
    time_scale = 1.0
  [../]
[]


[BCs]
  [./disp_x_1]
    type = DirichletBC
    variable = disp_x
    boundary = '1'
    value = 0
  [../]
  [./disp_y_1]
    type = DirichletBC
    variable = disp_y
    boundary = '1'
    value = 0
  [../]
  [./disp_z_1]
    type = DirichletBC
    variable = disp_z
    boundary = '1'
    value = 0
  [../]
  [./disp_x_2]
    type = DirichletBC
    variable = disp_x
    boundary = '2'
    value = 0
  [../]
  [./disp_y_2]
    type = DirichletBC
    variable = disp_y
    boundary = '2'
    value = 0
  [../]
  [./disp_z_2]
    type = DirichletBC
    variable = disp_z
    boundary = '2'
    value = 0
  [../]
  [./potential_int_1]
    type = DirichletBC
    variable = potential_int
    boundary = '1'
    value = 0.0001
  [../]
  [./potential_int_2]
    type = DirichletBC
    variable = potential_int
    boundary = '2'
    value = 0.0001
  [../]
[]

[Materials]
  [./spheres_ferroelectric]
    type=LinearFerroelectricMaterial
    block = '1'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
    #in GPA. from N. Pandech et al. Ceramic. Internat., times 1e9 to convert to N/m^2
    # C11 C12 C13 C22 C23 C33 C44 C55 C66
    #C_ijkl = '3.80e-7 1.5e-7 1.50e-7 3.80e-7 1.50e-7 3.80e-7 1.1e-7 1.1e-7 1.1e-7'
    #in m^4/C^2. from http://arxiv.org/pdf/1205.5640.pdf
    # Q11 Q12 Q13 Q22 Q23 Q33 Q44 Q55 Q66
    #euler_angle_1 = 0.0 #currently will only rotate C_ijkl
    #euler_angle_2 = 0.0
    #euler_angle_3 = 0.0
  [../]
  [./STO_dielectric]
    type=LinearElasticMaterial
    block = '2'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    #in GPA. from Materials Project.org -- chose P3mm form?
    # C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '319 99.6 99.6 319 99.6 319 109.53 109.53 109.53'


  [../]
  [./elasticity_tensor1]
    type = ComputeElasticityTensor
    block = '1'
    fill_method = symmetric9
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
  [../]
  [./strain1]
    type = ComputeSmallStrain
    block = '1'
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
  [../]
  [./stress1]
    type = ComputeLinearElasticStress
    block = '1'
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
  [../]
  [./elasticity_tensor2]
    type = ComputeElasticityTensor
    block = '2'
    fill_method = symmetric9
    C_ijkl = '319 99.6 99.6 319 99.6 319 109.53 109.53 109.53'
  [../]
  [./strain2]
    type = ComputeSmallStrain
    block = '2'
  [../]
  [./stress2]
    type = ComputeLinearElasticStress
    block = '2'
  [../]
[]


[Postprocessors]
  [./bulk_energy]
   type = BulkEnergy
   execute_on = 'timestep_end'
  # initial_from_file_var = bulk_energy
   block = '1'
  [../]
  [./wall_energy]
   type = WallEnergy
   execute_on = 'timestep_end'
  # initial_from_file_var = wall_energy
   block = '1'
  [../]
  [./elastic_energy]
   type = ElasticEnergy
   execute_on = 'timestep_end'
  # initial_from_file_var = elastic_energy
   block = '1 2'
  [../]
  [./coupled_energy]
    type = CoupledEnergy
    execute_on = 'timestep_end'
    block = '1'
  [../]
  [./electrostatic_energy]
   type = ElectrostaticEnergy
   execute_on = 'timestep_end'
   block = '1 2'
   permittivity = 0.5843763
  # initial_from_file_var = electrostatic_energy
  [../]
 #[./total_energy]
 #  type = TotalEnergy
 #  bulk_energy = bulk_energy
 #  wall_energy = wall_energy
 #  bulk_energy_fourth = bulk_energy_fourth
 #  elastic_energy = elastic_energy
 #  coupled_energy = coupled_energy
 #  electrostatic_energy = electrostatic_energy
 #  execute_on = 'timestep_end'
 # # initial_from_file_var = total_energy
 # [../]
  [./total_energy_noelastic]
   type = TotalEnergyFlow
   bulk_energy = bulk_energy
   wall_energy = wall_energy
   bulk_energy_fourth = bulk_energy_fourth
   coupled_energy = coupled_energy
   electrostatic_energy = electrostatic_energy
   execute_on = 'timestep_end'
  # initial_from_file_var = total_energy
  [../]
  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = total_energy_noelastic
  [../]
  [./|R(i)|]
    type = Residual
  [../]
  [./dt]
    type = TimestepSize
  [../]
[]

[UserObjects]
  [./kill]
    type = Terminator
    expression = 'perc_change <= 4.0e-5'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -options_left '
    petsc_options_iname = '-ksp_gmres_restart -snes_rtol -ksp_rtol -pc_type   -pc_asm_overlap -sub_pc_type    -sub_pc_factor_zeropivot -pc_factor_zeropivot -pc_side '
    petsc_options_value = '    121            1e-8	 1e-12      gamg        7   ilu          1e-50    1e-50	  left        '
  [../]
[]

[Executioner]
  type = Transient
  [./TimeStepper] #iterative DT halfs the time it takes to find a solution? oh well, our time is fake in this simulation anyway...
    type = IterationAdaptiveDT
    dt = 0.25 #max seems to be about 1.0 but could depend on refinement...
    #there is also a cutback on this for 0.2*optimal and yes i think it does count the 0th one.
    #iteration_window = 10
    optimal_iterations = 3 #i think this is 3 or more then cut? less than 3 grow, does it count the 0th iteration? no the cutting has to do with the iteration ratio
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    #linear_iteration_ratio = 1000
    cutback_factor =  0.95
  [../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 0.95
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_PTOSTOcomposite_single_4
#    output_initial = true
    elemental_as_nodal = true
    interval = 12
  [../]
[]
