[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 35
  ny = 35
  nz = 16
  xmin = -15
  xmax = 15
  ymin = -15
  ymax = 15
  zmin = -7
  zmax = 7
[]

[GlobalParams]
len_scale = 1
alpha1 = 0.00372834 # (3.766(T-765.1)*10^5) C^{-2} nm^2 (T = 673 K)
alpha11 = -0.07253
alpha111 = 0.26
alpha12 = 0.75
alpha112 = 0.61
alpha123 = -3.6999999999999997
G110 = 0.06
G11/G110 = 0.6
G12/G110 = 0
G44/G110 = 0.3
G44P/G110 = 0.3
Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.0675 0.0675 0.0675'
permittivity = 38.568838572
polar_x = polar_x
polar_y = polar_y
polar_z = polar_z
potential_int = potential_int
#potential_ext = potential_ext
disp_x=disp_x
disp_y=disp_y
disp_z=disp_z
displacements='disp_x disp_y disp_z'
#use_displaced_mesh = false
C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
#initial_from_file_timestep = 120
[]



[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block = '0'
   # initial_from_file_var = polar_x
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block = '0'
  # initial_from_file_var = polar_y
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
    #initial_from_file_var = polar_y
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
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
    [../]
    block = '0'
    #scaling = 1e6
    #initial_from_file_var = disp_x
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
    block = '0'
    #scaling = 1e6
    #initial_from_file_var = disp_y
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    block = '0'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
    #initial_from_file_var = disp_z
  [../]
[]

[AuxVariables]
  #[./stress_xx]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  #initial_from_file_var = stress_xx
  #[../]
  #[./stress_yy]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  #initial_from_file_var = stress_yy
  #[../]
  #[./stress_xy]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  #initial_from_file_var = stress_xy
  #[../]
  #[./stress_xz]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  #initial_from_file_var = stress_xz
  #[../]
  #[./stress_zz]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  #initial_from_file_var = stress_zz
  #[../]
  #[./stress_yz]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  #initial_from_file_var = stress_yz
  #[../]
  #[./elastic_strain_xx]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  #initial_from_file_var = strain_xx
  #[../]
  [./elastic_strain_yy]
    order = CONSTANT
    family = MONOMIAL
    #initial_from_file_var = strain_yy
  [../]
  #[./elastic_strain_xy]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  #initial_from_file_var = strain_xy
  #[../]
  #[./elastic_strain_xz]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  #initial_from_file_var = strain_xz
  #[../]
  #[./elastic_strain_zz]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  #initial_from_file_var = strain_zz
  #[../]
  #[./elastic_strain_yz]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  #initial_from_file_var = strain_yz
  #[../]
  #[./rho_b]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  #initial_from_file_var = rho_b
  #[../]
#  [./screenfield]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./surf_charge]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
[]
#
[AuxKernels]
#  [./surfacecharge]
#    type = SurfaceChargeAux
#    variable = surf_charge
#    boundary = '3'
#  [../]
#  [./surfacescreenfield]
#    type = ScreenAux
#    variable = screenfield
#    boundary = '3 4'
#  [../]
#  [./boundcharge]
#    type = BoundCharge
#    variable = rho_b
#    block = '2'
#    execute_on = 'timestep_end'
#  [../]
  #[./matl_s11]
  #  type = RankTwoAux
  #  rank_two_tensor = stress
  #  index_i = 0
  #  index_j = 0
  #  block = '0'
  #  variable = stress_xx
  #  execute_on = 'timestep_end'
  #[../]
  #[./matl_s12]
  #  type = RankTwoAux
  #  rank_two_tensor = stress
  #  index_i = 0
  #  index_j = 1
  #  block = '0'
  #  variable = stress_xy
  #  execute_on = 'timestep_end'
  #[../]
  #[./matl_s13]
  #  type = RankTwoAux
  #  rank_two_tensor = stress
  #  index_i = 0
  #  index_j = 2
  #  block = '0'
  #  variable = stress_xz
  #  execute_on = 'timestep_end'
  #[../]
  #[./matl_s22]
  #  type = RankTwoAux
  #  rank_two_tensor = stress
  #  index_i = 1
  #  index_j = 1
  #  block = '0'
  #  variable = stress_yy
  #  execute_on = 'timestep_end'
  #[../]
  #[./matl_s23]
  #  type = RankTwoAux
  #  rank_two_tensor = stress
  #  index_i = 1
  #  index_j = 2
  #  block = '0'
  #  variable = stress_yz
  #  execute_on = 'timestep_end'
  #[../]
  #[./matl_s33]
  #  type = RankTwoAux
  #  rank_two_tensor = stress
  #  index_i = 2
  #  index_j = 2
  #  block = '0'
  #  variable = stress_zz
  #  execute_on = 'timestep_end'
  #[../]
  #[./matl_e11]
  #  type = RankTwoAux
  #  rank_two_tensor = elastic_strain
  #  index_i = 0
  #  index_j = 0
  #  block = '0'
  #  variable = elastic_strain_xx
  #  execute_on = 'timestep_end'
  #[../]
  #[./matl_e12]
  #  type = RankTwoAux
  #  rank_two_tensor = elastic_strain
  #  index_i = 0
  #  index_j = 1
  #  block = '0'
  #  variable = elastic_strain_xy
  #  execute_on = 'timestep_end'
  #[../]
  #[./matl_e13]
  #  type = RankTwoAux
  #  rank_two_tensor = elastic_strain
  #  index_i = 0
  #  index_j = 2
  #  block = '0'
  #  variable = elastic_strain_xz
  #  execute_on = 'timestep_end'
  #[../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    block = '0'
    variable = elastic_strain_yy
    execute_on = 'timestep_end'
  [../]
  #[./matl_e23]
  #  type = RankTwoAux
  #  rank_two_tensor = elastic_strain
  #  index_i = 1
  #  index_j = 2
  #  block = '0'
  #  variable = elastic_strain_yz
  #  execute_on = 'timestep_end'
  #[../]
  #[./matl_e33]
  #  type = RankTwoAux
  #  rank_two_tensor = elastic_strain
  #  index_i = 2
  #  index_j = 2
  #  block = '0'
  #  variable = elastic_strain_zz
  #  execute_on = 'timestep_end'
  #[../]
[]

[Kernels]
  #Elastic problem
  [./TensorMechanics] #This is an action block
  [../]

  #[./stressdiv_0]
  #  type = StressDivergenceTensorsScaled
  #  variable = disp_x
  #  component = 0
  #  block = '0'
  #[../]
  #[./stressdiv_1]
  #  type = StressDivergenceTensorsScaled
  #  variable = disp_y
  #  component = 1
  #  block = '0'
  #[../]
  #[./stressdiv_2]
  #  type = StressDivergenceTensorsScaled
  #  variable = disp_z
  #  component = 2
  #  block = '0'
  #[../]

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
  #[./ferroelectriccouplingu_x]
  #   type = FerroelectricCouplingU
  #   variable = disp_x
  #   component = 0
  #   block = '2'
  #[../]
  #[./ferroelectriccouplingu_y]
  #   type = FerroelectricCouplingU
  #   variable = disp_y
  #   component = 1
  #   block = '2'
  #[../]
  #[./ferroelectriccouplingu_z]
  #   type = FerroelectricCouplingU
  #   variable=disp_z
  #   component = 2
  #   block = '2'
  #[../]
  [./ferroelectriccouplingp_xx]
     type = FerroelectricCouplingP
     variable=polar_x
     component = 0
    # block = '2'
  [../]
  [./ferroelectriccouplingp_yy]
     type = FerroelectricCouplingP
     variable = polar_y
     component = 1
    # block = '2'
  [../]
  [./ferroelectriccouplingp_zz]
     type = FerroelectricCouplingP
     variable = polar_z
     component = 2
    # block = '2'
  [../]
  ##Electrostatics
  [./polar_x_electric_E]
     type=PolarElectricEStrong
     variable = potential_int
     block = '0'
    # implicit = false
  [../]
  #[./E_int]
  #   type=Electrostatics
  #   variable = potential_int
  #   block='1 3'
  #[../]
  [./FE_E_int]
     type=Electrostatics
     variable = potential_int
     block = '0'
  [../]
  #[./E_ext]
  #   type=Electrostatics
  #   variable = potential_ext
  #   block='1 3'
  #[../]
  #[./FE_E_ext]
  #   type=Electrostatics
  #   variable = potential_ext
  #   block='2'
  #[../]
  [./polar_electric_px]
     type=PolarElectricPStrong
     variable = polar_x
     component = 0
    # implicit = false
  [../]
  [./polar_electric_py]
     type=PolarElectricPStrong
     variable = polar_y
     component = 1
    # implicit = false
  [../]
  [./polar_electric_pz]
     type=PolarElectricPStrong
     variable = polar_z
     component = 2
    # implicit = false
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
  #[./potential_int_time]
  #   type=TimeDerivativeScaled
  #   variable = potential_int
  #   time_scale = 1.0
  #[../]
  #[./disp_x_time]
  #   type=TimeDerivativeScaled
  #   variable = disp_x
  #   time_scale = 1.0e-9 #this is chosen so that the elastic problem is solved in ~4 timesteps.
  #[../]
  #[./disp_y_time]
  #   type=TimeDerivativeScaled
  #   variable = disp_y
  #   time_scale = 1.0e-9
  #[../]
  #[./disp_z_time]
  #   type=TimeDerivativeScaled
  #   variable = disp_z
  #   time_scale = 1.0e-9
  #[../]
[]

[BCs]

  [./disp_x_slab5]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0.32
  [../]
  [./disp_x_slab7]
    type = DirichletBC
    variable = disp_x
    boundary = 'right'
    value = -0.32
  [../]

  [./disp_y_slab5]
    type = DirichletBC
    variable = disp_y
    boundary = 'top'
    value = -0.32
  [../]
  [./disp_y_slab7]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0.32
  [../]


  [./potential_cube5]
    type = DirichletBC
    boundary = 'front'
    value = 0.000002
    variable = potential_int
  [../]
  [./potential_cube6]
    type = DirichletBC
    boundary = 'back'
    value = 0.000002
    variable = potential_int
  [../]

[]

[Materials]
  [./slab_ferroelectric]
    type=LinearFerroelectricMaterial
    block = '0'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    #in GPA. from N. Pandech et al. Ceramic. Internat., times 1e9 to convert to N/m^2
    # C11 C12 C13 C22 C23 C33 C44 C55 C66
    #C_ijkl = '3.80e-7 1.5e-7 1.50e-7 3.80e-7 1.50e-7 3.80e-7 1.1e-7 1.1e-7 1.1e-7'
    #in m^4/C^2. from http://arxiv.org/pdf/1205.5640.pdf
    # Q11 Q12 Q13 Q22 Q23 Q33 Q44 Q55 Q66
    #euler_angle_1 = 0.0 #currently will only rotate C_ijkl
    #euler_angle_2 = 0.0
    #euler_angle_3 = 0.0
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = '0'
  #  block = '2'
    fill_method = symmetric9
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = '0'
  #  block = '2'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '0'
  #  block = '2'
  [../]
  #[./vacuum]
  #  type=GenericConstantMaterial
  #  block = '1 3'
  #[../]
  ##This seems to be what we want for a simple epitaxial test
  ## (note that most epitaxial conditions are a strain gradient from the interface)
  # Is this not seen by the simulation !!!!?!?!
  #[./eigen_strain_xx_yy] #Use for stress-free strain (ie epitaxial)
  #  type = ComputeEigenstrain
  #  block = '0'
  ##  block = '2'
  #  prefactor = -0.02
  #  # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
  #  eigen_base = '1 0 0 0 1 0 0 0 0'
  #[../]
[]

[Postprocessors]
  [./bulk_energy]
   type = BulkEnergy
   execute_on = 'timestep_end'
  # initial_from_file_var = bulk_energy
   block = '0'
  [../]
  [./wall_energy]
   type = WallEnergy
   execute_on = 'timestep_end'
  # initial_from_file_var = wall_energy
   block = '0'
  [../]
  [./elastic_energy]
   type = ElasticEnergy
   execute_on = 'timestep_end'
  # initial_from_file_var = elastic_energy
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
  # initial_from_file_var = electrostatic_energy
  [../]
  [./total_energy]
   type = TotalEnergy
   bulk_energy = bulk_energy
   wall_energy = wall_energy
   bulk_energy_fourth = bulk_energy_fourth
   elastic_energy = elastic_energy
   coupled_energy = coupled_energy
   electrostatic_energy = electrostatic_energy
   execute_on = 'timestep_end'
  # initial_from_file_var = total_energy
  [../]
  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = total_energy
  [../]
  [./|R(i)|]
    type = Residual
  [../]
[]

#[UserObjects]
#  [./kill]
#    type = Terminator
#    expression = 'perc_change <= 1.0e-15'
#  [../]
#[]


#[Preconditioning]
#  [./pbp]
#    type = PBP
#    #full = true
#    solve_order = 'potential_int disp_x disp_y disp_z polar_x polar_y polar_z '
#    preconditioner = 'ASM ASM ASM LU ASM ASM ASM'
#    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -options_left'
#    petsc_options_iname = '-ksp_gmres_restart -snes_rtol -ksp_rtol -pc_type  -sub_pc_type -sub_pc_factor_zeropivot -pc_factor_zeropivot -pc_side'
#    petsc_options_value = '    501             1e-8     1e-8        bjacobi       lu                   1e-50                 1e-50         left'
#  [../]
#[]



[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -options_left'
    petsc_options_iname = '-ksp_gmres_restart -snes_rtol -ksp_rtol -pc_type    -pc_factor_zeropivot -pc_side '
    petsc_options_value = '    501             1e-8   1e-14      bjacobi          1e-50          left        '
  [../]
[]

#Limits exist on -snes_rtol =< 1e-10.

[Executioner]
  type = Transient
  #[./TimeStepper]
    #type = IterationAdaptiveDT
    #dt = 2.0 #max seems to be about 1.0 but could depend on refinement...
    #optimal_iterations = 1
    #growth_factor = 1.0001
    #linear_iteration_ratio = 1000
    #cutback_factor =  0.5
  #[../]
  #l_max_its = 8000
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dt = 0.5
  dtmin = 1e-11
  dtmax = 0.5
  num_steps = 120
  #splitting = 'ferretsplit'
  #petsc_options_iname ='-pc_type -pc_factor_zeropivot'
  #petsc_options_value = 'lu          1e-50 '
[]
#
#[Splits]
#  [./ferretsplit]
#    type = Split
#    splitting = 'elastic ferroelectric' #split to two subproblems
#    splitting_type = schur #schur split somewhat bugged right now (only serial)
#    schur_type = full
#    schur_pre = A11
#  [../]
#  [./ferroelectric]
#    vars = 'polar_x polar_y polar_z potential_int'
#    petsc_options='-dm_view -ksp_monitor -inner_ksp_monitor' #'-ksp_monitor -inner_ksp_monitor'
#    petsc_options_iname=' -ksp_type   -ksp_gmres_restart  -ksp_rtol -inner_pc_type -pc_type   -pc_asm_overlap -inner_pc_factor_zeropivot  -pc_factor_zeropivot -pc_hypre_type -inner_pc_hypre_type'
#    petsc_options_value='    gmres            350              1e-8  lu        lu        2                         1e-50                1e-50  boomeramg boomeramg'
#  [../]
#  [./elastic]
#    vars = 'disp_x disp_y disp_z'
#    petsc_options='-dm_view -ksp_monitor -inner_ksp_monitor'  # ''-ksp_monitor'
#    petsc_options_iname=' -ksp_type  -ksp_gmres_restart -inner_pc_type -ksp_rtol -pc_type   -pc_asm_overlap -inner_pc_factor_zeropivot -pc_factor_zeropivot -pc_hypre_type -inner_pc_hypre_type'
#    petsc_options_value = ' gmres       350                lu        1e-8       lu          2           1e-50                        1e-50 boomeramg boomeramg'
#  [../]
#[]

#[./Adaptivity]
#    [./Indicators]
#      [./indicator_x]
#        type = GradientJumpIndicator
#        variable = polar_x
#      [../]
#      [./indicator_y]
#        type = GradientJumpIndicator
#        variable = polar_y
#      [../]
#      [./indicator_z]
#        type = GradientJumpIndicator
#        variable = polar_z
#      [../]
#    [../]
#[../]


[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_PbTiO3_50nm_scale-cnorm_T-9_strain-6
    output_initial = true
    elemental_as_nodal = true
    interval = 20
  [../]
  [./debug]
    type = VariableResidualNormsDebugOutput
  [../]
[]