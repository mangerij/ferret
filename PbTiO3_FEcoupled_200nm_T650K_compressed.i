
[Mesh]
  file = slab_exodus_coarse_200.e
  uniform_refine = 0
  #distribution = serial
[]

[GlobalParams]
  # #Can use a unit system where the differences between the coefficients
  # # in this table are minimized: [nm], [aC], [nN]
   #
  # #length scale
   len_scale = 1.0
  # #BulkEnergy coefficients
  # alpha1 = -0.03876 # (3.8 * (T-785) * 10^5) C^{-2} nm^2 (T = 650 K)
  # alpha11 = -0.073  #
  # alpha111 = 0.26
  # alpha12 = 0.75
  # alpha112 = 0.61
  # alpha123 = -3.7
  # #WallEnergy coefficients
  # G110 = 0.1
  # G11/G110 = 0.6
  # G12/G110 = 0.0
  # G44/G110 = 0.3
  # G44P/G110 = 0.3
  # #Electrostatics
   permittivity = 0.008854187
  # polar_x = polar_x
  # polar_y = polar_y
  # polar_z = polar_z
  # potential_int = potential_int
  # potential_ext = potential_ext
  # #elastic variables
   disp_x = disp_x
   disp_y = disp_y
   disp_z = disp_z
   displacements = 'disp_x disp_y disp_z'
   use_displaced_mesh = false

  # restart system
  # initial_from_file_timestep = 120
[]

[Variables]
  #[./polar_x]
  #  order = FIRST
  #  family = LAGRANGE
  #  block='2'
  #  [./InitialCondition]
  #    type = RandomIC
  #    min = -0.50
  #    max = 0.50
  #  [../]
  #  # initial_from_file_var = polar_x
  #[../]
  #[./polar_y]
  #  order = FIRST
  #  family = LAGRANGE
  #  block='2'
  #  [./InitialCondition]
  #    type = RandomIC
  #    min = -0.50
  #    max = 0.50
  #  [../]
  #  #initial_from_file_var = polar_y
  #[../]
  #[./polar_z]
  #  order = FIRST
  #  family = LAGRANGE
  #  block='2'
  #  #scaling = 1e5
  #  [./InitialCondition]
  #    type = RandomIC
  #    min = -0.50
  #    max = 0.50
  #  [../]
  #  #initial_from_file_var = polar_z
  #[../]
  [./potential_int]
    order=FIRST
    family = LAGRANGE
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
    block = '2'
    #initial_from_file_var = disp_x
    #[./InitialCondition] #Thought this was needed to get around the zero pivot
    #  type = RandomIC
    #  min = 0.1e-1
    #  max = 1.5e-1
    #[../]
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    block = '2'
    #initial_from_file_var = disp_y
    #[./InitialCondition]
    #  type = RandomIC
    #  min = 0.1e-1
    #  max = 1.5e-1
    #[../]
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    block = '2'
    #initial_from_file_var = disp_z
    #[./InitialCondition]
    #  type = RandomIC
    #  min = 0.1e-1
    #  max = 1.5e-1
    #[../]
  [../]
[]

#NOTE: AuxSystem causes a 20% slowdown in computation. Only use AuxVariables
#      if they are needed!

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
  #[./rho_b]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  initial_from_file_var = rho_b
  #[]
[]

[AuxKernels]
  #[./boundcharge]
  #  type = BoundCharge
  #  variable = rho_b
  #  execute_on = 'timestep_end'
  #  block = '2'
  #[../]
  [./matl_s11]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    block = '2'
    variable = stress_xx
    #execute_on = 'timestep_end'
  [../]
  [./matl_s12]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    block = '2'
    variable = stress_xy
    #execute_on = 'timestep_end'
  [../]
  [./matl_s13]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
    block = '2'
    variable = stress_xz
    #execute_on = 'timestep_end'
  [../]
  [./matl_s22]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    block = '2'
    variable = stress_yy
    #execute_on = 'timestep_end'
  [../]
  [./matl_s23]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    block = '2'
    variable = stress_yz
    #execute_on = 'timestep_end'
  [../]
  [./matl_s33]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    block = '2'
    variable = stress_zz
    #execute_on = 'timestep_end'
  [../]
  [./matl_e11]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    block = '2'
    variable = elastic_strain_xx
    #execute_on = 'timestep_end'
  [../]
  [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    block = '2'
    variable = elastic_strain_xy
    #execute_on = 'timestep_end'
  [../]
  [./matl_e13]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 2
    block = '2'
    variable = elastic_strain_xz
    #execute_on = 'timestep_end'
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    block = '2'
    variable = elastic_strain_yy
    #execute_on = 'timestep_end'
  [../]
  [./matl_e23]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 2
    block = '2'
    variable = elastic_strain_yz
    #execute_on = 'timestep_end'
  [../]
  [./matl_e33]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 2
    index_j = 2
    block = '2'
    variable = elastic_strain_zz
    #execute_on = 'timestep_end'
  [../]
[]

[Kernels]
  #Elastic problem
  [./TensorMechanics] #This is an action block
  [../]
  #Bulk energy density
  #[./bed_x]
  #  type = BulkEnergyDerivative
  #  variable = polar_x
  #  component = 0
  #[../]
  #[./bed_y]
  #  type = BulkEnergyDerivative
  #  variable = polar_y
  #  component = 1
  #[../]
  #[./bed_z]
  #  type = BulkEnergyDerivative
  #  variable = polar_z
  #  component = 2
  #[../]
  ##Wall energy penalty
  #[./walled_x]
  #   type=WallEnergyDerivative
  #   variable = polar_x
  #   component = 0
  #[../]
  #[./walled_y]
  #   type=WallEnergyDerivative
  #   variable = polar_y
  #   component = 1
  #[../]
  #[./walled_z]
  #   type=WallEnergyDerivative
  #   variable = polar_z
  #   component = 2
  #[../]
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
  #[./ferroelectriccouplingp_xx]
  #   type = FerroelectricCouplingP
  #   variable=polar_x
  #   component = 0
  #   block = '2'
  #[../]
  #[./ferroelectriccouplingp_yy]
  #   type = FerroelectricCouplingP
  #   variable=polar_y
  #   component = 1
  #   block = '2'
  #[../]
  #[./ferroelectriccouplingp_zz]
  #   type = FerroelectricCouplingP
  #   variable=polar_z
  #   component = 2
  #   block = '2'
  #[../]
  ##Electrostatics
  #[./polar_x_electric_E]
  #   type=PolarElectricEStrong
  #   variable=potential_int
  #   block='2'
  #[../]
  [./E_int]
     type=Electrostatics
     variable=potential_int
     block='1 3'
  [../]
  #[./FE_E_int]
  #   type=Electrostatics
  #   variable=potential_int
  #   block='2'
  #[../]
  #[./E_ext]
  #   type=Electrostatics
  #   variable=potential_ext
  #   block='1'
  #[../]
  #[./FE_E_ext]
  #   type=Electrostatics
  #   variable=potential_ext
  #   block='2'
  #[../]
  #[./polar_electric_px]
  #   type=PolarElectricPStrong
  #   variable = polar_x
  #   component = 0
  #[../]
  #[./polar_electric_py]
  #   type=PolarElectricPStrong
  #   variable = polar_y
  #   component = 1
  #[../]
  #[./polar_electric_pz]
  #   type=PolarElectricPStrong
  #   variable = polar_z
  #   component = 2
  #[../]
  ##Time dependence
  #[./polar_x_time]
  #   type=TimeDerivativeScaled #CoupledTimeDerivative?
  #   variable=polar_x
  #   time_scale = 1.0e2
  #[../]
  #[./polar_y_time]
  #   type=TimeDerivativeScaled
  #   variable=polar_y
  #   time_scale = 1.0e2
  #[../]
  #[./polar_z_time]
  #   type=TimeDerivativeScaled
  #   variable = polar_z
  #   time_scale = 1.0e2
  #[../]
  [./disp_x_time]
     type=TimeDerivativeScaled
     variable = disp_x
     time_scale = 0.5e-9
  [../]
  [./disp_y_time]
     type=TimeDerivativeScaled
     variable = disp_y
     time_scale = 0.5e-9
  [../]
  [./disp_z_time]
     type=TimeDerivativeScaled
     variable = disp_z
     time_scale = 0.5e-9
  [../]
[]

[BCs]
   [./potential_int_upz]
     type = DirichletBC
     variable = potential_int
     boundary = '1'
     value = 0.0
   [../]
   [./potential_int_downz]
     type = DirichletBC
     variable = potential_int
     boundary = '2'
     value = 0.0
   [../]
   #Top and bottom {3, 4}:
   [./disp_x_slab3]
     type = DirichletBC
     variable = disp_x
     boundary = '3'
     value = 0.0
   [../]
   [./disp_y_slab3]
     type = DirichletBC
     variable = disp_y
     boundary = '3'
     value = 0.0
   [../]
   [./disp_z_slab3]
     type = DirichletBC
     variable = disp_z
     boundary = '3'
     value = 0.0
   [../]
   [./disp_x_slab4]
     type = DirichletBC
     variable = disp_x
     boundary = '4'
     value = 0.0
   [../]
   [./disp_y_slab4]
     type = DirichletBC
     variable = disp_y
     boundary = '4'
     value = 0.0
   [../]
   [./disp_z_slab4]
     type = DirichletBC
     variable = disp_z
     boundary = '4'
     value = 0.0
   [../]

  #Sides {5, 6, 7, 8}:
  [./disp_x_slab5]
    type = DirichletBC
    variable = disp_x
    boundary = '5'
    value = 0.0
  [../]
  [./disp_y_slab5]
    type = DirichletBC
    variable = disp_y
    boundary = '5'
    value = 0.0
  [../]
  [./disp_z_slab5]
    type = DirichletBC
    variable = disp_z
    boundary = '5'
    value = 0.0
  [../]
  [./disp_x_slab6]
     type = DirichletBC
     variable = disp_x
     boundary = '6'
     value = 0.0
  [../]
  [./disp_y_slab6]
     type = DirichletBC
     variable = disp_y
     boundary = '6'
     value = 0.0
  [../]
  [./disp_z_slab6]
     type = DirichletBC
     variable = disp_z
     boundary = '6'
     value = 0.0
  [../]
  [./disp_x_slab7]
     type = DirichletBC
     variable = disp_x
     boundary = '7'
     value = 0.0
  [../]
  [./disp_y_slab7]
     type = DirichletBC
     variable = disp_y
     boundary = '7'
     value = 0.0
  [../]
  [./disp_z_slab7]
     type = DirichletBC
     variable = disp_z
     boundary = '7'
     value = 0.0
  [../]
  [./disp_x_slab8]
     type = DirichletBC
     variable = disp_x
     boundary = '8'
     value = 0.0
  [../]
  [./disp_y_slab8]
     type = DirichletBC
     variable = disp_y
     boundary = '8'
     value = 0.0
  [../]
  [./disp_z_slab8]
     type = DirichletBC
     variable = disp_z
     boundary = '8'
     value = 0.0
  [../]
  #[./potential_ext_upz]
  #  type = DirichletBC
  #  variable = potential_ext
  #  boundary = '1'
  #  value = 0.0
  #[../]
  #[./potential_ext_downz]
  #  type = DirichletBC
  #  variable = potential_ext
  #  boundary = '2'
  #  value = 0.0
  #[../]
[]

[Materials]
  #[./slab_ferroelectric]
  #  type=LinearFerroelectricMaterial
  #  block = '2'
  #  disp_x = disp_x
  #  disp_y = disp_y
  #  disp_z = disp_z
  #  #in GPA. from N. Pandech et al. Ceramic. Internat., times 1e9 to convert to N/m^2
  #  # C11 C12 C13 C22 C23 C33 C44 C55 C66
  #  C_ijkl = '3.80e-7 1.5e-7 1.50e-7 3.80e-7 1.50e-7 3.80e-7 1.1e-7 1.1e-7 1.1e-7'
  #  #in m^4/C^2. from http://arxiv.org/pdf/1205.5640.pdf
  #  # Q11 Q12 Q13 Q22 Q23 Q33 Q44 Q55 Q66
  #  Q_mnkl = '8.9e-2 -2.6e-2 -2.6e-2 8.9e-2 -2.6e-2 8.9e-2 3.4e-2 3.4e-2 3.4e-2'
  #  euler_angle_1 = 0.0 #currently will only rotate C_ijkl
  #  euler_angle_2 = 0.0
  #  euler_angle_3 = 0.0
  #[../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = '2'
    C_ijkl = '3.80e-7 1.5e-7 1.50e-7 3.80e-7 1.50e-7 3.80e-7 1.1e-7 1.1e-7 1.1e-7'
    fill_method = symmetric9
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = '2'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '2'
  [../]
  [./vacuum]
    type=GenericConstantMaterial
    block = '1 3'
  [../]
  ##This seems to be what we want for a simple epitaxial test
  ## (note that most epitaxial conditions are a strain gradient from the interface)
  [./eigen_strain_xx_yy] #Use for stress-free strain (ie epitaxial)
    type = ComputeEigenstrain
    block = '2'
    prefactor = -0.02
    # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
    eigen_base = '1 0 0 0 1 0 0 0 0'
  [../]
[]

[Postprocessors]
  #[./bulk_energy]
  # type = BulkEnergy
  # execute_on = 'timestep_end'
  ## initial_from_file_var = bulk_energy
  # block = '2'
  #[../]
  #[./wall_energy]
  # type = WallEnergy
  # execute_on = 'timestep_end'
  ## initial_from_file_var = wall_energy
  # block = '2'
  #[../]
  [./elastic_energy]
   type = ElasticEnergy
  # execute_on = 'timestep_end'
  # initial_from_file_var = elastic_energy
   block = '2'
  [../]
  #[./electrostatic_energy]
  # type = ElectrostaticEnergy
  # execute_on = 'timestep_end'
  ## initial_from_file_var = electrostatic_energy
  #[../]
  #[./total_energy]
  # type = TotalEnergy
  # bulk_energy = bulk_energy
  # wall_energy = wall_energy
  # elastic_energy = elastic_energy
  # electrostatic_energy = electrostatic_energy
  # execute_on = 'timestep_end'
  ## initial_from_file_var = total_energy
  #[../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -options_left'
    petsc_options_iname = '-ksp_gmres_restart -snes_rtol -ksp_rtol -pc_type -pc_factor_zeropivot -pc_side'
    petsc_options_value = '    675              1e-8      1e-14      hypre        1e-50              left'
  [../]
[]

[Executioner]
  type = Steady
  type = Transient
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 6
    optimal_iterations = 10
    growth_factor = 1.0001
    cutback_factor =  0.9999
  [../]
  solve_type = 'NEWTON'       #"PJNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  #dtmin = 1.0e-10
  dtmax = 6
  num_steps = 6
  #splitting = 'ferretsplit'
[]

#[Splits]
#  [./ferretsplit]
#    type = Split
#    splitting = 'ferroelectric elastic' #split to two subproblems
#    splitting_type = schur #schur split somewhat bugged right now (only serial)
#    schur_type = full
#    schur_pre = A11
#  [../]
#  [./ferroelectric]
#    vars = 'polar_x polar_y polar_z potential_int potential_ext'
#    petsc_options='-dm_view ' #'-ksp_monitor -inner_ksp_monitor'
#    petsc_options_iname=' -ksp_type   -ksp_gmres_restart  -ksp_rtol -inner_pc_type -pc_type   -pc_asm_overlap -inner_pc_factor_zeropivot  -pc_factor_zeropivot -pc_hypre_type -inner_pc_hypre_type'
#    petsc_options_value='    gmres            350              1e-10   hypre        hypre        2                         1e-50                1e-50  boomeramg boomeramg'
#  [../]
#  [./elastic]
#    vars = 'disp_x disp_y disp_z'
#    petsc_options='-dm_view'  # ''-ksp_monitor'
#    petsc_options_iname=' -ksp_type  -ksp_gmres_restart -inner_pc_type -ksp_rtol -pc_type   -pc_asm_overlap -inner_pc_factor_zeropivot -pc_factor_zeropivot -pc_hypre_type -inner_pc_hypre_type'
#    petsc_options_value = ' gmres       350                hypre         1e-10       hypre           2           1e-50                        1e-50 boomeramg boomeramg'
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
    file_base = out_PbTiO3_60nm_T650_Steady_noFE
    output_initial = true
    elemental_as_nodal = true
    interval = 1
  [../]
  #[./debug]
  #  type = VariableResidualNormsDebugOutput
  #[../]
[]
