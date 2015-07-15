[Mesh]
  file = slab_exodus_coarse_40.e  #if smaller mesh desired, use slab_exodus_coarse_150_cheap.e in /problems/coupled_system
  uniform_refine = 0
[]

[GlobalParams]
   #Can use a unit system where the differences between the coefficients
   # in this table are minimized: [nm], [fC], [hN]

   #length scale
   len_scale = 1.0
   #BulkEnergy coefficients
   alpha1 = -2.8576e-6 # 3.8(T-785)*10^5 C^{-2}nm^2 (T = 0 K)
   alpha11 = -7.3e-1
   alpha111 = 2.6e6
   alpha12 = 7.5e-1
   alpha112 = 6.1e6
   alpha123 = -3.7e7
   #WallEnergy coefficients
   G110 = 6.0e-7
   G11/G110 = 0.6
   G12/G110 = 0.0
   G44/G110 = 0.3
   G44P/G110 = 0.3
   #Electrostatics
   permittivity = 885.4187
   polar_x = polar_x
   polar_y = polar_y
   polar_z = polar_z
   potential_int = potential_int
   potential_ext = potential_ext
   #elastic variables
   disp_x = disp_x
   disp_y = disp_y
   disp_z = disp_z
   use_displaced_mesh = false
[]

[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block='2'
    [./InitialCondition]
      type = RandomIC
      min = 0.0001e-2
      max = 0.0003e-2
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block='2'
    [./InitialCondition]
      type = RandomIC
      min = 0.0001e-2
      max = 0.0003e-2
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block='2'
    #scaling = 1e5
    [./InitialCondition]
      type = RandomIC
      min = 0.0001e-3
      max = 0.0003e-3
    [../]
  [../]
  [./potential_int]
    order=FIRST
    family = LAGRANGE
  [../]
  [./potential_ext]
    order=FIRST
    family = LAGRANGE
  [../]
  #
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    block = '2'
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
    #[./InitialCondition]
    #  type = RandomIC
    #  min = 0.1e-1
    #  max = 1.5e-1
    #[../]
  [../]
[]

[AuxVariables]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./matl_s11]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    block = '2'
    variable = stress_xx
  [../]
  [./matl_s12]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    block = '2'
    variable = stress_xy
  [../]
  [./matl_s13]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
    block = '2'
    variable = stress_xz
  [../]
  [./matl_s22]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    block = '2'
    variable = stress_yy
  [../]
  [./matl_s23]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    block = '2'
    variable = stress_yz
  [../]
  [./matl_s33]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    block = '2'
    variable = stress_zz
  [../]
  [./matl_e11]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    block = '2'
    variable = strain_xx
  [../]
  [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    block = '2'
    variable = strain_xy
  [../]
  [./matl_e13]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 2
    block = '2'
    variable = strain_xz
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    block = '2'
    variable = strain_yy
  [../]
  [./matl_e23]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 2
    block = '2'
    variable = strain_yz
  [../]
  [./matl_e33]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 2
    index_j = 2
    block = '2'
    variable = strain_zz
  [../]
[]

[Kernels]
  #Elastic problem
  [./stressdiv_0]
    type = StressDivergenceTensorsScaled
    variable = disp_x
    component = 0
    block = '2'
  [../]
  [./stressdiv_1]
    type = StressDivergenceTensorsScaled
    variable = disp_y
    component = 1
    block = '2'
  [../]
  [./stressdiv_2]
    type = StressDivergenceTensorsScaled
    variable = disp_z
    component = 2
    block = '2'
  [../]
  ###Bulk energy density
  [./bed_x]
    type = BulkEnergyDerivative
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = BulkEnergyDerivative
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = BulkEnergyDerivative
    variable = polar_z
    component = 2
  [../]
  #
  ###Wall energy penalty
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
     block = '2'
  [../]
  [./ferroelectriccouplingp_yy]
     type = FerroelectricCouplingP
     variable=polar_y
     component = 1
     block = '2'
  [../]
  [./ferroelectriccouplingp_zz]
     type = FerroelectricCouplingP
     variable=polar_z
     component = 2
     block = '2'
  [../]
  #Electrostatics
  [./polar_x_electric_E]
     type=PolarElectricEStrong
     variable=potential_int
     block='2'
  [../]
  [./E_int]
     type=Electrostatics
     variable=potential_int
     block='1'
  [../]
  [./FE_E_int]
     type=Electrostatics
     variable=potential_int
     block='2'
  [../]
  [./E_ext]
     type=Electrostatics
     variable=potential_ext
     block='1'
  [../]
  [./FE_E_ext]
     type=Electrostatics
     variable=potential_ext
     block='2'
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
  #Time dependence
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
  #[./disp_x_time]
  #   type=TimeDerivativeScaled
  #   variable = disp_x
  #   time_scale = 1.0
  #[../]
  #[./disp_y_time]
  #   type=TimeDerivativeScaled
  #   variable = disp_y
  #   time_scale = 1.0
  #[../]
  #[./disp_z_time]
  #   type=TimeDerivativeScaled
  #   variable = disp_z
  #   time_scale = 1.0
  #[../]
[]

[Materials]
  [./slab_ferroelectric]
    type=LinearFerroelectricMaterial
    block = '2'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    #in GPA. from N. Pandech et al. Ceramic. Internat., times 1e9 to convert to N/m^2
    # C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '3.80e-17 1.5e-17 1.50e-17 3.80e-17 1.50e-17 3.80e-17 1.1e-17 1.1e-17 1.1e-17'
    #in m^4/C^2. from http://arxiv.org/pdf/1205.5640.pdf
    # Q11 Q12 Q13 Q22 Q23 Q33 Q44 Q55 Q66
    Q_mnkl = '8.9e4 -2.6e4 -2.6e4 8.9e4 -2.6e4 8.9e4 3.4e4 3.4e4 3.4e4'
    euler_angle_1 = 0.0 #currently will only rotate C_ijkl
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = '2'
    C_ijkl = '3.80e-17 1.5e-17 1.50e-17 3.80e-17 1.50e-17 3.80e-17 1.1e-17 1.1e-17 1.1e-17'
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
    block = '1'
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

   [./disp_x_slab3]
     type = PresetBC
     variable = disp_x
     boundary = '3'
     value = 0.0
   [../]
   [./disp_y_slab3]
     type = PresetBC
     variable = disp_y
     boundary = '3'
     value = 0.0
   [../]
   [./disp_z_slab3]
     type = PresetBC
     variable = disp_z
     boundary = '3'
     value = 0.0
   [../]
  # #
   [./disp_y_slab5]
     type = PresetBC
     variable = disp_y
     boundary = '5'
     value = -0.1 #probably ~2% strain
   [../]
   [./disp_y_slab7]
     type = PresetBC
     variable = disp_y
     boundary = '7'
     value = 0.1
   [../]
   [./potential_ext_upz]
    type = DirichletBC
    variable = potential_ext
    boundary = '1'
    value = 0.0 #this is near-zero
   [../]
   [./potential_ext_downz]
    type = DirichletBC
    variable = potential_ext
    boundary = '2'
    value = 0.0
   [../]
[]

[Postprocessors]
  [./bulk_energy]
   type = BulkEnergy
   block = '2'
  [../]
  [./wall_energy]
   type = WallEnergy
   block = '2'
  [../]
  [./elastic_energy]
   type = ElasticEnergy
   block = '2'
  [../]
  [./electrostatic_energy]
   type = ElectrostaticEnergy
  [../]
  [./total_energy]
   type = TotalEnergy
   bulk_energy = bulk_energy
   wall_energy = wall_energy
   elastic_energy = elastic_energy
   electrostatic_energy = electrostatic_energy
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-info -snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -options_left'
    petsc_options_iname = '-ksp_gmres_restart -snes_rtol -ksp_rtol -pc_type  -pc_asm_overlap -sub_pc_type  -sub_pc_factor_zeropivot -pc_factor_zeropivot '
    petsc_options_value = '    121              1e-8      1e-10      asm          2             lu             1e-50                    1e-50      '
  [../]
[]

[Executioner]
  type=Transient
  #nl_max_its = 5
  #nl_abs_tol = 6.90e-22
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0e1
    optimal_iterations = 10
    growth_factor = 1.01
    cutback_factor =  0.85
  [../]
  solve_type = 'NEWTON'       #"PJNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1.0e-10
  dtmax = 1.0e1
  num_steps = 1500000
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
#    petsc_options_iname=' -ksp_type   -ksp_gmres_restart  -ksp_rtol -inner_pc_type -pc_type   -pc_asm_overlap -inner_pc_factor_zeropivot  -pc_factor_zeropivot'
#    petsc_options_value='    gmres            350              1e-14   asm        asm        5                         1e-50                1e-50  '
#  [../]
#  [./elastic]
#    vars = 'disp_x disp_y disp_z'
#    petsc_options='-dm_view'  # ''-ksp_monitor'
#    petsc_options_iname=' -ksp_type  -ksp_gmres_restart -inner_pc_type -ksp_rtol -pc_type   -pc_asm_overlap -inner_pc_factor_zeropivot -pc_factor_zeropivot'
#    petsc_options_value = ' gmres       350                 asm         1e-14       asm           5           1e-50                        1e-50'
#  [../]
#[]

[./Adaptivity]
    [./Indicators]
      [./indicator_x]
        type = GradientJumpIndicator
        variable = polar_x
      [../]
      [./indicator_y]
        type = GradientJumpIndicator
        variable = polar_y
      [../]
      [./indicator_z]
        type = GradientJumpIndicator
        variable = polar_z
      [../]
    [../]
    #[./Markers]
    #  [./marker_x]
    #    type = ErrorFractionMarker
    #    indicator = indicator_x
    #    coarsen = 0.01
    #    refine = 0.01
    #  [../]
    #  [./marker_y]
    #    type = ErrorFractionMarker
    #    indicator = indicator_y
    #    coarsen = 0.01
    #    refine = 0.01
    #  [../]
    #[../]
[../]


[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_PBP_test_coupled
    output_initial = true
    elemental_as_nodal = false
    interval = 1
  [../]
  #[./debug]
  #  type = VariableResidualNormsDebugOutput
  #[../]
[]
