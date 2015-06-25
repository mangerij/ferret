[Mesh]
  file = slab_exodus_coarse_40.e
  uniform_refine=0
[]

[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block='2'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block='2'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block='2'
  [../]
  [./potential_int]
    order=FIRST
    family = LAGRANGE
  [../]
  [./potential_ext]
    order=FIRST
    family = LAGRANGE
  [../]

  [./disp_x]
    order = FIRST
    family = LAGRANGE
    block = '2'
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    block = '2'
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    block = '2'
  [../]

[]

[AuxVariables]
    [./rho_b]
      order = CONSTANT
      family = MONOMIAL
      block = '2'
    [../]

   [./stress_xx]
     order = CONSTANT
     family = MONOMIAL
     block = '2'
   [../]
   [./stress_yy]
     order = CONSTANT
     family = MONOMIAL
     block = '2'
   [../]
   [./stress_zz]
     order = CONSTANT
     family = MONOMIAL
     block = '2'
   [../]
   [./stress_xy]
     order = CONSTANT
     family = MONOMIAL
     block = '2'
   [../]
   [./stress_yz]
     order = CONSTANT
     family = MONOMIAL
     block = '2'
   [../]
   [./stress_zx]
     order = CONSTANT
     family = MONOMIAL
     block = '2'
   [../]
   [./strain_xx]
     order = CONSTANT
     family = MONOMIAL
     block = '2'
   [../]
   [./strain_yy]
     order = CONSTANT
     family = MONOMIAL
     block = '2'
   [../]
   [./strain_zz]
     order = CONSTANT
     family = MONOMIAL
     block = '2'
   [../]
   [./strain_xy]
     order = CONSTANT
     family = MONOMIAL
     block = '2'
   [../]
   [./strain_yz]
     order = CONSTANT
     family = MONOMIAL
     block = '2'
   [../]
   [./strain_zx]
     order = CONSTANT
     family = MONOMIAL
     block = '2'
   [../]
[]

[GlobalParams]
   len_scale=1e-9
   #BulkEnergy coefficients
   alpha1 = -1.8202e8 # 3.8(T-479)*10^5 C^{-2}m^2 (T=0 K)
   alpha11 = -7.3e7
   alpha111 = 2.6e8
   alpha12 = 7.5e8
   alpha112 = 6.1e8
   alpha123 = -3.7e9
   #WallEnergy coefficients
   G110 = 0.6e-10
   G11/G110=0.6
   G12/G110=0.0
   G44/G110=0.3
   G44P/G110=0.3
   #Electrostatics
   permittivity=8.85e-12
   polar_x = polar_x
   polar_y = polar_y
   polar_z = polar_z
   potential_int=potential_int
   potential_ext=potential_ext
   #
   disp_x = disp_x
   disp_y = disp_y
   disp_z = disp_z
   time_scale = 1.0e-31
   # time_scale ostensibly defines dt = 1e-16. Note if time_scale = 1.0 then dt=1e17 to see domain relaxation
   # Olle suggests that the dynamics should occur on the 10 ps time scale. This parameter allows the time scale to be
   # user-input and most likely depends on material and problem in question, maybe even temperature.
   # can just scale to experimentally relevant scale
[]

[AuxKernels]
  #[./auxk_electrostatic_energy_density_e]
  #  type = ElectrostaticEnergyDensityE
  #  variable =auxv_es_energy_density_e
  #  potential=potential_int
  #[../]
  #[./auxk_electrostatic_energy_density]
  #  type = ElectrostaticEnergyDensity
  #  variable =auxv_es_energy_density
  #[../]
  #[./auxk_bulk_energy_density]
  #  type = BulkEnergyDensity
  #  variable =auxv_bulk_energy_density
  #[../]
  [./bound_charge]
    type = BoundCharge
    variable = rho_b
    block = '2'
  [../]
  #[./Ez_fieldaux]
  #  type = Ez_fieldAux
  #  variable = Ez
  #[../]
  #[./Ex_fieldaux]
  #  type = Ex_fieldAux
  #  variable = Ex
  #[../]
  #[./Ey_fieldaux]
  #  type = Ey_fieldAux
  #  variable = Ey
  #[../]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
    block = '2'
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    block = '2'
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
    block = '2'
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
    block = '2'
  [../]
  [./stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 1
    index_j = 2
    block = '2'
  [../]
  [./stress_zx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zx
    index_i = 2
    index_j = 0
    block = '2'
  [../]
  [./strain_xx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_xx
    index_i = 0
    index_j = 0
    block = '2'
  [../]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
    block = '2'
  [../]
  [./strain_zz]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zz
    index_i = 2
    index_j = 2
    block = '2'
  [../]
  [./strain_xy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_xy
    index_i = 0
    index_j = 1
    block = '2'
  [../]
  [./strain_yz]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_yz
    index_i = 1
    index_j = 2
    block = '2'
  [../]
  [./strain_zx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zx
    index_i = 2
    index_j = 0
    block = '2'
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


  #Bulk energy density
  [./bed_x]
    type = BulkEnergyDerivative
    variable = polar_x
    component=0
  [../]
  [./bed_y]
    type = BulkEnergyDerivative
    variable = polar_y
    component=1
  [../]
  [./bed_z]
    type = BulkEnergyDerivative
    variable = polar_z
    component=2
  [../]

  #Wall energy penalty
  [./walled_x]
     type=WallEnergyDerivative
     variable=polar_x
     component=0
  [../]
  [./walled_y]
     type=WallEnergyDerivative
     variable=polar_y
     component=1
  [../]
  [./walled_z]
     type=WallEnergyDerivative
     variable=polar_z
     component=2
  [../]
  #Ferroelectric-strain coupling
  [./ferroelectriccouplingu_x]
     type = FerroelectricCouplingU
     variable=polar_x
     component=0
     block = '2'
  [../]
  [./ferroelectriccouplingu_y]
     type = FerroelectricCouplingU
     variable=polar_y
     component=1
     block = '2'
  [../]
  [./ferroelectriccouplingu_z]
     type = FerroelectricCouplingU
     variable=polar_z
     component=2
     block = '2'
  [../]

  [./ferroelectriccouplingp_x]
     type = FerroelectricCouplingP
     variable=polar_x
     component=0
     block = '2'
  [../]
  [./ferroelectriccouplingp_y]
     type = FerroelectricCouplingP
     variable=polar_y
     component=1
     block = '2'
  [../]
  [./ferroelectriccouplingp_z]
     type = FerroelectricCouplingP
     variable=polar_z
     component=2
     block = '2'
  [../]

  #Electrostatics
  [./polar_electric_E]
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
     variable=polar_x
     component=0
  [../]
  [./polar_electric_py]
     type=PolarElectricPStrong
     variable=polar_y
     component=1
  [../]
  [./polar_electric_pz]
     type=PolarElectricPStrong
     variable=polar_z
     component=2
  [../]
  [./polar_x_time]
     type=TimeDerivativeScaled
     variable=polar_x
     time_scale = 1.0e-31
  [../]
  [./polar_y_time]
     type=TimeDerivativeScaled
     variable=polar_y
     time_scale = 1.0e-31
  [../]
  [./polar_z_time]
     type=TimeDerivativeScaled
     variable=polar_z
     time_scale = 1.0e-31
  [../]

  #[./disp_x_time]
  #   type=TimeDerivativeScaled
  #   variable=disp_x
  #   time_scale = 1e-35
  #   block = '2'
  #[../]
  #[./disp_y_time]
  #   type=TimeDerivativeScaled
  #   variable=disp_y
  #   time_scale = 1e-35
  #   block = '2'
  #[../]
  #[./disp_z_time]
  #   type=TimeDerivativeScaled
  #   variable=disp_z
  #   time_scale = 1e-35
  #   block = '2'
  #[../]
[]

[ICs]
  [./polar_x_randic]
     type=RandomIC
     variable=polar_x
     min = -0.25
     max = 0.3
  [../]
  [./polar_y_randic]
     type=RandomIC
     variable=polar_y
     min = -0.25
     max = 0.3
  [../]
  [./polar_z_randic]
     type=RandomIC
     variable=polar_z
     min = -0.25
     max = 0.3
  [../]
  [./disp_x_randic]
     type=RandomIC
     variable=disp_x
     min = -0.25
     max = 0.3
  [../]
  [./disp_y_randic]
     type=RandomIC
     variable=disp_y
     min = -0.25
     max = 0.3
  [../]
  [./disp_z_randic]
     type=RandomIC
     variable=disp_z
     min = -0.25
     max = 0.3
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


   [./disp_x_slab4]
     type = DirichletBC
     variable = disp_x
     boundary = '4'
     value = 0
   [../]
   [./disp_y_slab4]
     type = DirichletBC
     variable = disp_y
     boundary = '4'
     value = 0
   [../]
   [./disp_z_slab4]
     type = DirichletBC
     variable = disp_z
     boundary = '4'
     value = 0
   [../]
  # [./disp_y_slab7]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = '7'
  #   value = -1.0
  # [../]
   #
  # [./disp_x_slab5]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = '6'
  #   value = 1.0
  # [../]
  # [./disp_x_slab7]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = '8'
  #   value = -1.0
  # [../]

     # Applied field: for zero field use NeumannBC on the external potential = 0. A
     # Note that \nabla^2 \Phi_{ext} = 0 is satisfied if \Phi_{ext} = 0, ie Dirichlet and Neumann BC classes are equivalent
   [./potential_ext_upz]
    type = NeumannBC
    variable = potential_ext
    boundary = '1'
    #value = 5.0e-11 #1e-12 gives a field of ~1e8 V/m ~ 1000 kV/cm for size 1, what about size 100?
                    # seems that this doesn't matter here? dt needs to change for the size however...
                    # lowering this value, to about 1e-16 is the limit here. 5e-11 gives 5.6e9 V/m for 100 size
    value = 0.0
  [../]
  [./potential_ext_downz]
    type = NeumannBC
    variable = potential_ext
    boundary = '2'
    #value = -5.0e-11
    value = 0.0
  [../]
[]
#
[Materials]
  [./slab_ferroelectric]
    type=LinearFerroelectricMaterial
    block = '2'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    #in GPA. from N. Pandech et al. Ceramic. Internat., times 1e9 to convert to N/m^2
    # C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '380.0e9 150.0e9 150.0e9 380.0e9 150.0e9 380.0e9 110.0e9 110.0e9 110.0e9'
    #in m^4/C^2. from http://arxiv.org/pdf/1205.5640.pdf
    # Q11 Q12 Q13 Q22 Q23 Q33 Q44 Q55 Q66
    Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.034 0.034 0.034'
    euler_angle_1 = 0.0 #currently will only rotate C_ijkl
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]
  [./slab_elastic]
    type=LinearElasticMaterial
    block = '2'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    #in GPA. from N. Pandech et al. Ceramic. Internat., times 1e9 to convert to N/m^2
    # C11 C12 C13 C22 C23 C33 C44 C55 C66
    C_ijkl = '380.0e9 150.0e9 150.0e9 380.0e9 150.0e9 380.0e9 110.0e9 110.0e9 110.0e9'
    euler_angle_1 = 0.0 #currently will only rotate C_ijkl
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]

  [./vacuum]
    type=GenericConstantMaterial
    block = '1'
  [../]
[]




[Postprocessors]
  [./bulk_energy]
   type=BulkEnergy
  [../]
  [./wall_energy]
   type=WallEnergy
  [../]
  [./electrostatic_energy]
   type=ElectrostaticEnergy
  [../]
  [./elastic_energy]
   type = ElasticEnergy
   block = '2'
  [../]
  [./total_energy]
   type=TotalEnergy
   bulk_energy = bulk_energy
   wall_energy = wall_energy
   elastic_energy = elastic_energy
   electrostatic_energy=electrostatic_energy
  [../]
 # [./_dt]
 #   type = TimestepSize
 # [../]
  [./_pps_percent]
    type = PercentChangePostprocessor
    postprocessor = total_energy
  [../]
  [./R(i)]
    type=Residual
  [../]
[]

#[UserObjects]
#  [./kill]
#    type = Terminator
#    expression = '_pps_percent <= 4.0e-7'
#  [../]
#[]

[Executioner]
  type=Transient

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.1e-15
    optimal_iterations = 3
    growth_factor = 1.001
    cutback_factor =  0.999
  [../]

  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"

  dtmin=1.0e-26
  dtmax=1.81e-11
  num_steps=1500
  petsc_options='-options_left -snes_converged_reason -ksp_converged_reason -snes_monitor -ksp_monitor_true_residual -snes_ls_monitor -snes_view'
  petsc_options_iname='-ksp_type  -snes_rtol -ksp_rtol -snes_type  -pc_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'gmres     1e-8      1e-6     newtonls        asm        ilu          2'
  splitting = 'ferretsplit'
[]

[Splits]
  [./ferretsplit]
    splitting = 'ferroelectric' #split to two subproblems
    splitting_type = schur
    schur_type = full
    schur_pre = A11
    petsc_options ='-pc_fieldsplit_detect_coupling'
  [../]

  [./ferroelectric]
    vars = 'polar_x polar_y polar_z potential_int potential_ext'
    petsc_options='-dm_view'
    petsc_options_iname='-ksp_type   -ksp_rtol -pc_type -pc_hypre_boomeramg_strong_threshold -snes_linesearch_type -pc_factor_zeropivot'
    petsc_options_value=' gmres        1e-8     hypre              0.5                          basic         1e-50  '
  [../]
  [./elastic]
    vars = 'disp_x disp_y disp_z'
    petsc_options='-dm_view'
    petsc_options_iname='-ksp_type   -ksp_rtol   -pc_type -sub_pc_type -pc_asm_overlap'
    petsc_options_value = 'gmres       1e-6        asm        lu          2'
[]



[Outputs]
  print_linear_residuals = true
  print_perf_log = true

  [./out]
    type = Exodus
    file_base = out_PbTiO3_T0_40nm_0field_G110-6e-1_coupling
    output_initial = true
    elemental_as_nodal = false
    interval = 1
  [../]
[]
