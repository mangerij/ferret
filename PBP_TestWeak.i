[Mesh]
  file = slab_exodus_coarse_40.e  #if smaller mesh desired, use slab_exodus_coarse_150_cheap.e in /problems/coupled_system
  uniform_refine = 0
[]

#[Problem]
#  dimNullSpace     = 6
#[]
#
#[UserObjects]
#  [./RigidModes3DNullSpace]
#     type=RigidBodyModes3D
#     variable = disp_x
#     subspace_name = NullSpace
#     subspace_indices = '0 1 2 3 4 5 '
#     modes = 'trans_x trans_y trans_z rot_x rot_y rot_z'
#     disp_x = disp_x
#     disp_y = disp_y
#     disp_z = disp_z
#  [../]
#  [./kill]
#    type = Terminator
#    expression = '_pps_percent <= 4.0e-7'
#  [../]
#[]

[GlobalParams]
   #length scale
   len_scale = 1e-9
   #BulkEnergy coefficients
   alpha1 = -1.8202e8 # 3.8(T-479)*10^5 C^{-2}m^2 (T = 0 K)
   alpha11 = -7.3e7
   alpha111 = 2.6e8
   alpha12 = 7.5e8
   alpha112 = 6.1e8
   alpha123 = -3.7e9
   #WallEnergy coefficients
   G110 = 0.6e-10
   G11/G110 = 0.6
   G12/G110 = 0.0
   G44/G110 = 0.3
   G44P/G110 = 0.3
   #Electrostatics
   permittivity = 8.85e-12
   polar_x = polar_x
   polar_y = polar_y
   polar_z = polar_z
   potential_int=potential_int
   potential_ext=potential_ext
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
    #scaling = 1e13
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block='2'
    #scaling = 1e13
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block='2'
    #scaling = 1e13
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

[AuxKernels]
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
  #[./TensorMechanics]
  #  disp_x = disp_x
  #  disp_y = disp_y
  #  disp_z = disp_z
  #[../]
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

  #Wall energy penalty
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
  #Polarization-strain coupling
  [./ferroelectriccouplingu_x]
     type = FerroelectricCouplingU
     variable = disp_x
     component = 0
     block = '2'
  [../]
  [./ferroelectriccouplingu_y]
     type = FerroelectricCouplingU
     variable = disp_y
     component = 1
     block = '2'
  [../]
  [./ferroelectriccouplingu_z]
     type = FerroelectricCouplingU
     variable=disp_z
     component = 2
     block = '2'
  [../]
  #
  [./ferroelectriccouplingp_xx]
     type = FerroelectricCouplingP
     variable=polar_x
     component = 0
     block = '2'
  [../]
  #
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
  [./polar_electric_E]
     type=PolarElectricEStrong
     variable=potential_int
    # implict = false
     block='2'
  [../]
  [./E_int]
     type=Electrostatics
     variable=potential_int
    # implicit = false
     block='1'
  [../]
  [./FE_E_int]
     type=Electrostatics
     variable=potential_int
    # implicit = false
     block='2'
  [../]
  [./E_ext]
     type=Electrostatics
     variable=potential_ext
    # implicit = false
     block='1'
  [../]
  [./FE_E_ext]
     type=Electrostatics
     variable=potential_ext
    # implicit = false
     block='2'
  [../]
  [./polar_electric_px]
     type=PolarElectricPStrong
     variable = polar_x
    # implict = false
     component = 0
  [../]
  [./polar_electric_py]
     type=PolarElectricPStrong
     variable=polar_y
    # implict = false
     component = 1
  [../]
  [./polar_electric_pz]
     type=PolarElectricPStrong
     variable=polar_z
    # implict = false
     component = 2
  [../]

  #Time dependence
  [./polar_x_time]
     type=TimeDerivativeScaled
     variable=polar_x
     time_scale = 1.0e-29
  [../]
  [./polar_y_time]
     type=TimeDerivativeScaled
     variable=polar_y
     time_scale = 1.0e-29
  [../]
  [./polar_z_time]
     type=TimeDerivativeScaled
     variable=polar_z
     time_scale = 1.0e-29
  [../]
  [./disp_x_time]
     type=TimeDerivativeScaled
     variable=disp_x
     time_scale = 1.0e-29
  [../]
  [./disp_y_time]
     type=TimeDerivativeScaled
     variable=disp_y
     time_scale = 1.0e-29
  [../]
  [./disp_z_time]
     type=TimeDerivativeScaled
     variable=disp_z
     time_scale = 1.0e-29
  [../]
[]

[ICs]
  [./polar_x_randic]
     type=RandomIC
     variable=polar_x
     min = 0.00006
     max = 0.000075
  [../]
  [./polar_y_randic]
     type=RandomIC
     variable=polar_y
     min = 0.00006
     max = 0.000075
  [../]
  [./polar_z_randic]
     type=RandomIC
     variable=polar_z
     min = 0.00001
     max = 0.000015
  [../]
  #[./disp_x_randic]
  #   type=RandomIC
  #   variable=disp_x
  #   min = 6e-10
  #   max = 7.5e-10
  #[../]
  #[./disp_y_randic]
  #   type=RandomIC
  #   variable=disp_y
  #   min = 6e-10
  #   max = 7.5e-10
  #[../]
  #[./disp_z_randic]
  #   type=RandomIC
  #   variable=disp_z
  #   min = 1e-10
  #   max = 1.5e-10
  #[../]
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

   #Displacements [5 6 7 8] [x y -x -y] boundaries, [3 4] [top bottom]
  # [./disp_x_slab3]
  #   type = PresetBC
  #   variable = disp_x
  #   boundary = '3'
  #   value = 0.0
  # [../]
  # [./disp_y_slab3]
  #   type = PresetBC
  #   variable = disp_y
  #   boundary = '3'
  #   value = 0.0
  # [../]
  # [./disp_z_slab3]
  #   type = PresetBC
  #   variable = disp_z
  #   boundary = '3'
  #   value = 0.0
  # [../]
   #
  # [./disp_x_slab4]
  #   type = PresetBC
  #   variable = disp_x
  #   boundary = '4'
  #   value = 0.0
  # [../]
  # [./disp_y_slab4]
  #   type = PresetBC
  #   variable = disp_y
  #   boundary = '4'
  #   value = 0.0
  # [../]
  # [./disp_z_slab4]
  #   type = PresetBC
  #   variable = disp_z
  #   boundary = '4'
  #   value = 0.0
  # [../]

  # [./disp_x_slab5]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = '5'
  #  # value = 1.0e-1
  #   value = 0.5
  # [../]
   [./disp_y_slab5]
     type = DirichletBC
     variable = disp_y
     boundary = '5'
     value = 0.5
   [../]
  # [./disp_z_slab5]
  #   type = PresetBC
  #   variable = disp_z
  #   boundary = '5'
  #  # value = 1.0e-1
  #   value = 0.0
  # [../]

  # [./disp_x_slab6]
  #   type = PresetBC
  #   variable = disp_x
  #   boundary = '6'
  #  # value = 1.0e-1
  #  value = 0.0
  # [../]
  # [./disp_y_slab6]
  #   type = PresetBC
  #   variable = disp_y
  #   boundary = '6'
  #  # value = 1.0e-1
  #  value = 0.0
  # [../]
  # [./disp_z_slab6]
  #   type = PresetBC
  #   variable = disp_z
  #   boundary = '6'
  #  # value = 1.0e-1
  #  value = 0.0
  # [../]

  # [./disp_x_slab7]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = '7'
  #  # value = -1.0e-1
  #   value = -0.5
  # [../]
   [./disp_y_slab7]
     type = DirichletBC
     variable = disp_y
     boundary = '7'
     value = -0.5
   [../]
  # [./disp_z_slab7]
  #   type = PresetBC
  #   variable = disp_z
  #   boundary = '7'
  #  # value = -1.0e-1
  #   value = 0.0
  # [../]

 #  [./disp_x_slab8]
 #    type = PresetBC
 #    variable = disp_x
 #    boundary = '8'
 #   # value = -1.0e-1
 #   value = 0.0
 #  [../]
 #  [./disp_y_slab8]
 #    type = PresetBC
 #    variable = disp_y
 #    boundary = '8'
 #   # value = -1.0e-1
 #   value = 0.0
 #  [../]
 #  [./disp_z_slab8]
 #    type = PresetBC
 #    variable = disp_z
 #    boundary = '8'
 #   # value = -1.0e-1
 #   value = 0.0
 #[../]
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
   block = '2'
  [../]
  [./wall_energy]
   type=WallEnergy
   block = '2'
  [../]
  [./electrostatic_energy]
   type=ElectrostaticEnergy
  [../]
  [./elastic_energy]
   type = ElasticEnergy
   block = '2'
  [../]
  #[./coupled_energy]
  #  type = CoupledEnergy
  #  block = '2'
  #[../]
  [./total_energy]
   type=TotalEnergy
   bulk_energy = bulk_energy
   wall_energy = wall_energy
   elastic_energy = elastic_energy
   electrostatic_energy=electrostatic_energy
  [../]
  #[./_pps_percent]
  #  type = PercentChangePostprocessor
  #  postprocessor = total_energy
  #[../]
  #[./R(i)]
  #  type = Residual
  #[../]
[]

[Executioner]
  type=Transient
  nl_max_its = 350
  #[./TimeStepper]
  #  type = IterationAdaptiveDT
  #  dt = 1.15e-15
  #  optimal_iterations = 10
  #  growth_factor = 1.01
  #  cutback_factor =  0.99
  #[../]
  solve_type = 'newton'
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1.0e-25
  dtmax = 1.81e15
  num_steps = 15000
  dt = 1.15e-15
  #splitting = 'ferretsplit'
  petsc_options ='-snes_linesearch_monitor -options_left'
  petsc_options_iname = '-snes_rtol -ksp_rtol -pc_type  -sub_pc_type -pc_asm_overlap -sub_pc_factor_zeropivot -pc_factor_zeropivot'
  petsc_options_value = '  1e-8       1e-14      asm         lu            2               1e-50                    1e-50'
[]

#
#[Splits]
#  [./ferretsplit]
#    type = Split
#    splitting = 'ferroelectric elastic' #split to two subproblems
#    splitting_type = schur #schur split somewhat bugged right now (only serial)
#    schur_type = full
#    schur_pre = A11
#    petsc_options ='-snes_linesearch_monitor '
#  [../]
#
#  [./ferroelectric]
#    vars = 'polar_x polar_y polar_z potential_int potential_ext'
#    petsc_options='-dm_view ' #'-ksp_monitor -inner_ksp_monitor'
#    petsc_options_iname='-ksp_type   -ksp_gmres_restart  -ksp_rtol -inner_pc_type -pc_type   -pc_asm_overlap -inner_pc_factor_zeropivot  -pc_factor_zeropivot'
#    petsc_options_value=' gmres            350              1e-14     hypre         hypre         5         1e-50                1e-50  '
#  [../]
#  [./elastic]
#    vars = 'disp_x disp_y disp_z'
#    petsc_options='-dm_view'  # ''-ksp_monitor'
#    petsc_options_iname='-ksp_type  -ksp_gmres_restart -inner_pc_type -ksp_rtol -pc_type   -pc_asm_overlap -inner_pc_factor_zeropivot -pc_factor_zeropivot'
#    petsc_options_value = 'gmres       350                 hypre         1e-14     hypre          5           1e-50                        1e-50'
#  [../]
#[]


[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_PBP_test
    output_initial = true
    elemental_as_nodal = false
    interval = 1
  [../]
  #[./debug]
  #  type = VariableResidualNormsDebugOutput
  #[../]
[]
