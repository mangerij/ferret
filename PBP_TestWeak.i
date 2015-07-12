[Mesh]
  file = slab_exodus_coarse_40.e  #if smaller mesh desired, use slab_exodus_coarse_150_cheap.e in /problems/coupled_system
  uniform_refine = 0
[]

[GlobalParams]
   #length scale
   len_scale = 1.0e-9
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
   permittivity = 8.854187e-12
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
      type = ConstantIC
      value = 6.5e-1
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block='2'
    [./InitialCondition]
      type = ConstantIC
      value = 6.5e-1
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block='2'
    [./InitialCondition]
      type = ConstantIC
      value = 6.5e-1
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

#[AuxVariables]
#  [./stress_xx]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./stress_yy]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#  [./stress_zz]
#    order = CONSTANT
#    family = MONOMIAL
#  [../]
#[]
#
#[AuxKernels]
#  [./matl_e11]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    index_i = 0
#    index_j = 0
#    variable = stress_xx
#  [../]
#  [./matl_e22]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    index_i = 1
#    index_j = 1
#    variable = stress_yy
#  [../]
#  [./matl_e33]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    index_i = 2
#    index_j = 2
#    variable = stress_zz
#  [../]
#[]

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
  ##Bulk energy density
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
     time_scale = 1.0e-32
  [../]
  [./disp_y_time]
     type=TimeDerivativeScaled
     variable=disp_y
     time_scale = 1.0e-32
  [../]
  [./disp_z_time]
     type=TimeDerivativeScaled
     variable=disp_z
     time_scale = 1.0e-32
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
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = '2'
    C_ijkl = '380.0e9 150.0e9 150.0e9 380.0e9 150.0e9 380.0e9 110.0e9 110.0e9 110.0e9'
    fill_method = symmetric9
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
     type = PresetBC
     variable = potential_int
     boundary = '2'
     value = 0.0
   [../]
   [./disp_y_slab5]
     type = DirichletBC
     variable = disp_y
     boundary = '5'
     value = -1.0
   [../]
   [./disp_y_slab7]
     type = DirichletBC
     variable = disp_y
     boundary = '7'
     value = 1.0
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
  [./total_energy]
   type=TotalEnergy
   bulk_energy = bulk_energy
   wall_energy = wall_energy
   elastic_energy = elastic_energy
   electrostatic_energy=electrostatic_energy
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-info -snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -options_left'
    petsc_options_iname = '-ksp_gmres_restart -snes_rtol -ksp_rtol -pc_type  -pc_asm_overlap -sub_pc_type  -sub_pc_factor_zeropivot -pc_factor_zeropivot -pc_hypre_type'
    petsc_options_value = '    101              1e-8      1e-13      hypre          2              lu             1e-50                    1e-50     boomeramg'
  [../]
[]

[Executioner]
  type=Transient
  nl_max_its = 350
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.61e-16
    optimal_iterations = 5
    growth_factor = 1.001
    cutback_factor =  0.999
  [../]
  solve_type = 'PJFNK'
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1.0e-38
  dtmax = 1.61e-16
  num_steps = 15000
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
[]
