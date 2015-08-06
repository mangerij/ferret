[Mesh]
  file = slab_exodus_coarse_50.e
  uniform_refine = 0
  #distribution = serial
[]

#[Mesh]
#  type = GeneratedMesh
#  dim = 3
#  nx = 11
#  ny = 11
#  nz = 11
#  xmin = -5
#  xmax = 5
#  ymin = -5
#  ymax = 5
#  zmin = -5
#  zmax = 5
#[]

#[NodalNormals]
#  boundary = '3 4'
#[]


[GlobalParams]
   # #Can use a unit system where the differences between the coefficients
   # # in this table are minimized: [nm], [aC], [nN]
   # #length scale
   len_scale = 1.0
   # #BulkEnergy coefficients
   alpha1 = -0.0346849 # (3.766 * (T - 765.1) * 10^5) C^{-2} nm^2 (T = 650 K)
   alpha11 = -0.07253  #
   alpha111 = 0.26
   alpha12 = 0.75
   alpha112 = 0.61
   alpha123 = -3.7
   # #WallEnergy coefficients
   G110 = 0.25
   G11/G110 = 0.6
   G12/G110 = 0.0
   G44/G110 = 0.3
   G44P/G110 = 0.3
   # #Ferroelectric materials property
    # Q_mnkl = '8.9e-2 -2.6e-2 -2.6e-2 8.9e-2 -2.6e-2 8.9e-2 6.75e-2 6.75e-2 6.75e-2'
   # #Electrostatics
   permittivity = 0.008854187
   polar_x = polar_x
   polar_y = polar_y
   polar_z = polar_z
   potential_int = potential_int
  # potential_ext = potential_ext
   # #elastic variables
  # disp_x = disp_x
  # disp_y = disp_y
  ## disp_z = disp_z
  # displacements = 'disp_x disp_y'
  # displacements = 'disp_x disp_y disp_z'
  # use_displaced_mesh = false
    #C_ijkl = '3.80e-7 1.5e-7 1.50e-7 3.80e-7 1.50e-7 3.80e-7 1.1e-7 1.1e-7 1.1e-7'
   # #restart system
  # initial_from_file_timestep = 120
[]

#[GlobalParams]
# len_scale = 1e-9
#alpha1 = -1.0 # (3.766(T-765.1)*10^5) C^{-2} nm^2 (T = 673 K)
#alpha11 = -1.198310847153484
#alpha111 = 2.4615968080672657
#alpha12 = 12.39119171880757
#alpha112 = 5.775284818927047
#alpha123 = -35.0304161148034
#G110 = 1
#G11/G110 = 0.6
#G12/G110 = 0
#G44/G110 = 0.3
#G44P/G110 = 0.3
# Q_mnkl = '1 -0.29213483146067415 -0.29213483146067415 1 -0.29213483146067415 1 0.7584269662921348 0.7584269662921348 0.7584269662921348'
#permittivity = 8.854186999999999e-12
#polar_x = polar_x
#polar_y = polar_y
#polar_z = polar_z
#potential_int = potential_int
##potential_ext = potential_ext
#disp_x=disp_x
#disp_y=disp_y
#disp_z=disp_z
#displacements='disp_x disp_y disp_z'
#use_displaced_mesh = false
#C_ijkl = '1 0.3947368421052632 0.3947368421052632 1 0.3947368421052632 1 0.2894736842105264 0.2894736842105264 0.2894736842105264'
##initial_from_file_timestep = 120
#strain_scale = 0.051 #\episilon_{sf} \epsilon_{kl} = Q_11 * P_0^2 \epsilon_{kl}
#[]


[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block='2'
    #scaling = 1e
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
    block='2'
    #scaling = 1e11
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
    block='2'
    #scaling = 1e5
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
      seed = 1
    [../]
    #initial_from_file_var = polar_z
  [../]
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
  #[./disp_x]
  #  order = FIRST
  #  family = LAGRANGE
  #  [./InitialCondition]
  #    type = RandomIC
  #    min = -0.5e-5
  #    max = 0.5e-5
  #    seed = 1
  #  [../]
  #  #scaling = 1e5
  #  block = '2'
  #  #scaling = 1e6
  #  #[./InitialCondition]
  #  #  type = RandomIC
  #  #  min = -0.5e-6
  #  #  max = 0.5e-6
  #  #[../]
  #  #initial_from_file_var = disp_x
  #[../]
  #[./disp_y]
  #  order = FIRST
  #  family = LAGRANGE
  #  [./InitialCondition]
  #    type = RandomIC
  #    min = -0.5e-5
  #    max = 0.5e-5
  #    seed = 1
  #  [../]
  #  #scaling = 1e5
  #  block = '2'
  #  #scaling = 1e6
  #  #[./InitialCondition]
  #  #  type = RandomIC
  #  #  min = -0.5e-6
  #  #  max = 0.5e-6
  #  #[../]
  #  #initial_from_file_var = disp_y
  #[../]
  #[./disp_z]
  #  order = FIRST
  #  family = LAGRANGE
  #  block = '2'
  #  [./InitialCondition]
  #    type = RandomIC
  #    min = -0.5e-5
  #    max = 0.5e-5
  #    seed = 1
  #  [../]
  #  #scaling = 1e5
  #  #[./InitialCondition]
  #  #  type = RandomIC
  #  #  min = -0.5e-6
  #  #  max = 0.5e-6
  #  #[../]
  #  #initial_from_file_var = disp_z
  #[../]
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
#  [./stress_xy]
#    order = CONSTANT
#    family = MONOMIAL
#    #initial_from_file_var = stress_xy
#  [../]
#  [./stress_xz]
#    order = CONSTANT
#    family = MONOMIAL
#    #initial_from_file_var = stress_xz
#  [../]
#  [./stress_zz]
#    order = CONSTANT
#    family = MONOMIAL
#    #initial_from_file_var = stress_zz
#  [../]
#  [./stress_yz]
#    order = CONSTANT
#    family = MONOMIAL
#    #initial_from_file_var = stress_yz
#  [../]
  #[./elastic_strain_xx]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  #initial_from_file_var = strain_xx
  #[../]
  #[./elastic_strain_yy]
  #  order = CONSTANT
  #  family = MONOMIAL
  #  #initial_from_file_var = strain_yy
  #[../]
#  [./elastic_strain_xy]
#    order = CONSTANT
#    family = MONOMIAL
#    #initial_from_file_var = strain_xy
#  [../]
#  [./elastic_strain_xz]
#    order = CONSTANT
#    family = MONOMIAL
#    #initial_from_file_var = strain_xz
#  [../]
#  [./elastic_strain_zz]
#    order = CONSTANT
#    family = MONOMIAL
#    #initial_from_file_var = strain_zz
#  [../]
#  [./elastic_strain_yz]
#    order = CONSTANT
#    family = MONOMIAL
#    #initial_from_file_var = strain_yz
#  [../]
#  [./rho_b]
#    order = CONSTANT
#    family = MONOMIAL
#    #initial_from_file_var = rho_b
#  [../]
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
  #  block = '2'
  #  variable = stress_xx
  #  execute_on = 'timestep_end'
  #[../]
#  [./matl_s12]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    index_i = 0
#    index_j = 1
#    block = '2'
#    variable = stress_xy
#    execute_on = 'timestep_end'
#  [../]
#  [./matl_s13]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    index_i = 0
#    index_j = 2
#    block = '2'
#    variable = stress_xz
#    execute_on = 'timestep_end'
#  [../]
#  [./matl_s22]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    index_i = 1
#    index_j = 1
#    block = '2'
#    variable = stress_yy
#    execute_on = 'timestep_end'
#  [../]
#  [./matl_s23]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    index_i = 1
#    index_j = 2
#    block = '2'
#    variable = stress_yz
#    execute_on = 'timestep_end'
#  [../]
#  [./matl_s33]
#    type = RankTwoAux
#    rank_two_tensor = stress
#    index_i = 2
#    index_j = 2
#    block = '2'
#    variable = stress_zz
#    execute_on = 'timestep_end'
#  [../]
  #[./matl_e11]
  #  type = RankTwoAux
  #  rank_two_tensor = elastic_strain
  #  index_i = 0
  #  index_j = 0
  #  block = '2'
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
#  [./matl_e13]
#    type = RankTwoAux
#    rank_two_tensor = elastic_strain
#    index_i = 0
#    index_j = 2
#    block = '2'
#    variable = elastic_strain_xz
#    execute_on = 'timestep_end'
#  [../]
  #[./matl_e22]
  #  type = RankTwoAux
  #  rank_two_tensor = elastic_strain
  #  index_i = 1
  #  index_j = 1
  #  block = '2'
  #  variable = elastic_strain_yy
  #  execute_on = 'timestep_end'
  #[../]
#  [./matl_e23]
#    type = RankTwoAux
#    rank_two_tensor = elastic_strain
#    index_i = 1
#    index_j = 2
#    block = '2'
#    variable = elastic_strain_yz
#    execute_on = 'timestep_end'
#  [../]
#  [./matl_e33]
#    type = RankTwoAux
#    rank_two_tensor = elastic_strain
#    index_i = 2
#    index_j = 2
#    block = '2'
#    variable = elastic_strain_zz
#    execute_on = 'timestep_end'
#  [../]
[]

[Kernels]
  #Elastic problem
  #[./TensorMechanics] #This is an action block
  #[../]

  #[./stressdiv_0]
  #  type = StressDivergenceTensorsScaled
  #  variable = disp_x
  #  component = 0
  #  block = '2'
  #[../]
  #[./stressdiv_1]
  #  type = StressDivergenceTensorsScaled
  #  variable = disp_y
  #  component = 1
  #  block = '2'
  #[../]
  #[./stressdiv_2]
  #  type = StressDivergenceTensorsScaled
  #  variable = disp_z
  #  component = 2
  #  block = '2'
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
  #[./ferroelectriccouplingp_xx]
  #   type = FerroelectricCouplingP
  #   variable=polar_x
  #   component = 0
  #  # block = '2'
  #[../]
  #[./ferroelectriccouplingp_yy]
  #   type = FerroelectricCouplingP
  #   variable = polar_y
  #   component = 1
  #  # block = '2'
  #[../]
  #[./ferroelectriccouplingp_zz]
  #   type = FerroelectricCouplingP
  #   variable = polar_z
  #   component = 2
  #  # block = '2'
  #[../]
  ##Electrostatics
  [./polar_x_electric_E]
     type=PolarElectricEStrong
     variable = potential_int
     block='2'
    # implicit = false
  [../]
  [./E_int]
     type=Electrostatics
     variable = potential_int
     block='1 3'
  [../]
  [./FE_E_int]
     type=Electrostatics
     variable = potential_int
     block='2'
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
  # "infinity condition"
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
   #screening?
  # [./potential_int_upz]
  #   type = DirichletBC
  #   variable = potential_int
  #   boundary = '3'
  #   value = 0.001
  # [../]
  # [./potential_int_downz]
  #   type = DirichletBC
  #   variable = potential_int
  #   boundary = '4'
  #   value = 0.001
  # [../]
  # [./potential_int_top]
  #   type = NeumannBC
  #   variable = potential_int
  #   boundary = '3'
  #   value = 0.0
  # [../]
  # [./potential_int_bottom]
  #   type = NeumannBC
  #   variable = potential_int
  #   boundary = '4'
  #   value = 0.0
  # [../]
  # [./potential_int_right]
  #   type = DirichletBC
  #   variable = potential_int
  #   boundary = '5'
  #   value = 0.0
  # [../]
  # [./potential_int_front]
  #   type = DirichletBC
  #   variable = potential_int
  #   boundary = '6'
  #   value = 0.0
  # [../]
  # [./potential_int_left]
  #   type = DirichletBC
  #   variable = potential_int
  #   boundary = '7'
  #   value = 0.0
  # [../]
  # [./potential_int_back]
  #   type = DirichletBC
  #   variable = potential_int
  #   boundary = '8'
  #   value = 0.0
  # [../]

  # [./depol_bottom]
  #   type = DepolScreenBC
  #   variable = potential_int
  #   boundary = '4'
  #   value = 0.99
  # [../]


  #[./disp_x_cube1]
  #  type = DirichletBC
  #  boundary = 'top'
  #  value = 0.0
  #  variable = disp_x
  #[../]
  #[./disp_x_cube2]
  #  type = DirichletBC
  #  boundary = 'bottom'
  #  value = 0.0
  #  variable = disp_x
  #[../]
  #[./disp_x_cube3]
  #  type = DirichletBC
  #  boundary = 'front'
  #  value = 0.0
  #  variable = disp_x
  #[../]
  #[./disp_x_cube4]
  #  type = DirichletBC
  #  boundary = 'back'
  #  value = 0.0
  #  variable = disp_x
  #[../]
  #[./disp_x_cube5]
  #  type = DirichletBC
  #  boundary = 'left'
  #  value = 0.0
  #  variable = disp_x
  #[../]
  #[./disp_x_cube6]
  #  type = DirichletBC
  #  boundary = 'right'
  #  value = 0.0
  #  variable = disp_x
  #[../]
  #
  #[./disp_y_cube1]
  #  type = DirichletBC
  #  boundary = 'top'
  #  value = 0.0
  #  variable = disp_y
  #[../]
  #[./disp_y_cube2]
  #  type = DirichletBC
  #  boundary = 'bottom'
  #  value = 0.0
  #  variable = disp_y
  #[../]
  #[./disp_y_cube3]
  #  type = DirichletBC
  #  boundary = 'front'
  #  value = 0.0
  #  variable = disp_y
  #[../]
  #[./disp_y_cube4]
  #  type = DirichletBC
  #  boundary = 'back'
  #  value = 0.0
  #  variable = disp_y
  #[../]
  #[./disp_y_cube5]
  #  type = DirichletBC
  #  boundary = 'left'
  #  value = 0.0
  #  variable = disp_y
  #[../]
  #[./disp_y_cube6]
  #  type = DirichletBC
  #  boundary = 'right'
  #  value = 0.0
  #  variable = disp_y
  #[../]
  #
  #[./disp_z_cube1]
  #  type = DirichletBC
  #  boundary = 'top'
  #  value = 0.0
  #  variable = disp_z
  #[../]
  #[./disp_z_cube2]
  #  type = DirichletBC
  #  boundary = 'bottom'
  #  value = 0.0
  #  variable = disp_z
  #[../]
  #[./disp_z_cube3]
  #  type = DirichletBC
  #  boundary = 'front'
  #  value = 0.0
  #  variable = disp_z
  #[../]
  #[./disp_z_cube4]
  #  type = DirichletBC
  #  boundary = 'back'
  #  value = 0.0
  #  variable = disp_z
  #[../]
  #[./disp_z_cube5]
  #  type = DirichletBC
  #  boundary = 'left'
  #  value = 0.0
  #  variable = disp_z
  #[../]
  #[./disp_z_cube6]
  #  type = DirichletBC
  #  boundary = 'right'
  #  value = 0.0
  #  variable = disp_z
  #[../]


  #[./polar_y_cube1]
  #  type = DirichletBC
  #  boundary = 'top'
  #  value = 0.3
  #  variable = polar_y
  #[../]
  #[./polar_y_cube2]
  #  type = DirichletBC
  #  boundary = 'bottom'
  #  value = -0.3
  #  variable = polar_y
  #[../]
  #[./disp_y_cube3]
  #  type = DirichletBC
  #  boundary = 'front'
  #  value = 0.0
  #  variable = disp_y
  #[../]
  #[./disp_y_cube4]
  #  type = DirichletBC
  #  boundary = 'back'
  #  value = 0.0
  #  variable = disp_y
  #[../]
  #[./potential_cube5]
  #  type = DirichletBC
  #  boundary = 'left'
  #  value = 0.00
  #  variable = potential_int
  #[../]
  #[./potential_cube6]
  #  type = DirichletBC
  #  boundary = 'right'
  #  value = 0.00
  #  variable = potential_int
  #[../]
   #Top and bottom {3, 4}:
  # [./disp_x_slab3]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = '3'
  #   value = 0.0
  # [../]
  # [./disp_y_slab3]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = '3'
  #   value = 0.0
  # [../]
  # [./disp_z_slab3]
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = '3'
  #   value = 0.0
  # [../]
  # [./disp_x_slab4]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = '4'
  #   value = 0.0
  # [../]
  # [./disp_y_slab4]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = '4'
  #   value = 0.0
  # [../]
  # [./disp_z_slab4]
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = '4'
  #   value = 0.0
  # [../]

  #Sides {5+/-y, 6+/-x, 7-/+y, 8-/+x}:
  #[./disp_y_slab5]
  #  type = DirichletBC
  #  variable = disp_y
  #  boundary = '5'
  #  value = 0.4
  #[../]
  #[./disp_y_slab5]
  #  type = DirichletBC
  #  variable = disp_y
  #  boundary = '5'
  #  value = 0.0
  #[../]
  #[./disp_z_slab5]
  #  type = DirichletBC
  #  variable = disp_z
  #  boundary = '5'
  #  value = 0.0
  #[../]
  #[./disp_x_slab5]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = '5'
  #   value = 0.0
  #[../]
  #[./disp_y_slab5]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = '5'
  #   value = 0.0
  #[../]
  #[./disp_z_slab5]
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = '5'
  #   value = 0.0
  #[../]
  #
  #
  #
  #[./disp_x_slab6]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = '6'
  #   value = 0.0
  #[../]
  #[./disp_y_slab6]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = '6'
  #   value = 0.0
  #[../]
  #[./disp_z_slab6]
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = '6'
  #   value = 0.0
  #[../]
  #[./disp_x_slab7]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = '7'
  #   value = 0.0
  #[../]
  #[./disp_y_slab7]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = '7'
  #   value = 0.0
  #[../]
  #[./disp_z_slab7]
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = '7'
  #   value = 0.0
  #[../]
  #[./disp_x_slab8]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = '8'
  #   value = 0.0
  #[../]
  #[./disp_y_slab8]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = '8'
  #   value = 0.0
  #[../]
  #[./disp_z_slab8]
  #   type = DirichletBC
  #   variable = disp_z
  #   boundary = '8'
  #   value = 0.0
  #[../]
  #[./potential_int_upz]
  #  type = DirichletBC
  #  variable = potential_int
  #  boundary = '1'
  #  value = 0.0
  #[../]
  #[./potential_int_downz]
  #  type = DirichletBC
  #  variable = potential_int
  #  boundary = '2'
  #  value = 0.0
  #[../]
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
  #  #C_ijkl = '3.80e-7 1.5e-7 1.50e-7 3.80e-7 1.50e-7 3.80e-7 1.1e-7 1.1e-7 1.1e-7'
  #  #in m^4/C^2. from http://arxiv.org/pdf/1205.5640.pdf
  #  # Q11 Q12 Q13 Q22 Q23 Q33 Q44 Q55 Q66
  #  #euler_angle_1 = 0.0 #currently will only rotate C_ijkl
  #  #euler_angle_2 = 0.0
  #  #euler_angle_3 = 0.0
  #[../]
  #[./elasticity_tensor]
  #  type = ComputeElasticityTensor
  #  block = '2'
  ##  block = '2'
  #  fill_method = symmetric9
  #[../]
  #[./strain]
  #  type = ComputeSmallStrain
  #  block = '2'
  ##  block = '2'
  #[../]
  #[./stress]
  #  type = ComputeLinearElasticStress
  #  block = '2'
  ##  block = '2'
  #[../]
  #[./vacuum]
  #  type=GenericConstantMaterial
  #  block = '1 3'
  #[../]
  ###This seems to be what we want for a simple epitaxial test
  ### (note that most epitaxial conditions are a strain gradient from the interface)
  #[./eigen_strain_xx_yy] #Use for stress-free strain (ie epitaxial)
  #  type = ComputeEigenstrain
  #  block = '2'
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
   block = '2'
  [../]
  [./wall_energy]
   type = WallEnergy
   execute_on = 'timestep_end'
  # initial_from_file_var = wall_energy
   block = '2'
  [../]
  #[./elastic_energy]
  # type = ElasticEnergy
  # execute_on = 'timestep_end'
  ## initial_from_file_var = elastic_energy
  # block = '2'
  #[../]
  #[./coupled_energy]
  #  type = CoupledEnergy
  #  execute_on = 'timestep_end'
  #  block = '2'
  #[../]
  [./electrostatic_energy]
   type = ElectrostaticEnergy
   execute_on = 'timestep_end'
   block = '2'
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

[UserObjects]
  [./kill]
    type = Terminator
    expression = 'perc_change <= 8.5e-6'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason -options_left'
    petsc_options_iname = '-ksp_gmres_restart -snes_rtol -ksp_rtol -pc_type -sub_pc_type -sub_pc_factor_zeropivot -pc_factor_zeropivot -pc_side'
    petsc_options_value = '    101             1e-8     1e-10        bjacobi       lu                   1e-50                 1e-50              left'
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
  #dtmin = 1.0e-10
  dt = 3.0
  dtmin = 1e-30
  dtmax = 3.0
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
    file_base = out_PbTiO3_50nm_T673K_implicit_FE_prob
    output_initial = true
    elemental_as_nodal = true
    interval = 1
  [../]
  [./debug]
    type = VariableResidualNormsDebugOutput
  [../]
[]
