
[Mesh]
  type = GeneratedMesh
  dim = 3
 
  # number of elements along each of the cartesian directions
  # NOTE: the gradient term needs ~1-2 elements across domain walls
  nx = 11
  ny = 11
  nz = 11


  # dimensions of the total computational volume in
  # nanometers
  xmin = -4.0
  xmax = 4.0

  ymin = -4.0
  ymax = 4.0

  zmin = -4.0
  zmax = 4.0
[]


[MeshModifiers]
  [./subdomain1]
    type = SubdomainBoundingBox
    bottom_left = '-4.0 -4.0 0.0'
    top_right = '4.0 4.0 4.0'
    block_id = 1
  [../]
  [./break_boundary]
    depends_on = subdomain1
    type = BreakBoundaryOnSubdomain
  [../]
  [./interface]
    type = SideSetsBetweenSubdomains
    depends_on = break_boundary
    master_block = '0'
    paired_block = '1'
    new_boundary = 'master0_interface'
  [../]

  # need to generate sidesets with SideSetsFromNormals
  # [./boundary_1]
  # [../]
  # ect...
[]


[GlobalParams]
  # default length scale
  # Note that all units are in nanometers
  #                            nanoNewtons
  #                            attoCoulombs
  # By doing this, this allows for optimal
  # convergence of the LGD solve
  len_scale = 1.0

  # Landau expansion parameters
  alpha1 = -0.1722883 # Room temperature
  alpha11 = -0.07253
  alpha111 = 0.26
  alpha12 = 0.75
  alpha112 = 0.61
  alpha123 = -3.67

  # Gradient energy coefficients
  G110 = 0.173
  G11/G110 = 0.6
  G12/G110 = 0
  G44/G110 = 0.3
  G44P/G110 = 0.3


  # Some default stuff. 
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int
  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  displacements = 'disp_x disp_y disp_z'

  # Epitaxial eigenstrain
  # prefactor = 0.01 #negative = tension, positive = compression
[]



[Variables]
  # Here we define the variables that we want to solve over the mesh
  # polar field polar_k, displacement field disp_k, and electrostatic
  # potential_int (for k = x,y,z).
  [./polar_x]
    block = '1'
    order = FIRST
    family = LAGRANGE
    # Note that the initial condition of the problem is near the paraelectric
    # This ensures that the domain structure that forms is not biased.
    # See L.-Q. Chen's APL for a discussion.
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./polar_y]
    block = '1'
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./polar_z]
    block = '1'
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
    block = '0 1'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-5
      max = 0.5e-5
    [../]
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]


[AuxVariables]
  # These are auxilliary variables of the problem
  # This is how you store stress/strain (locally)
  # You can also calculate anything you want (locally) 
  # using the above variables (P, u, \phi)

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
 # The above block stores the variables, while
 # this block calculates them.
  [./matl_e11]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    variable = strain_xx_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e12]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    variable = strain_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e13]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 2
    variable = strain_xz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e22]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    variable = strain_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e23]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 2
    variable = strain_yz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_e33]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    index_i = 2
    index_j = 2
    variable = strain_zz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s11]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = stress_xx_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s12]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = stress_xy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s13]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
    variable = stress_xz_elastic
    execute_on = 'timestep_end'
  [../]
 [./matl_s22]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = stress_yy_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s23]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    variable = stress_yz_elastic
    execute_on = 'timestep_end'
  [../]
  [./matl_s33]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    variable = stress_zz_elastic
    execute_on = 'timestep_end'
  [../]
[]

[Materials]
  [./eigen_strain_zz] #Use for stress-free strain (ie epitaxial)
   type = ComputeEigenstrain
   block = '0 1'
  # eigen_base = 'exx exy exz eyx eyy eyz ezx ezy ezz'
   eigen_base = '1 0 0 0 1 0 0 0 0'
   eigenstrain_name = 'eigenstrain'
 [../]


  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
    block = '1'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    block = '1'
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
    block = '1'
  [../]

  [./slab_ferroelectric]
    block = '1'
    type = ComputeElectrostrictiveTensor
    Q_mnkl = '0.089 -0.026 -0.026 0.089 -0.026 0.089 0.03375 0.03375 0.03375'
    C_ijkl = '380. 150. 150. 380. 150. 380. 110. 110. 110.'
  [../]

  [./elasticity_tensor_2]
    type = ComputeElasticityTensor
    C_ijkl = '319 99.6 99.6 319 99.6 319 109.53 109.53 109.53'
    fill_method = symmetric9
    block = '0'
  [../]
  [./strain_2]
    type = ComputeSmallStrain
    block = '0'
  [../]
  [./stress_2]
    type = ComputeLinearElasticStress
    block = '0'
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


  [./ferroelectriccouplingX_xx]
    type = FerroelectricCouplingX
    block = '1'
    variable = disp_x
    component = 0
  [../]
  [./ferroelectriccouplingX_yy]
    type = FerroelectricCouplingX
    block = '1'
    variable = disp_y
    component = 1
  [../]
  [./ferroelectriccouplingX_zz]
    type = FerroelectricCouplingX
    block = '1'
    variable = disp_z
    component = 2
  [../]
  ##Electrostatics
  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
     permittivity = 0.08854187
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     block = '1'
     permittivity = 0.08854187
  [../]

  [./DIE_E_int]
     type = Electrostatics
     variable = potential_int
     block  = '0'
     permittivity = 2.6562561
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


#[BCs]
#  [./bot_disp_x]
#    variable = disp_x
#    type = DirichletBC
#    value = 0
#    boundary = '7'
#  [../]
#  [./bot_disp_y]
#    variable = disp_y
#    type = DirichletBC
#    value = 0
#    boundary = '7'
#  [../]
#  [./bot_disp_z]
#    variable = disp_z
#    type = DirichletBC
#    value = 0
#    boundary = '7'
#  [../]

#  [./bot_potential_int]
#    variable = disp_x
#    type = DirichletBC
#    value = 0
#    boundary = '7'
#  [../]

#  [./Periodic]
#    [./TB_disp_x_pbc]
#      variable = disp_x
#      primary = '3'
#      secondary = '5'
#      translation = '0 20 0'
#    [../]
#    [./TB_disp_y_pbc]
#      variable = disp_y
#      primary = '3'
#      secondary = '5'
#      translation = '0 20 0'
#    [../]
#    [./TB_disp_z_pbc]
#      variable = disp_z
#      primary = '3'
#      secondary = '5'
#      translation = '0 20 0'
#    [../]

#    [./TB_polar_x_pbc]
#      variable = polar_x
#      primary = '3'
#      secondary = '5'
#      translation = '0 20 0'
#    [../]
#    [./TB_polar_y_pbc]
#      variable = polar_y
#      primary = '3'
#      secondary = '5'
#      translation = '0 20 0'
#    [../]
#    [./TB_polar_z_pbc]
#      variable = polar_z
#      primary = '3'
#      secondary = '5'
#      translation = '0 20 0'
#    [../]
#    [./TB_potential_int_pbc]
#      variable = potential_int
#      primary = '3'
#      secondary = '5'
#      translation = '0 20 0'
#    [../]
  #

#    [./RL_disp_x_pbc]
#      variable = disp_x
#      primary = '4'
#      secondary = '6'
#      translation = '20 0 0'
#    [../]
#    [./RL_disp_y_pbc]
#      variable = disp_y
#      primary = '4'
#      secondary = '6'
#      translation = '20 0 0'
#    [../]
#    [./RL_disp_z_pbc]
#      variable = disp_z
#      primary = '4'
#      secondary = '6'
#      translation = '20 0 0'
#    [../]

#    [./RL_polar_x_pbc]
#      variable = polar_x
#      primary = '4'
#      secondary = '6'
#      translation = '20 0 0'
#    [../]
#    [./RL_polar_y_pbc]
#      variable = polar_y
#      primary = '4'
#      secondary = '6'
#      translation = '20 0 0'
#    [../]
#    [./RL_polar_z_pbc]
#      variable = polar_z
#      primary = '4'
#      secondary = '6'
#      translation = '20 0 0'
#    [../]
#    [./RL_potential_int_pbc]
#      variable = potential_int
#      primary = '4'
#      secondary = '6'
#      translation = '20 0 0'
#    [../]
#  [../]
#[]


[Postprocessors]
   [./Fbulk]
      type = BulkEnergy
      block = '1'
      execute_on = 'timestep_end'
    [../]
    [./Fwall]
      type = WallEnergy
      block = '1'
      execute_on = 'timestep_end'
    [../]
    [./Felastic]
      type = ElasticEnergy
      block = '0 1'
      execute_on = 'timestep_end'
    [../]
    [./Fcoupled]
      block = '1'
      type = CoupledEnergy
      execute_on = 'timestep_end'
    [../]
    [./Felec]
      block = '1'
      type = ElectrostaticEnergy
      permittivity = 0.08854187
      execute_on = 'timestep_end'
    [../]
    [./Ftotal]
      type = TotalEnergyFlow
      Fbulk = Fbulk
      Fwall = Fwall
      Fcoupled = Fcoupled
      Felec = Felec
      execute_on = 'timestep_end'
    [../]
    [./perc_change]
     type = PercentChangePostprocessor
     postprocessor = Ftotal
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
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121                1e-6      1e-8    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
    [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.8
    #iteration_window = 3
    optimal_iterations = 6 #should be 5 probably
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.8
[../]
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  #dt = 0.5
  dtmin = 1e-13
  dtmax = 0.8
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = outPTO_test
    elemental_as_nodal = true
    interval = 1
  [../]
[]
