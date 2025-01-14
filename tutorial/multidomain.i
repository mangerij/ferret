
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 32
  ny = 32
  nz = 6
  xmin = -15
  xmax = 15
  ymin = -15
  ymax = 15
  zmin = -2
  zmax = 2
  elem_type = HEX8
[]



[GlobalParams]
  len_scale = 1.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  potential_E_int = potential_int
[]

[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.0001
      max = 0.0001
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.0001
      max = 0.0001
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.0001
      max = 0.0001
    [../]
  [../]

  [./potential_int]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]

  #---------------------------------------#
  #                                       #
  #     Bulk (homogeneous) free energy    #
  #                                       #
  #---------------------------------------#

  [./bed_x]
    type = BulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = BulkEnergyDerivativeEighth
    variable = polar_z
    component = 2
  [../]

  #---------------------------------------#
  #                                       #
  #     Gradient free energy              #
  #                                       #
  #---------------------------------------#

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

  [./walled2_x]
    type = Wall2EnergyDerivative
    variable = polar_x
    component = 0
  [../]
  [./walled2_y]
    type = Wall2EnergyDerivative
    variable = polar_y
    component = 1
  [../]
  [./walled2_z]
    type = Wall2EnergyDerivative
    variable = polar_z
    component = 2
  [../]

  #---------------------------------------#
  #                                       #
  #     Poissons equation and P*E         #
  #                                       #
  #---------------------------------------#

  [./polar_x_electric_E]
    type = PolarElectricEStrong
    variable = potential_int
  [../]
  [./FE_E_int]
    type = Electrostatics
    variable = potential_int
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

  #---------------------------------------#
  #                                       #
  #     Time dependence of Pj             #
  #                                       #
  #---------------------------------------#

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

[Materials]

  ## Use coefficients for lead-titanate in this convienient unit system (aC, kg, sec, nm)

  [./Landau_P_FE]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '-0.1722883 -0.073 0.75 0.26 0.61 -3.67 0.0 0.0 0.0 0.0'
    block = '0'
  [../]
  [./Landau_G_FE]
    type = GenericConstantMaterial
    prop_names = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '0.173 0.6 0.0 0.3 0.3'
    block = '0'
  [../]

  [./permitivitty_1]

    ###############################################
    ##
    ##  so-called background dielectric constant
    ##  (it encapsulates the motion of core electrons
    ##  at high frequency) = e_b*e_0 (here we use
    ##  e_b = 10), see PRB. 74, 104014, (2006)
    ##  This does not influence results much.
    ##
    ###############################################

    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '0.08854187'
  [../]
[]

[BCs]

[]

[Postprocessors]
  [./avePz]
    type = ElementAverageValue
    variable = polar_z
    execute_on = 'timestep_end'
  [../]
  [./FbP]
    type = BulkEnergyEighth
    execute_on = 'initial timestep_end'
  [../]
  [./FgP]
    type = WallEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'FbP FgP'
    pp_coefs = ' 1 1'
    execute_on = 'timestep_end'

    ##########################################
    #
    # NOTE: Ferret output is in attojoules
    #
    ##########################################
  [../]
  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = Ftot
  [../]
[]

[UserObjects]
  [./kill]
    type = Terminator
    expression = 'perc_change <= 1.0e-4'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    120               1e-10      1e-8     1e-4    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  scheme = 'implicit-euler'
  dtmax = 0.7
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    execute_on = 'timestep_end'
    file_base = out_multidomain
    elemental_as_nodal = true
    interval = 1
  [../]
  [./outcsv]
    type = CSV
    file_base = out_multidomain
    execute_on = 'timestep_end'
  [../]
[]
