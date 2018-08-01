
#alpha1 alpha3
alpha1 = '-0.0827836 -0.250056'
#alpha11 alpha33
alpha11 = '0.42229 0.04990909090909092'
#alpha12 alpha13
alpha12 = '0.734181277056277 0.4521818181818182'
alpha111 = 0.0#0.26
alpha112 = 0.0#0.61
alpha123 = 0.0#-3.7
alpha1_para = '0.1854662 0.1854662'
alpha11_para = '0.0 0.0'
alpha12_para = '0.0 0.0'
alpha111_para = 0.0
alpha112_para = 0.0
alpha123_para = 0.0
G110 = 0.173
G11_G110 = 1.6
G12_G110 = 0
G44_G110 = 1.6
G44P_G110 = 1.6
G110_para = 0.173
G11_G110_para = 1.6
G12_G110_para = 0
G44_G110_para = 1.6
G44P_G110_para = 1.6
permittivity_electrostatic = 0.0885
permittivity_electrostatic_para = 0.0885

permitivitty_depol = 0.00885
permitivitty_depol_para = 0.00885

[Mesh]
  type = FileMesh
  file = template_mesh.msh
[]

[GlobalParams]
  len_scale = 1.0
[]

[Variables]
  [./FE_polar_x]
    order = FIRST
    family = LAGRANGE
    block = 'ferro_volume'
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-5
      max = 0.01e-5
      seed = 5
    [../]
  [../]
  [./FE_polar_y]
    order = FIRST
    family = LAGRANGE
    block = 'ferro_volume'
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-5
      max = 0.01e-5
      seed = 5
    [../]
  [../]
  [./FE_polar_z]
    order = FIRST
    family = LAGRANGE
    block = 'ferro_volume'
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-5
      max = 0.01e-5
      seed = 5
    [../]
  [../]
  [./PE_polar_x]
    order = FIRST
    family = LAGRANGE
    block = 'para_bottom_volume para_top_volume'
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-5
      max = 0.01e-5
      seed = 5
    [../]
  [../]
  [./PE_polar_y]
    order = FIRST
    family = LAGRANGE
    block = 'para_bottom_volume para_top_volume'
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-5
      max = 0.01e-5
      seed = 5
    [../]
  [../]
  [./PE_polar_z]
    order = FIRST
    family = LAGRANGE
    block = 'para_bottom_volume para_top_volume'
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-5
      max = 0.01e-5
      seed = 5
    [../]
  [../]
  [./potential_E_int]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.01e-5
      max = 0.01e-5
      seed = 5
    [../]
  [../]
[]

[Kernels]
  #FERROELECTRIC BLOCK
  [./bed_x_ferro]
    type = BulkEnergyDerivativeSixth
    block = 'ferro_volume'
    variable = FE_polar_x
    polar_x = FE_polar_x
    polar_y = FE_polar_y
    polar_z = FE_polar_z
    alpha1 = ${alpha1}
    alpha11 = ${alpha11}
    alpha12 = ${alpha12}
    alpha111 = ${alpha111}
    alpha112 = ${alpha112}
    alpha123 = ${alpha123}
    component = 0
  [../]
  [./bed_y_ferro]
    type = BulkEnergyDerivativeSixth
    block = 'ferro_volume'
    variable = FE_polar_y
    polar_x = FE_polar_x
    polar_y = FE_polar_y
    polar_z = FE_polar_z
    alpha1 = ${alpha1}
    alpha11 = ${alpha11}
    alpha12 = ${alpha12}
    alpha111 = ${alpha111}
    alpha112 = ${alpha112}
    alpha123 = ${alpha123}
    component = 1
  [../]
  [./bed_z_ferro]
    type = BulkEnergyDerivativeSixth
    block = 'ferro_volume'
    variable = FE_polar_z
    polar_x = FE_polar_x
    polar_y = FE_polar_y
    polar_z = FE_polar_z
    alpha1 = ${alpha1}
    alpha11 = ${alpha11}
    alpha12 = ${alpha12}
    alpha111 = ${alpha111}
    alpha112 = ${alpha112}
    alpha123 = ${alpha123}
    component = 2
  [../]
  [./polar_x_electric_E_ferro]
     type = PolarElectricEStrong
       block = 'ferro_volume'
       polar_x = FE_polar_x
       polar_y = FE_polar_y
       polar_z = FE_polar_z
       variable = potential_E_int
  [../]
  [./FE_E_int_ferro]
       type = Electrostatics
       block = 'ferro_volume'
       variable = potential_E_int
       permittivity = ${permittivity_electrostatic}
  [../]
  [./polar_electric_px_ferro]
     type = PolarElectricPStrong
       block = 'ferro_volume'
       variable = FE_polar_x
       potential_E_int = potential_E_int
       component = 0
  [../]
  [./polar_electric_py_ferro]
     type = PolarElectricPStrong
       block = 'ferro_volume'
       variable = FE_polar_y
       potential_E_int = potential_E_int
       component = 1
  [../]
  [./polar_electric_pz_ferro]
     type = PolarElectricPStrong
       block = 'ferro_volume'
       variable = FE_polar_z
       potential_E_int = potential_E_int
       component = 2
  [../]
  [./walled_x_ferro]
    type = WallEnergyDerivativeAlt
    block = 'ferro_volume'
    variable = FE_polar_x
    polar_x = FE_polar_x
    polar_y = FE_polar_y
    polar_z = FE_polar_z
    G110 = ${G110}
    G11_G110 = ${G11_G110}
    G12_G110 = ${G12_G110}
    G44_G110 = ${G44_G110}
    G44P_G110 = ${G44P_G110}
    component = 0
  [../]
  [./walled_y_ferro]
    type = WallEnergyDerivativeAlt
    block = 'ferro_volume'
    variable = FE_polar_y
    polar_x = FE_polar_x
    polar_y = FE_polar_y
    polar_z = FE_polar_z
    G110 = ${G110}
    G11_G110 = ${G11_G110}
    G12_G110 = ${G12_G110}
    G44_G110 = ${G44_G110}
    G44P_G110 = ${G44P_G110}
    component = 1
  [../]
  [./walled_z_ferro]
    type = WallEnergyDerivativeAlt
    block = 'ferro_volume'
    variable = FE_polar_z
    polar_x = FE_polar_x
    polar_y = FE_polar_y
    polar_z = FE_polar_z
    G110 = ${G110}
    G11_G110 = ${G11_G110}
    G12_G110 = ${G12_G110}
    G44_G110 = ${G44_G110}
    G44P_G110 = ${G44P_G110}
    component = 2
  [../]
  [./polar_x_time_ferro]
     type = TimeDerivativeScaled
     block = 'ferro_volume'
     variable = FE_polar_x
  [../]
  [./polar_y_time_ferro]
     type = TimeDerivativeScaled
     block = 'ferro_volume'
     variable = FE_polar_y
  [../]
  [./polar_z_time_ferro]
     type = TimeDerivativeScaled
     block = 'ferro_volume'
     variable = FE_polar_z
  [../]
  
  #PARAELECTRIC BLOCK
  [./bed_xp_para]
    alpha1111 = 0
    alpha1112 = 0
    alpha1122 = 0
    alpha1123 = 0
    type = BulkEnergyDerivativeSixth
    block = 'para_bottom_volume para_top_volume'
    variable = PE_polar_x
    polar_x = PE_polar_x
    polar_y = PE_polar_y
    polar_z = PE_polar_z
    alpha1 = ${alpha1_para}
    alpha11 = ${alpha11_para}
    alpha12 = ${alpha12_para}
    alpha111 = ${alpha111_para}
    alpha112 = ${alpha112_para}
    alpha123 = ${alpha123_para}
    component = 0
  [../]
  [./bed_yp_para]
    alpha1111 = 0
    alpha1112 = 0
    alpha1122 = 0
    alpha1123 = 0
    type = BulkEnergyDerivativeSixth
    block = 'para_bottom_volume para_top_volume'
    variable = PE_polar_y
    polar_x = PE_polar_x
    polar_y = PE_polar_y
    polar_z = PE_polar_z
    alpha1 = ${alpha1_para}
    alpha11 = ${alpha11_para}
    alpha12 = ${alpha12_para}
    alpha111 = ${alpha111_para}
    alpha112 = ${alpha112_para}
    alpha123 = ${alpha123_para}
    component = 1
  [../]
  [./bed_zp_para]
    alpha1111 = 0
    alpha1112 = 0
    alpha1122 = 0
    alpha1123 = 0
    type = BulkEnergyDerivativeSixth
    block = 'para_bottom_volume para_top_volume'
    variable = PE_polar_z
    polar_x = PE_polar_x
    polar_y = PE_polar_y
    polar_z = PE_polar_z
    alpha1 = ${alpha1_para}
    alpha11 = ${alpha11_para}
    alpha12 = ${alpha12_para}
    alpha111 = ${alpha111_para}
    alpha112 = ${alpha112_para}
    alpha123 = ${alpha123_para}
    component = 2
  [../]

  [./polar_x_electric_Ep_para]
     type = PolarElectricEStrong
     block = 'para_bottom_volume para_top_volume'
     polar_x = PE_polar_x
     polar_y = PE_polar_y
     polar_z = PE_polar_z
     variable = potential_E_int
  [../]
  [./FE_E_intp_para]
     type = Electrostatics
     block = 'para_bottom_volume para_top_volume'
     variable = potential_E_int
     permittivity = ${permittivity_electrostatic_para}
  [../]
  
  [./polar_electric_pxp_para]
     type = PolarElectricPStrong
     block = 'para_bottom_volume para_top_volume'
     variable = PE_polar_x
     potential_E_int = potential_E_int
     component = 0
  [../]
  [./polar_electric_pyp_para]
     type = PolarElectricPStrong
     block = 'para_bottom_volume para_top_volume'
     variable = PE_polar_y
     potential_E_int = potential_E_int
     component = 1
  [../]
  [./polar_electric_pzp_para]
     type = PolarElectricPStrong
     block = 'para_bottom_volume para_top_volume'
     variable = PE_polar_z
     potential_E_int = potential_E_int
     component = 2
  [../]
  
  [./walled_xp_para]
    type = WallEnergyDerivativeAlt
    block = 'para_bottom_volume para_top_volume'
    variable = PE_polar_x
    polar_x = PE_polar_x
    polar_y = PE_polar_y
    polar_z = PE_polar_z
    G110 = ${G110_para}
    G11_G110 = ${G11_G110_para}
    G12_G110 = ${G12_G110_para}
    G44_G110 = ${G44_G110_para}
    G44P_G110 = ${G44P_G110_para}
    component = 0
  [../]
  [./walled_yp_para]
    type = WallEnergyDerivativeAlt
    block = 'para_bottom_volume para_top_volume'
    variable = PE_polar_y
    polar_x = PE_polar_x
    polar_y = PE_polar_y
    polar_z = PE_polar_z
    G110 = ${G110_para}
    G11_G110 = ${G11_G110_para}
    G12_G110 = ${G12_G110_para}
    G44_G110 = ${G44_G110_para}
    G44P_G110 = ${G44P_G110_para}
    component = 1
  [../]
  [./walled_zp_para]
    type = WallEnergyDerivativeAlt
    block = 'para_bottom_volume para_top_volume'
    variable = PE_polar_z
    polar_x = PE_polar_x
    polar_y = PE_polar_y
    polar_z = PE_polar_z
    G110 = ${G110_para}
    G11_G110 = ${G11_G110_para}
    G12_G110 = ${G12_G110_para}
    G44_G110 = ${G44_G110_para}
    G44P_G110 = ${G44P_G110_para}
    component = 2
  [../]
  
  [./polar_x_timep_para]
     type = TimeDerivativeScaled
     block = 'para_bottom_volume para_top_volume'
     variable = PE_polar_x
  [../]
  [./polar_y_timep_para]
     type = TimeDerivativeScaled
     block = 'para_bottom_volume para_top_volume'
     variable = PE_polar_y
  [../]
  [./polar_z_timep_para]
     type = TimeDerivativeScaled
     block = 'para_bottom_volume para_top_volume'
     variable = PE_polar_z
  [../]
[]

[BCs]
   [./top_phi]
      type = DirichletBC
      variable = potential_E_int
      value = 0.0001
      boundary = 'para_top_surface_top'
   [../]
   [./bottom_phi]
      type = DirichletBC
      variable = potential_E_int
      value = 0.0
      boundary = 'para_bottom_surface_bottom'
   [../]   
   
  [./Periodic]
   [./per1]
     variable = 'PE_polar_x PE_polar_y PE_polar_z potential_E_int'
     primary = 'para_top_side_1'
     secondary = 'para_top_side_3'
     translation = '0 1 0'
   [../]
   [./per2]
     variable = 'PE_polar_x PE_polar_y PE_polar_z potential_E_int'
     primary = 'para_bottom_side_1'
     secondary = 'para_bottom_side_3'
     translation = '0 1 0'
   [../]
   [./per3]
     variable = 'PE_polar_x PE_polar_y PE_polar_z potential_E_int'
     primary = 'para_top_side_4'
     secondary = 'para_top_side_2'
     translation = '1 0 0'
   [../]
   [./per4]
     variable = 'PE_polar_x PE_polar_y PE_polar_z potential_E_int'
     primary = 'para_bottom_side_4'
     secondary = 'para_bottom_side_2'
     translation = '1 0 0'
   [../]
   [./per5]
     variable = 'FE_polar_x FE_polar_y FE_polar_z potential_E_int'
     primary = 'ferro_side_4'
     secondary = 'ferro_side_2'
     translation = '1 0 0'
   [../]
   [./per6]
     variable = 'FE_polar_x FE_polar_y FE_polar_z potential_E_int'
     primary = 'ferro_side_1'
     secondary = 'ferro_side_3'
     translation = '0 1 0'
   [../]
  [../]
[]

[Postprocessors]
  [./aveFEPz]
    type = ElementAverageValue
    variable = FE_polar_z
    execute_on = 'initial timestep_end'
    block = 'ferro_volume'
    #'initial linear nonlinear timestep_begin timestep_end'
  [../]
  [./avePEPz]
    type = ElementAverageValue
    variable = PE_polar_z
    execute_on = 'initial timestep_end'
    block = 'para_bottom_volume para_top_volume'
    #'initial linear nonlinear timestep_begin timestep_end'
  [../]
  [./Fbulk_ferro]
    type = BulkEnergy
    block = 'ferro_volume'
    polar_x = FE_polar_x
    polar_y = FE_polar_y
    polar_z = FE_polar_z
    alpha1 = ${alpha1}
    alpha11 = ${alpha11}
    alpha12 = ${alpha12}
    alpha111 = ${alpha111}
    alpha112 = ${alpha112}
    alpha123 = ${alpha123}
    execute_on = 'initial timestep_end'
  [../]
  [./Fwall_ferro]
    type = WallEnergy
    polar_x = FE_polar_x
    polar_y = FE_polar_y
    polar_z = FE_polar_z
    G110 = ${G110}
    G11_G110 = ${G11_G110}
    G12_G110 = ${G12_G110}
    G44_G110 = ${G44_G110}
    G44P_G110 = ${G44P_G110}
    
    block = 'ferro_volume'
    execute_on = 'initial timestep_end'
  [../]
  [./Felec_ferro]
    type = ElectrostaticEnergy
    block = 'ferro_volume'
    polar_x = FE_polar_x
    polar_y = FE_polar_y
    polar_z = FE_polar_z
    potential_E_int = potential_E_int
    execute_on = 'initial timestep_end'
  [../]
    
  [./Fbulk_para]
    type = BulkEnergy
    block = 'para_bottom_volume para_top_volume'
    polar_x = PE_polar_x
    polar_y = PE_polar_y
    polar_z = PE_polar_z
    alpha1 = ${alpha1_para}
    alpha11 = ${alpha11_para}
    alpha12 = ${alpha12_para}
    alpha111 = ${alpha111_para}
    alpha112 = ${alpha112_para}
    alpha123 = ${alpha123_para}
    execute_on = 'initial timestep_end'
  [../]
  [./Fwall_para]
    type = WallEnergy
    polar_x = PE_polar_x
    polar_y = PE_polar_y
    polar_z = PE_polar_z
    G110 = ${G110_para}
    G11_G110 = ${G11_G110_para}
    G12_G110 = ${G12_G110_para}
    G44_G110 = ${G44_G110_para}
    G44P_G110 = ${G44P_G110_para}
    
    block = 'para_bottom_volume para_top_volume'
    execute_on = 'initial timestep_end'
  [../]
  [./Felec_para]
    type = ElectrostaticEnergy
    block = 'para_bottom_volume para_top_volume'
    polar_x = PE_polar_x
    polar_y = PE_polar_y
    polar_z = PE_polar_z
    potential_E_int = potential_E_int
    execute_on = 'initial timestep_end'
  [../]
  [./perc_change_ferro]
    type = PercentChangePostprocessor
    postprocessor = Fbulk_ferro
  [../]
[]

[UserObjects]
  [./kill]
    type = Terminator
    expression = 'perc_change_ferro <= 2.5e-6'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol  -pc_type '     
    petsc_options_value = '    120               1e-10      1e-8     1e-5        bjacobi'    
  [../]
[]


[Executioner]
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5
    growth_factor = 1.4
    cutback_factor =  0.9
    optimal_iterations = 7
  [../]
  type = Transient
  solve_type = 'PJFNK'       #"PJFNK, JFNK, NEWTON"
  scheme = 'bdf2'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmax = 1.0
  num_steps = 25
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    execute_on = 'timestep_end'
    file_base = out_split_polar
    elemental_as_nodal = true
    interval = 1
  [../]
[]
