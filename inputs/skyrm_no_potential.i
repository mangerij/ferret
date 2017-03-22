
[Mesh]
  file = exodus_disk_r8_h1.e
[]

[MeshModifiers]
  [./centernodeset_1]
    type = AddExtraNodeset
    new_boundary = 'center_node_1'
    coord = '0.0 0.0 -0.5'
  [../]
  [./centernodeset_2]
    type = AddExtraNodeset
    new_boundary = 'center_node_2'
    coord = '0.0 0.0 0.5'
  [../]
  [./centernodeset_3]
    type = AddExtraNodeset
    new_boundary = 'center_node_3'
    coord = '0.0 0.0 0.166667'
  [../]
  [./centernodeset_4]
    type = AddExtraNodeset
    new_boundary = 'center_node_4'
    coord = '0.0 0.0 -0.166667'
  [../]
[]

[GlobalParams]
  len_scale = 1.0

  alpha1 = -0.09179 #room temp PTO
  alpha11 = 0.0706
  alpha111 = 0.0
  alpha12 = 0.1412
  alpha112 = 0.0
  alpha123 = 0.0

  G110 = 0.141
  G11/G110 = 0.0 #this is here to somehow prevent P_z "ringing" problems on the side...
  G12/G110 = 0.0 #perhaps this allows for divP =0?
  G44/G110 = 1.0
  G44P/G110 = 1.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  potential_int = potential_int
[]

[Functions]
  [./parsed_function_x]
    type = ParsedFunction
    value = '-(0.738217-0.00686984*(x^2+y^2)^(0.5)+0.00644497*(x^2+y^2)-0.0188174*(x^2+y^2)^(1.5)+0.00441745*(x^2+y^2)^2-0.000274842*(x^2+y^2)^(5/2))*sin(-0.028395+0.267482*(x^2+y^2)^(0.5)-0.146762*(x^2+y^2)+0.0632932*(x^2+y^2)^(1.5)-0.00790942*(x^2+y^2)^(2)+0.000294936*(x^2+y^2)^(5/2))*sin(atan(y/x))'
  [../]
  [./parsed_function_y]
    type = ParsedFunction
    value = '(0.738217-0.00686984*(x^2+y^2)^(0.5)+0.00644497*(x^2+y^2)-0.0188174*(x^2+y^2)^(1.5)+0.00441745*(x^2+y^2)^2-0.000274842*(x^2+y^2)^(5/2))*sin(-0.028395+0.267482*(x^2+y^2)^(0.5)-0.146762*(x^2+y^2)+0.0632932*(x^2+y^2)^(1.5)-0.00790942*(x^2+y^2)^(2)+0.000294936*(x^2+y^2)^(5/2))*cos(atan(y/x))'
  [../]
  [./parsed_function_z]
    type = ParsedFunction
    value = '(0.738217-0.00686984*(x^2+y^2)^(0.5)+0.00644497*(x^2+y^2)-0.0188174*(x^2+y^2)^(1.5)+0.00441745*(x^2+y^2)^2-0.000274842*(x^2+y^2)^(5/2))*cos(-0.028395+0.267482*(x^2+y^2)^(0.5)-0.146762*(x^2+y^2)+0.0632932*(x^2+y^2)^(1.5)-0.00790942*(x^2+y^2)^(2)+0.000294936*(x^2+y^2)^(5/2))'
  [../]
[]

[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = parsed_function_x
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = parsed_function_y
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = parsed_function_z
    [../]
  [../]
[]


[Kernels]
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


  [./anis_x]
    type = AnisotropyEnergy
    variable = polar_x
    component = 0
    K = 0.0565 #minus means IN-PLANE
  [../]
  [./anis_y]
    type = AnisotropyEnergy
    variable = polar_y
    component = 1
    K = 0.0565
  [../]

  ##Electrostatics

  #PdotE
  #[./polar_electric_px]
  #   type = PolarElectricPStrong
  #   variable = polar_x
  #   component = 0
  #[../]
  #[./polar_electric_py]
  #   type = PolarElectricPStrong
  #   variable = polar_y
  #   component = 1
  #[../]
  #[./polar_electric_pz] #  when K > 0, turning this block off will get the two band state
  #   type = PolarElectricPStrong
  #   variable = polar_z
  #   component = 2
  #[../]

  #Poisson problem:
  #[./FE_E_int]
  #   type = Electrostatics
  #   variable = potential_int
  #   block = '1'
  #   permittivity = 3.54
  #[../]
  #[./polar_x_electric_E]
  #   type = PolarElectricEStrong
  #   variable = potential_int
  #[../]

  [./depol_z]
    type = DepolEnergy
    permitivitty = 0.008854187
    lambda = 0.00175 #this needs to be adjusted as a function of thickness. Can we remove the ringing with G11 = 0?
    variable = polar_z
    avePz = avePz
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


[BCs]
  #[./center_pol_x]
  #  type = DirichletBC
  #  variable = 'polar_x'
  #  value = 0.0
  #  boundary = 'center_node_1 center_node_2 center_node_3 center_node_4'
  #[../]
  #[./center_pol_y]
  #  type = DirichletBC
  #  variable = 'polar_y'
  #  value = 0.0
  #  boundary = 'center_node_1 center_node_2 center_node_3 center_node_4'
  #[../]

 # [./center_pol_z]
 #   type = DirichletBC
 #   variable = 'polar_z'
 #   value = 0.65
 #   boundary = 'center_node_1 center_node_2 center_node_3 center_node_4'
 # [../]

  #[./side_neumann_x]
  #  variable = 'polar_x'
  #  type = NeumannBC
  #  value = 0.0
  #  boundary = '1'
  #[../]
  #[./side_neumann_y]
  #  variable = 'polar_y'
  #  type = NeumannBC
  #  value = 0.0
  #  boundary = '1'
  #[../]
  #[./side_Neumann_z]
  #  variable = 'polar_z'
  #  type = NeumannBC
  #  value = 0.0
  #  boundary = '1'
  #[../]
[]



[Postprocessors]
   [./avePz]
     type = ElementAverageValue
     variable = polar_z
     execute_on = 'initial linear nonlinear timestep_begin timestep_end'
   [../]
   [./Fbulk]
      type = BulkEnergy
      execute_on = 'timestep_end'
    [../]
    [./Fwall]
      type = WallEnergy
      execute_on = 'timestep_end'
    [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    250              1e-10      1e-8      1e-6      bjacobi   '
  [../]
[]

[Executioner]
  type = Transient
    [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    #iteration_window = 3
    optimal_iterations = 6 #should be 5 probably
    growth_factor = 1.4
    linear_iteration_ratio = 1000
    cutback_factor =  0.9
[../]
  solve_type = 'PJFNK'       #"PJFNK, JFNK, NEWTON"
  scheme = 'implicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin = 1e-13
  dtmax = 0.7
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    execute_on = 'timestep_end'
    file_base = out_skyrm_L01_no_potential_noG11
    elemental_as_nodal = true
  [../]
[]
