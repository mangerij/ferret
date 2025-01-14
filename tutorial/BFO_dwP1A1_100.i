Nx = 250
Ny = 1
Nz = 1

xMax = 31.41592653589793
yMax = 1.0
zMax = 1.0

freq = 0.2

g11 = 12e-3
g12 = -3.0e-3
g44 = 3.0e-3

h11 = 2.0e-4
h12 = -0.2e-3
h44 = 0.8e-3


[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${Nx}
    ny = ${Ny}
    nz = ${Nz}
    xmin = 0.0
    xmax = ${xMax}
    ymin = 0.0
    ymax = ${yMax}
    zmin = 0.0
    zmax = ${zMax}
    elem_type = HEX8
  []
  [./cnode]
    input = gen

    ############################################
    ##
    ##   additional boundary sideset (one node)
    ##   to zero one of the elastic displacement vectors
    ##   vectors and eliminates rigid body translations
    ##   from the degrees of freedom
    ##
    ##   NOTE: This must conform with the about
    ##         [Mesh] block settings
    ##
    ############################################

    type = ExtraNodesetGenerator
    coord = '0.0 0.0 0.0'
    new_boundary = 100
  [../]
[]

[GlobalParams]
  len_scale = 1.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  antiphase_A_x = antiphase_A_x
  antiphase_A_y = antiphase_A_y
  antiphase_A_z = antiphase_A_z

  displacements = 'u_x u_y u_z'

  potential_E_int = potential_E_int
[]


[Functions]
  [./stripeP1]
    type = ParsedFunction
    value = 0.54*cos(${freq}*(x))
  [../]
  [./stripeP2]
    type = ParsedFunction
    value = -0.54*cos(${freq}*(x))
  [../]
  [./stripeA1]
    type = ParsedFunction
    value = 7.37*cos(${freq}*(x))
  [../]
  [./stripeA2]
    type = ParsedFunction
    value = -7.37*cos(${freq}*(x))
  [../]

  [./constPm]
    type = ParsedFunction
    value = -0.54
  [../]
  [./constPp]
    type = ParsedFunction
    value = 0.54
  [../]
  [./constAm]
    type = ParsedFunction
    value = -7.37
  [../]
  [./constAp]
    type = ParsedFunction
    value = 7.37
  [../]
[]


[Variables]
  [./u_x]
  [../]
  [./u_y]
  [../]
  [./u_z]
  [../]
  [./global_strain]
    order = SIXTH
    family = SCALAR
  [../]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = constPp
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = constPp
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = stripeP1
    [../]
  [../]
  [./antiphase_A_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = constAp
    [../]
  [../]
  [./antiphase_A_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = constAp
    [../]
  [../]
  [./antiphase_A_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = stripeA1
    [../]
  [../]

  [./potential_E_int]
    order = FIRST
    family = LAGRANGE
  [../]
[]


[AuxVariables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
  [./s00]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s01]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./s11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e00]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e01]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e10]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e21]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e02]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e20]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./TensorMechanics]
  [../]

  [./rotostr_ux]
    type = RotostrictiveCouplingDispDerivative
    variable = u_x
    component = 0
  [../]
  [./rotostr_uy]
    type = RotostrictiveCouplingDispDerivative
    variable = u_y
    component = 1
  [../]
  [./rotostr_uz]
    type = RotostrictiveCouplingDispDerivative
    variable = u_z
    component = 2
  [../]

  [./electrostr_ux]
    type = ElectrostrictiveCouplingDispDerivative
    variable = u_x
    component = 0
  [../]
  [./electrostr_uy]
    type = ElectrostrictiveCouplingDispDerivative
    variable = u_y
    component = 1
  [../]
  [./electrostr_uz]
    type = ElectrostrictiveCouplingDispDerivative
    variable = u_z
    component = 2
  [../]

  ### Operators for the polar field: ###
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


  [./walled_a_x]
    type = AFDWallEnergyDerivative
    variable = antiphase_A_x
    component = 0
  [../]
  [./walled_a_y]
    type = AFDWallEnergyDerivative
    variable = antiphase_A_y
    component = 1
  [../]
  [./walled_a_z]
     type = AFDWallEnergyDerivative
     variable = antiphase_A_z
     component = 2
  [../]

  [./walled2_a_x]
    type = AFDWall2EnergyDerivative
    variable = antiphase_A_x
    component = 0
  [../]
  [./walled2_a_y]
    type = AFDWall2EnergyDerivative
    variable = antiphase_A_y
    component = 1
  [../]
  [./walled2_a_z]
     type = AFDWall2EnergyDerivative
     variable = antiphase_A_z
     component = 2
  [../]

  [./roto_polar_coupled_x]
    type = RotoPolarCoupledEnergyPolarDerivativeAlt
    variable = polar_x
    component = 0
  [../]
  [./roto_polar_coupled_y]
    type = RotoPolarCoupledEnergyPolarDerivativeAlt
    variable = polar_y
    component = 1
  [../]
  [./roto_polar_coupled_z]
    type = RotoPolarCoupledEnergyPolarDerivativeAlt
    variable = polar_z
    component = 2
  [../]
  [./roto_dis_coupled_x]
    type = RotoPolarCoupledEnergyDistortDerivativeAlt
    variable = antiphase_A_x
    component = 0
  [../]
  [./roto_dis_coupled_y]
    type = RotoPolarCoupledEnergyDistortDerivativeAlt
    variable = antiphase_A_y
    component = 1
  [../]
  [./roto_dis_coupled_z]
    type = RotoPolarCoupledEnergyDistortDerivativeAlt
    variable = antiphase_A_z
    component = 2
  [../]

  [./electrostr_polar_coupled_x]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_x
    component = 0
    u_x = disp_x
    u_y = disp_y
    u_z = disp_z
  [../]
  [./electrostr_polar_coupled_y]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_y
    component = 1
    u_x = disp_x
    u_y = disp_y
    u_z = disp_z
  [../]
  [./electrostr_polar_coupled_z]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_z
    component = 2
    u_x = disp_x
    u_y = disp_y
    u_z = disp_z
  [../]

  #Operators for the AFD field

  [./rbed_x]
    type = RotoBulkEnergyDerivativeEighthAlt
    variable = antiphase_A_x
    component = 0
  [../]
  [./rbed_y]
    type = RotoBulkEnergyDerivativeEighthAlt
    variable = antiphase_A_y
    component = 1
  [../]
  [./rbed_z]
    type = RotoBulkEnergyDerivativeEighthAlt
    variable = antiphase_A_z
    component = 2
  [../]

  [./rotostr_dis_coupled_x]
    type = RotostrictiveCouplingDistortDerivative
    variable = antiphase_A_x
    component = 0
    u_x = disp_x
    u_y = disp_y
    u_z = disp_z
  [../]
  [./rotostr_dis_coupled_y]
    type = RotostrictiveCouplingDistortDerivative
    variable = antiphase_A_y
    component = 1
    u_x = disp_x
    u_y = disp_y
    u_z = disp_z
  [../]
  [./rotostr_dis_coupled_z]
    type = RotostrictiveCouplingDistortDerivative
    variable = antiphase_A_z
    component = 2
    u_x = disp_x
    u_y = disp_y
    u_z = disp_z
  [../]

  [./polar_x_electric_E]
    type = PolarElectricEStrong
    variable = potential_E_int
  [../]
  [./FE_E_int]
    type = Electrostatics
    variable = potential_E_int
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

  [./polar_x_time]
    type = TimeDerivativeScaled
    variable=polar_x
    time_scale = 1.0
    block = '0'
  [../]
  [./polar_y_time]
    type = TimeDerivativeScaled
    variable=polar_y
    time_scale = 1.0
    block = '0'
  [../]
  [./polar_z_time]
    type = TimeDerivativeScaled
    variable = polar_z
    time_scale = 1.0
    block = '0'
  [../]

  [./a_x_time]
    type = TimeDerivativeScaled
    variable = antiphase_A_x
    time_scale = 0.01
    block = '0'
  [../]
  [./a_y_time]
    type = TimeDerivativeScaled
    variable = antiphase_A_y
    time_scale = 0.01
    block = '0'
  [../]
  [./a_z_time]
    type = TimeDerivativeScaled
    variable = antiphase_A_z
    time_scale = 0.01
    block = '0'
  [../]

  [./u_x_time]
    type = TimeDerivativeScaled
    variable = u_x
    time_scale = 1.0
  [../]
  [./u_y_time]
    type = TimeDerivativeScaled
    variable = u_y
    time_scale = 1.0
  [../]
  [./u_z_time]
    type = TimeDerivativeScaled
    variable = u_z
    time_scale = 1.0
  [../]
[]


[AuxKernels]
  [./disp_x]
    type = GlobalDisplacementAux
    variable = disp_x
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
    component = 0
  [../]
  [./disp_y]
    type = GlobalDisplacementAux
    variable = disp_y
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
    component = 1
  [../]
  [./disp_z]
    type = GlobalDisplacementAux
    variable = disp_z
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
    component = 2
  [../]
  [./s00]
    type = RankTwoAux
    variable = s00
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
  [../]
  [./s01]
    type = RankTwoAux
    variable = s01
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
  [../]
  [./s10]
    type = RankTwoAux
    variable = s10
    rank_two_tensor = stress
    index_i = 1
    index_j = 0
  [../]
  [./s11]
    type = RankTwoAux
    variable = s11
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
  [../]
  [./e00]
    type = RankTwoAux
    variable = e00
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 0
  [../]
  [./e01]
    type = RankTwoAux
    variable = e01
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 1
  [../]
  [./e10]
    type = RankTwoAux
    variable = e10
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 0
  [../]
  [./e11]
    type = RankTwoAux
    variable = e11
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 1
  [../]
  [./e12]
    type = RankTwoAux
    variable = e12
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 2
  [../]
  [./e21]
    type = RankTwoAux
    variable = e21
    rank_two_tensor = total_strain
    index_i = 2
    index_j = 1
  [../]
  [./e20]
    type = RankTwoAux
    variable = e20
    rank_two_tensor = total_strain
    index_i = 2
    index_j = 0
  [../]
  [./e02]
    type = RankTwoAux
    variable = e02
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 2
  [../]
  [./e22]
    type = RankTwoAux
    variable = e22
    rank_two_tensor = total_strain
    index_i = 2
    index_j = 2
  [../]
[]

[ScalarKernels]
  [./global_strain]
    type = GlobalStrain
    variable = global_strain
    global_strain_uo = global_strain_uo
  [../]
[]

[Materials]
  [./Landau_P]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '-2.81296 1.72351 2.24147 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
  [../]

  [./Landau_A]
    type = GenericConstantMaterial
    prop_names = 'beta1 beta11 beta12 beta111 beta112 beta123 beta1111 beta1112 beta1122 beta1123'
    prop_values = '-0.0137763 0.0000349266 0.0000498846 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
  [../]

  [./P_A_couple]
    type = GenericConstantMaterial
    prop_names = 't1111 t1122 t1212 t42111111 t24111111 t42111122 t24112222 t42112233 t24112233 t42112211 t24111122 t42111212   t42123312 t24121112 t24121233 t6211111111 t2611111111 t6211111122 t2611222222 t4411111111 t4411112222'
    prop_values = '0.012516 0.0180504 -0.036155 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
  [../]

  [./Landau_G]
    type = GenericConstantMaterial
    prop_names = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '1.0 ${g11} ${g12} ${g44} 0.0'
  [../]
  [./Landau_H]
    type = GenericConstantMaterial
    prop_names = 'H110 H11_H110 H12_H110 H44_H110 H44P_H110'
    prop_values = '1.0 ${h11} ${h12} ${h44} 0.0'
  [../]

  [./mat_C]
    type = GenericConstantMaterial
    prop_names = 'C11 C12 C44'
    prop_values = '295.179 117.567 74.0701'
  [../]
  [./mat_Q]
    type = GenericConstantMaterial
    prop_names = 'Q11 Q12 Q44'
    prop_values = '-0.0603833 0.0111245 -0.0175686'
  [../]
  [./mat_R]
    type = GenericConstantMaterial
    prop_names = 'R11 R12 R44'
    prop_values = '-0.0000878064 0.0000295306 0.0000627962'
  [../]
  [./mat_q]
    type = GenericConstantMaterial
    prop_names = 'q11 q12 q44'
    prop_values = '-30.4162 -5.01496 -10.4105'
  [../]
  [./mat_r]
    type = GenericConstantMaterial
    prop_names = 'r11 r12 r44'
    prop_values = '-0.0379499 0.00373096 0.0372105'
  [../]
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '295.179 117.567 117.567 295.179 117.567 295.179 74.0701 74.0701 74.0701'
  [../]

  [./strain]
    type = ComputeSmallStrain
    global_strain = global_strain
  [../]
  [./global_strain]
    type = ComputeGlobalStrain
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
  [../]

  [./stress]
    type = ComputeLinearElasticStress
  [../]

  [./permitivitty_1]

    ###############################################
    ##
    ##  so-called background dielectric constant
    ##  (it encapsulates the motion of core electrons
    ##  at high frequency) = e_b*e_0 (here we use
    ##  e_b = 10), see PRB. 74, 104014, (2006)
    ##
    ###############################################

    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '0.08854187'
  [../]

[]

[Postprocessors]
  [./dt]
     type = TimestepSize
  [../]
  [./FbP]
    type = BulkEnergyEighth
    execute_on = 'timestep_end'
  [../]
  [./FbA]
    type = RotoBulkEnergyEighth
    execute_on = 'timestep_end'
  [../]
  [./FcPA]
    type = RotoPolarCoupledEnergyEighth
    execute_on = 'timestep_end'
  [../]
  [./FgP]
    type = WallEnergy
    execute_on = 'timestep_end'
  [../]
  [./FgA]
    type = AFDWallEnergy
    execute_on = 'timestep_end'
  [../]
  [./FcPu]
    type = ElectrostrictiveCouplingEnergy
    execute_on = 'timestep_end'
    u_x = disp_x
    u_y = disp_y
    u_z = disp_z
  [../]
  [./FcAu]
    type = RotostrictiveCouplingEnergy
    execute_on = 'timestep_end'
    u_x = disp_x
    u_y = disp_y
    u_z = disp_z
  [../]

  [./Felu]
    type = ElasticEnergy
    execute_on = 'timestep_end'
  [../]
  [./Fele]
    type = ElectrostaticEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'FbP FbA FgP FgA FcPA FcPu FcAu Felu Fele'
    pp_coefs = ' 1 1 1 1 1 1 1 1 1'
    execute_on = 'timestep_end'

    ##########################################
    #
    # NOTE: Ferret output is in attojoules
    #
    ##########################################
  [../]
  [./perc_change]
    type = EnergyRatePostprocessor
    postprocessor = Ftot
    execute_on = 'timestep_end'
    dt = dt
  [../]

  [./elapsed]
    type = PerfGraphData
    section_name = "Root"  # for profiling the problem
    data_type = total
  [../]
[]


[BCs]
  [./Periodic]
    [./x]
      auto_direction = 'x'
      variable = 'u_x u_y u_z polar_x polar_y polar_z antiphase_A_x antiphase_A_y antiphase_A_z'
    [../]
    [./xyz]
      auto_direction = 'x y z'
      variable = 'potential_E_int'
    [../]
  [../]
  # fix center point location
  [./centerfix_x]
    type = DirichletBC
    boundary = 100
    variable = u_x
    value = 0
  [../]
  [./centerfix_y]
    type = DirichletBC
    boundary = 100
    variable = u_y
    value = 0
  [../]
  [./centerfix_z]
    type = DirichletBC
    boundary = 100
    variable = u_z
    value = 0
  [../]
[]


[UserObjects]
  [./global_strain_uo]
    type = GlobalBFOMaterialRVEUserObject
    execute_on = 'Initial Linear Nonlinear'
  [../]
  [./kill]
   type = Terminator
   expression = 'perc_change <= 5.0e-7'
  [../]
[]

#=

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -build_twosided'
    petsc_options_value = '    121            1e-8          1e-7       1e-6     bjacobi    allreduce'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  scheme = 'bdf2'
  dtmin = 1e-13
  dtmax = 10.0

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 25  #usually 10
    linear_iteration_ratio = 100
    dt = 0.08
    growth_factor = 1.1
  [../]

  num_steps = 1000
[]

[Outputs]
  print_linear_residuals = false
  perf_graph_live = false
  [./out]
    type = Exodus
    file_base = BFO_dwP1A1_100
    elemental_as_nodal = true
  [../]
[]
