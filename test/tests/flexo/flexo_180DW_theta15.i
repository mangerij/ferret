g11 = 0.51
g12 = -0.02
g44 = 0.02

L = 0.25
n_xy = 1
Z = 11.0
n_z = 88.0

period = 0.28559933214452665

thetadeg = 15
theta0 = ${fparse thetadeg*pi/180}

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${n_z}
    ny = ${n_xy}
    nz = ${n_xy}
    xmin = -${Z}
    xmax = ${Z}
    ymin = -${L}
    ymax = ${L}
    zmin = -${L}
    zmax = ${L}
    elem_type = HEX8
  []
  [./cnode]
    input = gen
    type = ExtraNodesetGenerator
    coord = '-${Z} -${L} -${L}'
    new_boundary = 100
  [../]

[]

[GlobalParams]

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  displacements = 'u_x u_y u_z'

  potential_E_int = potential_E_int
[]

[Functions]
  [./DW_func]
    type = ParsedFunction
    symbol_names = 'A p'
    symbol_values = '0.2 ${period}'
    expression = 'A*cos(p*x)'
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
      type = RandomIC
      min = -1e-5
      max = 1e-5
      seed = 1
    [../]
    block = '0'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -1e-5
      max = 1e-5
      seed = 2
    [../]
    block = '0'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = DW_func
    [../]
    block = '0'
  [../]
  [./potential_E_int]
    order = FIRST
    family = LAGRANGE
    block = '0'
  [../]
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [./all_FE]
        strain = SMALL
        add_variables = true
        incremental = false
        eigenstrain_names = 'total_eigenstrain'
        global_strain = global_strain
        generate_output = 'strain_xx strain_yy strain_zz strain_yz strain_xz strain_xy'
	block = '0'
      [../]
    []

    [./GlobalStrain]
      [./global_strain]
        scalar_global_strain = global_strain
        displacements = 'u_x u_y u_z'
        auxiliary_displacements = 'disp_x disp_y disp_z'
        global_displacements = 'ug_x ug_y ug_z'
      [../]
    [../]
  []
[]

[Kernels]

  [./bed_x]
    type = RotatedBulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
    block = '0'
  [../]
  [./bed_y]
    type = RotatedBulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
    block = '0'
  [../]
  [./bed_z]
    type = RotatedBulkEnergyDerivativeEighth
    variable = polar_z
    component = 2
    block = '0'
  [../]

  [./walled_x]
    type = RotatedWallEnergyDerivative
    variable = polar_x
    component = 0
    block = '0'
  [../]
  [./walled_y]
    type = RotatedWallEnergyDerivative
    variable = polar_y
    component = 1
    block = '0'
  [../]
  [./walled_z]
    type = RotatedWallEnergyDerivative
    variable = polar_z
    component = 2
    block = '0'
  [../]

  [./flexo_polar_coupled_x]
    type = RotatedCubicParentFlexoelectricFullPDerivative
    variable = polar_x
    component = 0
    block = '0'
  [../]
  [./flexo_polar_coupled_y]
    type = RotatedCubicParentFlexoelectricFullPDerivative
    variable = polar_y
    component = 1
    block = '0'
  [../]
  [./flexo_polar_coupled_z]
    type = RotatedCubicParentFlexoelectricFullPDerivative
    variable = polar_z
    component = 2
    block = '0'
  [../]

  [./polar_x_electric_E]
    type = PolarElectricEStrong
    variable = potential_E_int
    block = '0'
  [../]
  [./FE_E_int]
    type = Electrostatics
    variable = potential_E_int
    block = '0'
  [../]

  [./polar_electric_px]
    type = PolarElectricPStrong
    variable = polar_x
    component = 0
    block = '0'
  [../]
  [./polar_electric_py]
    type = PolarElectricPStrong
    variable = polar_y
    component = 1
    block = '0'
  [../]
  [./polar_electric_pz]
    type = PolarElectricPStrong
    variable = polar_z
    component = 2
    block = '0'
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
[]

[Materials]

  [./twist_angle_0]
    type = GenericConstantMaterial
    prop_names = 'twist'
    prop_values = '${theta0}'
    block = '0'
  [../]

  [./eigen_strain]
    type = ComputeEigenstrain
    eigen_base = '0.0 0.0 0.0 0.0 0.0 0.0'
    eigenstrain_name = 'epitaxy'
    block = '0'
  [../]

  [./Landau_P]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '-0.04277 -0.11540 0.65300 -2.10600 4.09100 -6.68800 75.90000 -21.93000 -22.21000 24.16000'
    block = '0'
  [../]

  [./Landau0_G_FE]
    type = GenericConstantMaterial
    prop_names = 'G11 G12 G44'
    prop_values = '${g11} ${g12} ${g44}'
    block = '0'
  [../]

  [./mat_C]
    type = GenericConstantMaterial
    prop_names = 'C11 C12 C44'
    prop_values = '178.0 96.4 122.0'
    block = '0'
  [../]

  [./mat_Q]
    type = GenericConstantMaterial
    prop_names = 'Q11 Q12 Q44'
    prop_values = '0.11 -0.045 0.0295'
    block = '0'
  [../]

  [./mat_f]
    type = GenericConstantMaterial
    prop_names = 'mu11 mu12 mu44'
    prop_values = '0.012822 -0.467312 -0.162872'
    block = '0'
  [../]

  [./mat_S]
    type = GenericConstantMaterial
    prop_names = 's11 s12 s44'
    prop_values = '0.00906889 -0.003186 0.00819672'
    block = '0'
  [../]

  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '178.0 96.4 96.4 178.0 96.4 178.0 122.0 122.0 122.0'
    block = '0'
    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = ${fparse -1*thetadeg}
  [../]

  [./film_eigenstrain]
    type = CompositeEigenstrain
    tensors = 'ferro flexo epitaxy'
    weights = 'weight1 weight2 weight3'
    eigenstrain_name = total_eigenstrain
    coupled_variables = 'polar_x polar_y polar_z'
    block = '0'
  [../]
  [./weight1]
    type = DerivativeParsedMaterial
    block = '0'
    expression = '1'
    property_name = weight1
    coupled_variables = 'polar_x polar_y polar_z'
  [../]
  [./weight2]
    type = DerivativeParsedMaterial
    block = '0'
    expression = '1'
    property_name = weight2
  [../]
  [./weight3]
    type = DerivativeParsedMaterial
    block = '0'
    expression = '1'
    property_name = weight3
  [../]

  [./stress_1]
    type = ComputeLinearElasticStress
    block = '0'
  [../]

  [./electrostrictive_eigenstrain]
    type = ComputeRotatedCubicParentElectrostrictiveStrain
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
    eigenstrain_name = 'ferro'
    block = '0'
  [../]
  [./flexoelectric_eigenstrain]
    type = ComputeRotatedCubicParentFlexoelectricStrain
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
    eigenstrain_name = 'flexo'
    block = '0'
  [../]
  [./permitivitty]

    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '0.39825'
    block = '0'
  [../]
[]

[Postprocessors]
  [./dt]
     type = TimestepSize
  [../]
  [./FbP]
    type = RotatedBulkEnergyEighth
    execute_on = 'timestep_end'
    block = '0'
  [../]
  [./FgP]
    type = RotatedWallEnergy
    execute_on = 'timestep_end'
    block = '0'
  [../]

  [./Fele]
    type = ElectrostaticEnergy
    execute_on = 'timestep_end'
    block = '0'
  [../]

  [./Fela]
    type = RotatedCubicParentFlexoelectricFullEnergy
    execute_on = 'timestep_end'
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
    block = '0'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'FbP FgP Fela Fele'
    pp_coefs = ' 1 1 1 1'
    execute_on = 'timestep_end'

  [../]

  [./maxPy]
     type = NodalExtremeValue
     variable = polar_y
  [../]
  [./minPx]
     type = NodalExtremeValue
     variable = polar_x
     value_type = min
  [../]
  [./minPy]
     type = NodalExtremeValue
     variable = polar_y
     value_type = min
  [../]
  [./maxPx]
     type = NodalExtremeValue
     variable = polar_x
  [../]

  [./perc_change]
    type = EnergyRatePostprocessor
    postprocessor = Ftot
    execute_on = 'timestep_end'
    dt = dt
  [../]
[]

[BCs]
  [./Periodic]
    [./x]
      auto_direction = 'x y z'
      variable = 'u_x u_y u_z polar_x polar_y polar_z potential_E_int'
    [../]
  [../]
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
  [./kill]
   type = Terminator
   expression = 'perc_change <= -1.0'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121            1e-8          1e-6      1e-5    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  scheme = 'bdf2'
  num_steps = 12
  dtmin = 1e-13
  dtmax = 10.0

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 25
    linear_iteration_ratio = 100
    dt = 0.1
    growth_factor = 1.1
  [../]
[]

[Outputs]
  print_linear_residuals = false
  perf_graph_live = false
  [./csv]
    type = CSV
    file_base = flexo_180DW_theta15_out
  [../]
[]
