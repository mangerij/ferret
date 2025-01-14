Nx = 3
Ny = 3
Nz = 3

xMax = 1.0
yMax = 1.0
zMax = 1.0


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


[]


[Functions]
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
      function = constPp
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
      function = constAp
    [../]
  [../]
[]

[AuxVariables]

[]

[AuxKernels]

[]

[Kernels]
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

  [./Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'FbP FbA FcPA'
    pp_coefs = ' 1 1 1'
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

  [./nodes]
    type = NumNodes
  [../]
[]


[BCs]
  [./Periodic]
    [./xy]
      auto_direction = 'x y z'
      variable = 'polar_x polar_y polar_z antiphase_A_x antiphase_A_y antiphase_A_z'
    [../]
  [../]
[]

[UserObjects]
  [./kill]
   type = Terminator
   expression = 'perc_change <= 5.0e-5'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type -build_twosided'
    petsc_options_value = '    121            1e-10          1e-10       1e-6     bjacobi    allreduce'
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.08
  solve_type = 'NEWTON'
  scheme = 'bdf2'
  dtmin = 1e-13
  dtmax = 10.0

  num_steps = 10
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_P0A0_f
    elemental_as_nodal = true
  [../]
[]
