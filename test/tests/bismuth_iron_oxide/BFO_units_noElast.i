[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 3
    ny = 3
    nz = 3
    xmin = -0.5
    xmax = 0.5
    ymin = -0.5
    ymax = 0.5
    zmin = -0.5
    zmax = 0.5
    elem_type = HEX8
  []
  [./cnode]
    input = gen
    type = ExtraNodesetGenerator
    coord = '-0.5 -0.5 -0.5'
    new_boundary = 100
  [../]
[]

[GlobalParams]
  len_scale = 1.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  antiferrodis_A_x = antiferrodis_A_x
  antiferrodis_A_y = antiferrodis_A_y
  antiferrodis_A_z = antiferrodis_A_z


[]

[Variables]

  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.5
      max = 0.51
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.5
      max = 0.51
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.5
      max = 0.51
    [../]
  [../]
  [./antiferrodis_A_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 7.4
      max = 7.41
    [../]
  [../]
  [./antiferrodis_A_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 7.4
      max = 7.41
    [../]
  [../]
  [./antiferrodis_A_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 7.4
      max = 7.41
    [../]
  [../]
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
    variable = antiferrodis_A_x
    component = 0
  [../]
  [./roto_dis_coupled_y]
    type = RotoPolarCoupledEnergyDistortDerivativeAlt
    variable = antiferrodis_A_y
    component = 1
  [../]
  [./roto_dis_coupled_z]
    type = RotoPolarCoupledEnergyDistortDerivativeAlt
    variable = antiferrodis_A_z
    component = 2
  [../]






  #Operators for the AFD field

  [./rbed_x]
    type = RotoBulkEnergyDerivativeEighthAlt
    variable = antiferrodis_A_x
    component = 0
  [../]
  [./rbed_y]
    type = RotoBulkEnergyDerivativeEighthAlt
    variable = antiferrodis_A_y
    component = 1
  [../]
  [./rbed_z]
    type = RotoBulkEnergyDerivativeEighthAlt
    variable = antiferrodis_A_z
    component = 2
  [../]


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

  [./a_x_time]
     type = TimeDerivativeScaled
     variable = antiferrodis_A_x
     time_scale = 1.0
  [../]
  [./a_y_time]
     type = TimeDerivativeScaled
     variable = antiferrodis_A_y
     time_scale = 1.0
  [../]
  [./a_z_time]
     type = TimeDerivativeScaled
     variable = antiferrodis_A_z
     time_scale = 1.0
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

#some strange divergence happens here....

[]

[Postprocessors]
  [./dt]
     type = TimestepSize
  [../]
  [./FbP]
    type = BulkEnergyEighth
    execute_on = 'initial timestep_end'
  [../]
  [./FbA]
    type = RotoBulkEnergyEighth
    execute_on = 'initial timestep_end'
  [../]


  [./Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'FbP FbA'
    pp_coefs = ' 1 1 '
    execute_on = 'initial timestep_end'
  [../]
  [./perc_change]
    type = EnergyRatePostprocessor
    postprocessor = Ftot
    execute_on = 'initial timestep_end'
    dt = dt
  [../]
[]


[BCs]
  [./Periodic]
    [./xy]
      auto_direction = 'x y z'
      variable = 'polar_x polar_y polar_z antiferrodis_A_x antiferrodis_A_y antiferrodis_A_z'
    [../]
  [../]
[]

[UserObjects]
  [./kill]
   type = Terminator
   expression = 'perc_change <= 1.0e-8'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121            1e-10          1e-10       1e-8    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.05
  solve_type = 'NEWTON'
  scheme = 'bdf2'
  dtmin = 1e-13
  dtmax = 10.0
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_bfo
    elemental_as_nodal = true
  [../]
[]
