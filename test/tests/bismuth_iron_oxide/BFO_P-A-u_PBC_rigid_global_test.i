[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 4
  ny = 4
  nz = 4
  xmin = -1.0
  xmax = 1.0
  ymin = -1.0
  ymax = 1.0
  zmin = -1.0
  zmax = 1.0
  elem_type = HEX8
[]

[MeshModifiers]
  [./cnode]
    type = AddExtraNodeset
    coord = '0.0 0.0 0.0'
    new_boundary = 100
  [../]
[]

[GlobalParams]
  len_scale = 1.0

  alpha1 = -3.29328296e-1
  alpha11 = 7.15515014e-2
  alpha12 = 8.85202744e-2
  alpha111 = -4.44849403e-3
  alpha112 = 1.96626638e-3
  alpha123 = -5.32358548e-2
  alpha1111 = 1.86312506e-4
  alpha1112 = -4.70479587e-4
  alpha1122 = 9.55233783e-4
  alpha1123 = 3.0759255e-3

  beta1 = -3.88826748e-1
  beta11 = 9.14933665e-2
  beta12 = 1.07046279e-1
  beta111 = -7.73266661e-3
  beta112 = -6.09716685e-3
  beta123 = -6.926e-3
  beta1111 = 3.961e-4
  beta1112 = 1.294e-4
  beta1122 = 9.67e-4
  beta1123 = 8.115e-4

  t1111 = 1.165e-1
  t1122 = 1.539e-1
  t1212 = -1.925e-1
  t42111111 = 4.432e-3
  t24111111 = -2.662e-2
  t42111122 = -1.695e-2
  t24112222 = -2.157e-2
  t42112233 = -1.577e-2
  t24112233 = -1.133e-2
  t42112211 = 1.272e-2
  t24111122 = 1.214e-2
  t42111212 = -3.652e-2
  t42123312 = 4.972e-2
  t24121112 = -3.891e-3
  t24121233 = -2.554e-2
  t6211111111 = -1.327e-3
  t2611111111 = 3.663e-3
  t6211111122 = 7.066e-4
  t2611222222 = 1.7e-3
  t4411111111 = 6.133e-3
  t4411112222 = 2.438e-3

  G110 = 0.173
  G11_G110 = 0.5
  G12_G110 = 0
  G44_G110 = 0.5 
  G44P_G110 = 0.5

  H110 = 0.173
  H11_H110 = 0.5
  H12_H110 = 0
  H44_H110 = 0.5
  H44P_H110 = 0.5

  #from mathematica

  C11 = 1.762732e-2
  C12 = 9.4905087e-3
  C44 = 5.24373333e-5

  Q11 = -1.35386
  Q12 = 0.238295
  Q44 = -4.17474
  R11 = -2.35149
  R12 = 1.1143805
  R44 = 3.67475

  q11 = -1.93418062e-2
  q12 = -6.3875442e-3
  q44 = -8.75649995e-4

  r11 = -7.7644e-3
  r12 = 1.63357e-3 
  r44 = 7.70775329e-4

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  antiferrodis_A_x = antiferrodis_A_x
  antiferrodis_A_y = antiferrodis_A_y
  antiferrodis_A_z = antiferrodis_A_z

  displacements = 'u_x u_y u_z'
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
      min = 0.08
      max = 0.09
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.08
      max = 0.09
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.08
      max = 0.09
    [../]
  [../]
  [./antiferrodis_A_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.08
      max = 0.09
    [../]
  [../]
  [./antiferrodis_A_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.08
      max = 0.09
    [../]
  [../]
  [./antiferrodis_A_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.08
      max = 0.09
    [../]
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

  [./electrostr_polar_coupled_x]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_x
    component = 0
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
  [./electrostr_polar_coupled_y]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_y
    component = 1
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
  [./electrostr_polar_coupled_z]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_z
    component = 2
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
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

  [./roto_walled_x]
    type = AFDAntiphaseEnergyDerivative
    variable = antiferrodis_A_x
    component = 0
  [../]
  [./roto_walled_y]
    type = AFDAntiphaseEnergyDerivative
    variable = antiferrodis_A_y
    component = 1
  [../]
  [./roto_walled_z]
    type = AFDAntiphaseEnergyDerivative
    variable = antiferrodis_A_z
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

  [./rotostr_dis_coupled_x]
    type = RotostrictiveCouplingDistortDerivative
    variable = antiferrodis_A_x
    component = 0
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
  [./rotostr_dis_coupled_y]
    type = RotostrictiveCouplingDistortDerivative
    variable = antiferrodis_A_y
    component = 1
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
  [./rotostr_dis_coupled_z]
    type = RotostrictiveCouplingDistortDerivative
    variable = antiferrodis_A_z
    component = 2
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]

  ###Time dependence
  [./polar_x_time]
    type = TimeDerivativeScaled
    variable = polar_x
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

  [./antiferrodis_A_x_time]
    type = TimeDerivativeScaled
    variable = antiferrodis_A_x
    time_scale = 1.0
  [../]
  [./antiferrodis_A_y_time]
     type = TimeDerivativeScaled
     variable = antiferrodis_A_y
    time_scale = 1.0
  [../]
  [./antiferrodis_A_z_time]
     type = TimeDerivativeScaled
     variable = antiferrodis_A_z
    time_scale = 1.0
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
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '1.762732e-2 9.4905087e-3 9.4905087e-3 1.762732e-2 9.4905087e-3 1.762732e-2 5.24373333e-5 5.24373333e-5 5.24373333e-5'
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
  [./FgP]
    type = WallEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./FgA]
    type = AFDWallEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./FcPA]
    type = RotoPolarCoupledEnergyEighth
    execute_on = 'initial timestep_end'
  [../]
  [./FcPu]
    type = ElectrostrictiveCouplingEnergy
    execute_on = 'initial timestep_end'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
  [./FcAu]
    type = RotostrictiveCouplingEnergy
    execute_on = 'initial timestep_end'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
  [./Felu]
    type = ElasticEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor 
    pp_names = 'FbP FbA FgP FgA FcPA FcPu FcAu Felu' 
    pp_coefs = ' 1 1 1 1 1 1 1 1' 
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
      variable = 'u_x u_y u_z polar_x polar_y polar_z antiferrodis_A_x antiferrodis_A_y antiferrodis_A_z'
    [../]
  [../]
  # fix center point location
  [./centerfix_x]
    type = PresetBC
    boundary = 100
    variable = u_x
    value = 0
  [../]
  [./centerfix_y]
    type = PresetBC
    boundary = 100
    variable = u_y
    value = 0
  [../]
  [./centerfix_z]
    type = PresetBC
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
   expression = 'perc_change <= 1.0e-7'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121            1e-10          1e-8       1e-8     bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.18
  solve_type = 'NEWTON'
  scheme = 'bdf2'
  dtmin = 1e-13
  dtmax = 0.18
  num_steps = 4
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_global_test
    elemental_as_nodal = true
  [../]
[]
