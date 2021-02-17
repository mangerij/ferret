[Mesh]
  [gen]
    type = GeneratedMeshGenerator
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
  [./cnode]
    input = gen
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

  antiferrodis_A_x = antiferrodis_A_x
  antiferrodis_A_y = antiferrodis_A_y
  antiferrodis_A_z = antiferrodis_A_z

  displacements = 'u_x u_y u_z'

  u_x = u_x
  u_y = u_y
  u_z = u_z
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
      legacy_generator = true
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.08
      max = 0.09
      legacy_generator = true
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.08
      max = 0.09
      legacy_generator = true
    [../]
  [../]
  [./antiferrodis_A_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.08
      max = 0.09
      legacy_generator = true
    [../]
  [../]
  [./antiferrodis_A_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.08
      max = 0.09
      legacy_generator = true
    [../]
  [../]
  [./antiferrodis_A_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.08
      max = 0.09
      legacy_generator = true
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
  [../]
  [./electrostr_polar_coupled_y]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_y
    component = 1
  [../]
  [./electrostr_polar_coupled_z]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_z
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

  [./roto_walled_x]
    type = WallEnergyDerivative
    variable = antiferrodis_A_x
    component = 0
    polar_x = antiferrodis_A_x
    polar_y = antiferrodis_A_y
    polar_z = antiferrodis_A_z
  [../]
  [./roto_walled_y]
    type = WallEnergyDerivative
    variable = antiferrodis_A_y
    component = 1
    polar_x = antiferrodis_A_x
    polar_y = antiferrodis_A_y
    polar_z = antiferrodis_A_z
  [../]
  [./roto_walled_z]
    type = WallEnergyDerivative
    variable = antiferrodis_A_z
    component = 2
    polar_x = antiferrodis_A_x
    polar_y = antiferrodis_A_y
    polar_z = antiferrodis_A_z
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
  [../]
  [./rotostr_dis_coupled_y]
    type = RotostrictiveCouplingDistortDerivative
    variable = antiferrodis_A_y
    component = 1
  [../]
  [./rotostr_dis_coupled_z]
    type = RotostrictiveCouplingDistortDerivative
    variable = antiferrodis_A_z
    component = 2
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
  [./Landau_P]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '-3.29328296e-1 7.15515014e-2 8.85202744e-2 -4.44849403e-3 1.96626638e-3 -5.32358548e-2 1.86312506e-4 -4.70479587e-4 9.55233783e-4 3.0759255e-3'
  [../]

  [./Landau_A]
    type = GenericConstantMaterial
    prop_names = 'beta1 beta11 beta12 beta111 beta112 beta123 beta1111 beta1112 beta1122 beta1123'
    prop_values = '-3.88826748e-1 9.14933665e-2 1.07046279e-1 -7.73266661e-3 -6.09716685e-3 -6.926e-3 3.961e-4 1.294e-4 9.67e-4 8.115e-4'
  [../]

  [./P_A_couple]
    type = GenericConstantMaterial
    prop_names = 't1111 t1122 t1212 t42111111 t24111111 t42111122 t24112222 t42112233 t24112233 t42112211 t24111122 t42111212 t42123312 t24121112 t24121233 t6211111111 t2611111111 t6211111122 t2611222222 t4411111111 t4411112222'
    prop_values = '1.165e-1 1.539e-1 -1.925e-1 4.432e-3 -2.662e-2 -1.695e-2 -2.157e-2 -1.577e-2 -1.133e-2 1.272e-2 1.214e-2 -3.652e-2 4.972e-2 -3.891e-3 -2.554e-2 -1.327e-3 3.663e-3 7.066e-4 1.7e-3 6.133e-3 2.438e-3'
  [../]

  [./grad_P]
    type = GenericConstantMaterial
    prop_names = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '0.173 0.5 0 0.5 0.5'
  [../]

  [./grad_A]
    type = GenericConstantMaterial
    prop_names = 'H110 H11_H110 H12_H110 H44_H110 H44P_H110'
    prop_values = '0.173 0.5 0 0.5 0.5'
  [../]

  [./mat_C]
    type = GenericConstantMaterial
    prop_names = 'C11 C12 C44'
    prop_values = '1.762732e-2 9.4905087e-3 5.24373333e-5'
  [../]

  [./mat_Q]
    type = GenericConstantMaterial
    prop_names = 'Q11 Q12 Q44'
    prop_values = '-1.35386 0.238295 -4.17474'
  [../]

  [./mat_R]
    type = GenericConstantMaterial
    prop_names = 'R11 R12 R44'
    prop_values = '-2.35149 1.1143805 3.67475'
  [../]

  [./mat_q]
    type = GenericConstantMaterial
    prop_names = 'q11 q12 q44'
    prop_values = '-1.93418062e-2 -6.3875442e-3 -8.75649995e-4'
  [../]

  [./mat_r]
    type = GenericConstantMaterial
    prop_names = 'r11 r12 r44'
    prop_values = '-7.7644e-3 1.63357e-3 7.70775329e-4'
  [../]

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
    type = WallEnergy
    execute_on = 'initial timestep_end'
    polar_x = antiferrodis_A_x
    polar_y = antiferrodis_A_y
    polar_z = antiferrodis_A_z
    G110 = 0.173
    G11_G110 = 0.5
    G12_G110 = 0
    G44_G110 = 0.5
    G44P_G110 = 0.5
  [../]
  [./FcPA]
    type = RotoPolarCoupledEnergyEighth
    execute_on = 'initial timestep_end'
  [../]

  [./FcPu]
    type = ElectrostrictiveCouplingEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./FcAu]
    type = RotostrictiveCouplingEnergy
    execute_on = 'initial timestep_end'
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
  line_search = 'bt'
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
