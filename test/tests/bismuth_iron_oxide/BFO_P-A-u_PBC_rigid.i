[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 3
  ny = 3
  nz = 3
  xmin = -1.0
  xmax = 1.0
  ymin = -1.0
  ymax = 1.0
  zmin = -1.0
  zmax = 1.0
  elem_type = HEX8
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



  #factor of two here? q = 2 C Q?. Last set gives shear ~ 1.0 whereis new set gives normal and shear ~1. Note a factor of 8 difference!
  #Q11 = 2.03079   #4.06157    # 12.1847    #0.676929
  #Q12 = -0.357442 #-0.714883  # -2.14465  #-0.119147
  #Q44 = 2.08738    # 4.17475    # 16.699      #4.17474
  #R11 = 1.21666   #2.43332    # 7.29995     #0.405553
  #R12 = -0.51628  #-1.03256   # -3.09768    #-0.172093 
  #R44 = -1.83738   # -3.67475   # -14.699     #-3.67475


  #from mathematica

  C11 = 1.762732e-2
  C12 = 9.4905087e-3
  C44 = 5.24373333e-5

  Q11 = -1.35386
  Q12 = 0.238295
  Q44 = -4.17474 #16.698978569333416
  R11 = -2.35149
  R12 = 1.1143805
  R44 = 3.67475 #-14.69898038101976

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

  disp_x = disp_x
  disp_y = disp_y
  disp_z = disp_z
  
  displacements = 'disp_x disp_y disp_z'

  #potential_int = potential_int
[]


[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.85
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.85
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.85
    [../]
  [../]
  [./antiferrodis_A_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.85
    [../]
  [../]
  [./antiferrodis_A_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.85
    [../]
  [../]
  [./antiferrodis_A_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ConstantIC
      value = 0.85
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
  [./stress_xx_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xx_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xy_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yz_elastic]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
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

[Kernels]

  #Elastic problem
  [./TensorMechanics]
  #This is an action block
  [../]

  [./rotostr_ux]
    type = RotostrictiveCouplingDispDerivative
    variable = disp_x
    component = 0
  [../]
  [./rotostr_uy]
    type = RotostrictiveCouplingDispDerivative
    variable = disp_y
    component = 1
  [../]
  [./rotostr_uz]
    type = RotostrictiveCouplingDispDerivative
    variable = disp_z
    component = 2
  [../]

  [./electrostr_ux]
    type = ElectrostrictiveCouplingDispDerivative
    variable = disp_x
    component = 0
  [../]
  [./electrostr_uy]
    type = ElectrostrictiveCouplingDispDerivative
    variable = disp_y
    component = 1
  [../]
  [./electrostr_uz]
    type = ElectrostrictiveCouplingDispDerivative
    variable = disp_z
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

[Materials]
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '1.762732e-2 9.4905087e-3 9.4905087e-3 1.762732e-2 9.4905087e-3 1.762732e-2 5.24373333e-5 5.24373333e-5 5.24373333e-5'
  [../]
  [./strain_1]
    type = ComputeSmallStrain
  [../]
  [./stress_1]
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
    type = SumEightPostprocessors
    var0 = FbP
    var1 = FbA
    var2 = FgP
    var3 = FgA
    var4 = FcPA
    var5 = FcPu
    var6 = FcAu
    var7 = Felu
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

[Problem]
  null_space_dimension = 6
[]

[UserObjects]
 [./rigidbodymodes_x]
    type = RigidBodyModes3D
    subspace_name = NullSpace
    subspace_indices = '0 1 2 3 4 5'
    modes = 'trans_x trans_y trans_z rot_x rot_y rot_z'
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
  dt = 0.25
  solve_type = 'NEWTON'
  scheme = 'bdf2'
  dtmin = 1e-13
  dtmax = 0.25
[]

[Outputs]
  print_linear_residuals = false
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_BFO_P-A-u
    elemental_as_nodal = true
  [../]
[]
