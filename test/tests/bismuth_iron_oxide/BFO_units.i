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

  [./eigs00]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eigs11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eigs22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eigs01]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eigs12]
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

  [./eigs00]
    type = LocalBFOEigenstressAux
    variable = eigs00
    index_i = 0
    index_j = 0
  [../]
  [./eigs11]
    type = LocalBFOEigenstressAux
    variable = eigs11
    index_i = 1
    index_j = 1
  [../]
  [./eigs22]
    type = LocalBFOEigenstressAux
    variable = eigs22
    index_i = 2
    index_j = 2
  [../]
  [./eigs01]
    type = LocalBFOEigenstressAux
    variable = eigs01
    index_i = 0
    index_j = 1
  [../]
  [./eigs12]
    type = LocalBFOEigenstressAux
    variable = eigs12
    index_i = 1
    index_j = 2
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


  [./rotostr_dis_coupled_x]
    type = RotostrictiveCouplingDistortDerivative
    variable = antiferrodis_A_x
    component = 0
  u_x = disp_x
  u_y = disp_y
  u_z = disp_z
  [../]
  [./rotostr_dis_coupled_y]
    type = RotostrictiveCouplingDistortDerivative
    variable = antiferrodis_A_y
    component = 1
  u_x = disp_x
  u_y = disp_y
  u_z = disp_z
  [../]
  [./rotostr_dis_coupled_z]
    type = RotostrictiveCouplingDistortDerivative
    variable = antiferrodis_A_z
    component = 2
  u_x = disp_x
  u_y = disp_y
  u_z = disp_z
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


#hmmm
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

#the point is the following: use a slightly different definition of Q_ij than Hlinka

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


  [./FcPu]
    type = ElectrostrictiveCouplingEnergy
    execute_on = 'initial timestep_end'
  u_x = u_x
  u_y = u_y
  u_z = u_z
  [../]
  [./FcAu]
    type = RotostrictiveCouplingEnergy
    execute_on = 'initial timestep_end'
  u_x = u_x
  u_y = u_y
  u_z = u_z
  [../]

  [./Felu]
    type = ElasticEnergy
    execute_on = 'initial timestep_end'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'FbP FbA FcPu FcAu Felu'
    pp_coefs = ' 1 1 1 1 1'
    execute_on = 'initial timestep_end'
  [../]
 # [./perc_change]
 #   type = EnergyRatePostprocessor
 #   postprocessor = Ftot
 #   execute_on = 'initial timestep_end'
 #   dt = dt
 # [../]
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
 # [./kill]
 #  type = Terminator
 #  expression = 'perc_change <= 1.0e-8'
 # [../]
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
  type = Steady
  #dt = 0.03
  solve_type = 'NEWTON'
 # scheme = 'bdf2'
 # dtmin = 1e-13
 # dtmax = 10.0
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_bfo
    elemental_as_nodal = true
  [../]
[]
