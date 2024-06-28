[Mesh]
  file = out_BTO_8-120nm_0.5_TWD_T298K_E0.e
[]

[GlobalParams]
  len_scale = 1.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_E_int = potential_E_int

  displacements = 'u_x u_y u_z'

  ##############################################
  ##=
  ##  IMPORTANT(!): Units in Ferret are nm, kg,
  ##                seconds, and attocoulombs
  ##
  ##############################################

  u_x = u_x
  u_y = u_y
  u_z = u_z
[]

[Variables]

  #################################
  ##
  ##  Variable definitions
  ##    P, u, phi, e^global_ij
  ##  and their initial conditions
  ##
  #################################

  [./polar_x]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_x
    block = '1'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_y
    block = '1'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_z
    block = '1'
  [../]

  [./potential_E_int]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
  [../]

  [./potential_H_int]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
  [../]

  [./u_x]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
    initial_from_file_var = u_x
  [../]
  [./u_y]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
    initial_from_file_var = u_y
  [../]
  [./u_z]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
    initial_from_file_var = u_z
  [../]
[]

[Functions]
  [./e_harmonic]
    type = ParsedFunction
    vars = 'A0 omega'
    vals = '0.02 2E10'
    value = 'A0*sin(omega*t)'
  [../]
[]

[AuxVariables]

  ######################################
  ##
  ##  Auxiarilly variable definitions
  ##   (can be intermediate variables
  ##   or for postprocessed quantities)
  ##
  ######################################


  ######################################
  ##
  ##  Stress/strain tensor components
  ##
  ######################################

  [./e00]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./e01]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./e10]
    order = CONSTANT
    family = MONOMIAL
     block = '1 2'
  [../]
  [./e11]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./e12]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./e22]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]

  [./s00]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./s01]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./s10]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./s11]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./s12]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./s22]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]

  [./E_x]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./E_y]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./E_z]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]

  [./H_x]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./H_y]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]
  [./H_z]
    order = CONSTANT
    family = MONOMIAL
    block = '1 2'
  [../]


  [./ind_P_x]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./ind_P_y]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./ind_P_z]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]

  [./P0_x]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./P0_y]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./P0_z]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]

  [./u0_x]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./u0_y]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./u0_z]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]

  [./ind_u_x]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./ind_u_y]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]
  [./ind_u_z]
    order = FIRST
    family = LAGRANGE
    block = '1'
  [../]


  [./m_x]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  [./m_y]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]
  [./m_z]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./dP_dt_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dP_dt_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dP_dt_z]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./pont]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./divP]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

  [./ard]
    order = CONSTANT
    family = MONOMIAL
    block = '1'
  [../]

[]

[AuxKernels]


  [./dp_dt_x]
    type = TimeDerivativeAux
    variable = dP_dt_x
    functor = ind_P_x
    factor = 1
    execute_on = 'TIMESTEP_END'
  [../]
  [./dp_dt_y]
    type = TimeDerivativeAux
    variable = dP_dt_y
    functor = ind_P_y
    factor = 1
    execute_on = 'TIMESTEP_END'
  [../]
  [./dp_dt_z]
    type = TimeDerivativeAux
    variable = dP_dt_z
    functor = ind_P_z
    factor = 1
    execute_on = 'TIMESTEP_END'
  [../]


  ######################################
  ##
  ##  Auxiarilly Kernel definitions
  ##   (can be intermediate "operations"
  ##   or for postprocessed quantities)
  ##
  ######################################

  [./e00]
    type = RankTwoAux
    variable = e00
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 0
    block = '1 2'
  [../]
  [./e01]
    type = RankTwoAux
    variable = e01
    rank_two_tensor = total_strain
    index_i = 0
    index_j = 1
    block = '1 2'
  [../]
  [./e10]
    type = RankTwoAux
    variable = e10
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 0
    block = '1 2'
  [../]
  [./e12]
    type = RankTwoAux
    variable = e12
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 2
    block = '1 2'
  [../]
  [./e11]
    type = RankTwoAux
    variable = e11
    rank_two_tensor = total_strain
    index_i = 1
    index_j = 1
    block = '1 2'
  [../]
  [./e22]
    type = RankTwoAux
    variable = e22
    rank_two_tensor = total_strain
    index_i = 2
    index_j = 2
    block = '1 2'
  [../]

  [./s00]
    type = RankTwoAux
    variable = s00
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    block = '1 2'
  [../]
  [./s01]
    type = RankTwoAux
    variable = s01
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    block = '1 2'
  [../]
  [./s10]
    type = RankTwoAux
    variable = s10
    rank_two_tensor = stress
    index_i = 1
    index_j = 0
    block = '1 2'
  [../]
  [./s12]
    type = RankTwoAux
    variable = s12
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    block = '1 2'
  [../]
  [./s11]
    type = RankTwoAux
    variable = s11
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    block = '1 2'
  [../]
  [./s22]
    type = RankTwoAux
    variable = s22
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    block = '1 2'
  [../]


  [./ex]
    type = ElecFieldAux
    variable = E_x
    component = 0
    block = '1 2'
  [../]
  [./ey]
    type = ElecFieldAux
    variable = E_y
    component = 1
    block = '1 2'
  [../]
  [./ez]
    type = ElecFieldAux
    variable = E_z
    component = 2
    block = '1 2'
  [../]

  [./hx]
    type = QuasistaticFieldAux
    variable = H_x
    component = 0
    potential_int = potential_H_int
    block = '1 2'
  [../]
  [./hy]
    type = QuasistaticFieldAux
    variable = H_y
    component = 1
    potential_int = potential_H_int
    block = '1 2'
  [../]
  [./hz]
    type = QuasistaticFieldAux
    variable = H_z
    component = 2
    potential_int = potential_H_int
    block = '1 2'
  [../]


  [./compute_spont_Px]
    type = ValueAux
    variable = P0_x
    var1 = polar_x
    execute_on = 'INITIAL'
    block = '1'
  [../]
  [./compute_spont_Py]
    type = ValueAux
    variable = P0_y
    var1 = polar_y
    execute_on = 'INITIAL'
    block = '1'
  [../]
  [./compute_spont_Pz]
    type = ValueAux
    variable = P0_z
    var1 = polar_z
    execute_on = 'INITIAL'
    block = '1'
  [../]

  [./compute_spont_ux]
    type = ValueAux
    variable = u0_x
    var1 = u_x
    execute_on = 'INITIAL'
    block = '1'
  [../]
  [./compute_spont_uy]
    type = ValueAux
    variable = u0_y
    var1 = u_y
    execute_on = 'INITIAL'
    block = '1'
  [../]
  [./compute_spont_uz]
    type = ValueAux
    variable = u0_z
    var1 = u_z
    execute_on = 'INITIAL'
    block = '1'
  [../]


  [./compute_ind_Px]
    type = VectorDiffOrSum
    variable = ind_P_x
    var2 = P0_x
    var1 = polar_x
    block = '1'
    diffOrSum = 0
  [../]
  [./compute_ind_Py]
    type = VectorDiffOrSum
    variable = ind_P_y
    var2 = P0_y
    var1 = polar_y
    block = '1'
    diffOrSum = 0
  [../]
  [./compute_ind_Pz]
    type = VectorDiffOrSum
    variable = ind_P_z
    var2 = P0_z
    var1 = polar_z
    block = '1'
    diffOrSum = 0
  [../]

  [./compute_ind_ux]
    type = VectorDiffOrSum
    variable = ind_u_x
    var2 = u0_x
    var1 = u_x
    block = '1'
    diffOrSum = 0
  [../]
  [./compute_ind_uy]
    type = VectorDiffOrSum
    variable = ind_u_y
    var2 = u0_y
    var1 = u_y
    block = '1'
    diffOrSum = 0
  [../]
  [./compute_ind_uz]
    type = VectorDiffOrSum
    variable = ind_u_z
    var2 = u0_z
    var1 = u_z
    block = '1'
    diffOrSum = 0
  [../]


  [./compute_mx]
    type = AFMSpinCurrentLLdot #dummy Kernel for P cross dPdt
    variable = m_x
    Neel_L_x = ind_P_x
    Neel_L_y = ind_P_y
    Neel_L_z = ind_P_z
    dL_dt_x = dP_dt_x
    dL_dt_y = dP_dt_y
    dL_dt_z = dP_dt_z
    component = 0
    block = '1'
    factor = 0.0012378 #output in nuclear magnetons per nm3
  [../]
  [./compute_my]
    type = AFMSpinCurrentLLdot #dummy Kernel for P cross dPdt
    variable = m_y
    Neel_L_x = ind_P_x
    Neel_L_y = ind_P_y
    Neel_L_z = ind_P_z
    dL_dt_x = dP_dt_x
    dL_dt_y = dP_dt_y
    dL_dt_z = dP_dt_z
    component = 1
    block = '1'
    factor = 0.0012378 #output in nuclear magnetons per nm3
  [../]
  [./compute_mz]
    type = AFMSpinCurrentLLdot #dummy Kernel for P cross dPdt
    variable = m_z
    Neel_L_x = ind_P_x
    Neel_L_y = ind_P_y
    Neel_L_z = ind_P_z
    dL_dt_x = dP_dt_x
    dL_dt_y = dP_dt_y
    dL_dt_z = dP_dt_z
    component = 2
    block = '1'
    factor = 0.0012378 #output in nuclear magnetons per nm3
  [../]


  [./pont]
    type = PontryaginDensity
    variable = pont
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
    block = '1'
  [../]

  [./dP]
    type = DivP
    variable = divP
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
    block = '1'
  [../]

  [./ard]
    type = AxionResponseDensity
    variable = ard
    E_x = E_x
    E_y = E_y
    E_z = E_z
    H_x = H_x
    H_y = H_y
    H_z = H_z
    block = '1'
  [../]


[]

[Materials]

  #################################################
  ##
  ## Landau coefficients from Li et al (2001)
  ##
  ##################################################

  [./Landau_P_FE]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '-0.027722 -0.64755 0.323 8.00399999 4.47 4.91 0.0 0.0 0.0 0.0'
    block = '1'
  [../]

  [./Landau_G_FE]
    type = GenericConstantMaterial
    prop_names = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '0.5 0.51 -0.02 0.02 0.0'
    block = '1'
  [../]

  [./mat_C_FE]
    type = GenericConstantMaterial
    prop_names = 'C11 C12 C44'
    prop_values = '275.0 179.0 54.3'
    block = '1'
  [../]
  [./mat_C_sub]
    type = GenericConstantMaterial
    prop_names = 'C11 C12 C44'
    prop_values = '220.0 34.4 161.1'
    block = '2'
  [../]

  ##################################################
  ##=
  ## NOTE: Sign convention in Ferret for the
  ##        electrostrictive coeff. is multiplied by
  ##        an overall factor of (-1)
  ##
  ##################################################

  [./mat_Q]
    type = GenericConstantMaterial
    prop_names = 'Q11 Q12 Q44'
    prop_values = '-0.11 0.045 -0.029'
    block = '1 2'
  [../]

  [./mat_q]
    type = GenericConstantMaterial
    prop_names = 'q11 q12 q44'
    prop_values = '-14.2 0.74 -1.57'
    block = '1 2'
  [../]

  [./Ms]
    type = GenericConstantMaterial
    prop_names = 'Ms'
    prop_values = '6.25e-7'
    block = '1'
  [../]

  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9

   ###############################################
   ##
   ## symmetric9 fill_method is (default)
   ##     C11 C12 C13 C22 C23 C33 C44 C55 C66
   ##
   ###############################################

    C_ijkl = '275.0 179.0 179.0 275.0 179.0 275.0 54.3 54.3 54.3'
    block = '1'
  [../]

  [./elasticity_tensor_2]
    type = ComputeElasticityTensor
    fill_method = symmetric9

   ###############################################
   ##
   ## symmetric9 fill_method is (default)
   ##     C11 C12 C13 C22 C23 C33 C44 C55 C66
   ##
   ###############################################

    C_ijkl = '220.0 34.4 34.4 220.0 34.4 220.0 161.1 161.1 161.1'
    block = '2'
  [../]

  [./stress_1]
    type = ComputeLinearElasticStress
    block = '1'
  [../]
  [./stress_2]
    type = ComputeLinearElasticStress
    block = '2'
  [../]

  [./strain_1]
    type = ComputeSmallStrain
    block = '1'
  [../]
  [./strain_2]
    type = ComputeSmallStrain
    block = '2'
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
    block = '1'
  [../]

  [./permitivitty_2]

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
    prop_values = '2.0'
    block = '2'
  [../]
[]


[Kernels]

  ###############################################
  ##
  ## Physical Kernel operators
  ## to enforce TDLGD evolution
  ##
  ###############################################


  #Elastic problem
  [./TensorMechanics]
    use_displaced_mesh = false
    eigenstrain_names = eigenstrain
  [../]

  [./bed_x]
    type = BulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
    block = '1'
  [../]
  [./bed_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
    block = '1'
  [../]
  [./bed_z]
    type = BulkEnergyDerivativeEighth
    variable = polar_z
    component = 2
    block = '1'
  [../]

  [./walled_x]
    type = WallEnergyDerivative
    variable = polar_x
    component = 0
    block = '1'
  [../]
  [./walled_y]
    type = WallEnergyDerivative
    variable = polar_y
    component = 1
    block = '1'
  [../]
  [./walled_z]
    type = WallEnergyDerivative
    variable = polar_z
    component = 2
    block = '1'
  [../]

  [./walled2_x]
    type = Wall2EnergyDerivative
    variable = polar_x
    component = 0
    block = '1'
  [../]
  [./walled2_y]
    type = Wall2EnergyDerivative
    variable = polar_y
    component = 1
    block = '1'
  [../]
  [./walled2_z]
    type = Wall2EnergyDerivative
    variable = polar_z
    component = 2
    block = '1'
  [../]

  [./electrostr_ux]
    type = ElectrostrictiveCouplingDispDerivative
    variable = u_x
    component = 0
    block = '1'
  [../]
  [./electrostr_uy]
    type = ElectrostrictiveCouplingDispDerivative
    variable = u_y
    component = 1
    block = '1'
  [../]
  [./electrostr_uz]
    type = ElectrostrictiveCouplingDispDerivative
    variable = u_z
    component = 2
    block = '1'
  [../]

  [./electrostr_polar_coupled_x]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_x
    component = 0
    block = '1'
  [../]
  [./electrostr_polar_coupled_y]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_y
    component = 1
    block = '1'
  [../]
  [./electrostr_polar_coupled_z]
    type = ElectrostrictiveCouplingPolarDerivative
    variable = polar_z
    component = 2
    block = '1'
  [../]


  [./polar_x_electric_E]
    type = PolarElectricEStrong
    variable = potential_E_int
    block = '1'
  [../]
  [./FE_E_int]
    type = Electrostatics
    variable = potential_E_int
    block = '1 2'
  [../]

  [./m_H]
    type = MagHStrongCart
    variable = potential_H_int
    block = '1'
    mag_x = m_x
    mag_y = m_y
    mag_z = m_z
  [../]
  [./H_int]
    type = Electrostatics
    variable = potential_H_int
    block = '1 2'
  [../]

  [./polar_electric_px]
    type = PolarElectricPStrong
    variable = polar_x
    component = 0
    block = '1'
  [../]
  [./polar_electric_py]
    type = PolarElectricPStrong
    variable = polar_y
    component = 1
    block = '1'
  [../]
  [./polar_electric_pz]
    type = PolarElectricPStrong
    variable = polar_z
    component = 2
    block = '1'
  [../]

  [./polar_x_time]
    type = TimeDerivativeScaled
    variable=polar_x
    time_scale = 5.0e-12
    block = '1'
  [../]
  [./polar_y_time]
    type = TimeDerivativeScaled
    variable = polar_y
    time_scale = 5.0e-12
    block = '1'
  [../]
  [./polar_z_time]
    type = TimeDerivativeScaled
    variable = polar_z
    time_scale = 5.0e-12
    block = '1'
  [../]

[]


[BCs]
  [./boundary_infinity_grounding]
    type = DirichletBC
    boundary = '1'
    variable = potential_E_int
    value = 0.0
  [../]
  [./boundary_field]
    type = FunctionDirichletBC
    boundary = '2'
    variable = potential_E_int
    function = e_harmonic
  [../]
  [./boundary_mech_grounding_u_x]
    type = DirichletBC
    boundary = '1 2 3 4 5 6'
    variable = u_x
    value = 0.0
  [../]
  [./boundary_mech_grounding_u_y]
    type = DirichletBC
    boundary = '1 2 3 4 5 6'
    variable = u_y
    value = 0.0
  [../]
  [./boundary_mech_grounding_u_z]
    type = DirichletBC
    boundary = '1 2 3 4 5 6'
    variable = u_z
    value = 0.0
  [../]


  [./boundary_M]
    type = DirichletBC
    boundary = '1 2 3 4 5 6'
    variable = potential_H_int
    value = 0.0
  [../]

 []

[Postprocessors]

  ###############################################
  ##
  ##  Postprocessors (integrations over the
  ##  computational domain) to calculate the total
  ##  energy decomposed into linear combinations of
  ##  the different physics.
  ##
  ###############################################


  # Volume-integrated quantities

  [./ave_indP_x]
    type = ElementAverageValue
    variable = ind_P_x
    execute_on = 'timestep_end'
    block = '1'
  [../]
  [./ave_indP_y]
    type = ElementAverageValue
    variable = ind_P_y
    execute_on = 'timestep_end'
    block = '1'
  [../]
  [./ave_indP_z]
    type = ElementAverageValue
    variable = ind_P_z
    execute_on = 'timestep_end'
    block = '1'
  [../]

  [./ave_m_x]
    type = ElementAverageValue
    variable = m_x
    execute_on = 'timestep_end'
    block = '1'
  [../]
  [./ave_m_y]
    type = ElementAverageValue
    variable = m_y
    execute_on = 'timestep_end'
    block = '1'
  [../]
  [./ave_m_z]
    type = ElementAverageValue
    variable = m_z
    execute_on = 'timestep_end'
    block = '1'
  [../]

  [./ave_ard]
    type = ElementAverageValue
    variable = ard
    execute_on = 'timestep_end'
    block = '1'
  [../]

  [./pont_ave]
    type = ElementAverageValue
    variable = pont
    execute_on = 'timestep_end'
    block = '1'
  [../]

  [./surfVc_indP_x]
    type = PointValue
    point = '0.0 0.0 3.95'
    variable = ind_P_x
    execute_on = 'timestep_end'
  [../]
  [./surfVc_indP_y]
    type = PointValue
    point = '0.0 0.0 3.95'
    variable = ind_P_y
    execute_on = 'timestep_end'
  [../]
  [./surfVc_indP_z]
    type = PointValue
    point = '0.0 0.0 3.95'
    variable = ind_P_z
    execute_on = 'timestep_end'
  [../]


  [./surf_m_x]
    type = SideAverageValue
    variable = m_x
    execute_on = 'timestep_end'
    boundary = '7'
  [../]
  [./surf_m_y]
    type = SideAverageValue
    variable = m_y
    execute_on = 'timestep_end'
    boundary = '7'
  [../]
  [./surf_m_z]
    type = SideAverageValue
    variable = m_z
    execute_on = 'timestep_end'
    boundary = '7'
  [../]

	# Energies

  [./Fbulk]
    type = BulkEnergyEighth
    execute_on = 'timestep_end'
    block = '1'
  [../]
  [./Fwall]
    type = WallEnergy
    execute_on = 'timestep_end'
    block = '1'
  [../]
  [./Felastic]
    type = ElasticEnergy
    execute_on = 'timestep_end'
    use_displaced_mesh = false
    block = '1 2'
  [../]
  [./Fcoupled]
    type = ElectrostrictiveCouplingEnergy
    execute_on = 'timestep_end'
    block = '1'
  [../]
  [./Felec]
    type = ElectrostaticEnergy
    execute_on = 'timestep_end'
    block = '1'
  [../]
  [./Ftotal]
    type = LinearCombinationPostprocessor
    pp_names = 'Fbulk Fwall Fcoupled Felec'
    pp_coefs = '0.160218 0.160218 0.160218 0.160218' #converted to eV
    execute_on = 'timestep_end'
  [../]

  # Other

  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = Ftotal
    execute_on = 'timestep_end'
  [../]
  [./elapsed]
    type = PerfGraphData
    section_name = "Root"  # for profiling the problem [on]
    data_type = total
  [../]
[]

[Preconditioning]

  ###############################################
  ##
  ##  Numerical preconditioning/solver options
  ##
  ###############################################

  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    80             1e-8        1e-6      1e-5       bjacobi'
  [../]
[]

[Executioner]

  ##########################################
  ##
  ##  Time integration/solver options
  ##
  ##########################################

  type = Transient
  solve_type = 'NEWTON'
  scheme = 'bdf2'
  dtmin = 1e-18
  dtmax = 1.0e-5

  l_max_its = 200

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    growth_factor = 1.2
    cutback_factor = 0.75
    linear_iteration_ratio = 1000
    dt = 5E-12
  [../]
  verbose = true
  nl_max_its = 20
  num_steps = 2000
[]

[Outputs]

  ###############################################
  ##==
  ##  Output options
  ##
  ###############################################

  print_linear_residuals = false
  perf_graph = false

  [./out_exodus]
    type = Exodus
   file_base = out_aBTO_8-200nm_0.5_T298K_VzEx_E2e-2_w2E10
    elemental_as_nodal = true
    interval = 250
  [../]

  [./out_csv]
    type = CSV
   file_base = out_aBTO_8-200nm_0.5_T298K_VzEx_E2e-2_w2E10
    interval = 1
  [../]
[]
