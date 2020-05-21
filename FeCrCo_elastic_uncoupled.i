[Mesh]
  type = GeneratedMesh  # https://mooseframework.inl.gov/source/mesh/GeneratedMesh.html
  dim = 2

  elem_type = QUAD4

  nx = 15
  ny = 15
  xmax = 15
  ymax = 15
[]

[Variables]             #c2 = Cr, c3 = Co
# Default for variable blocks that do not have anything inbetween them are:
# family = LAGRANGE
# order = FIRST
  [./c2]
  [../]

  [./w2]
  [../]

  [./c3]
  [../]

  [./w3]
  [../]

  [./disp_x]
  [../]

  [./disp_y]
  [../]

  [./mag_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi1
      theta = polar_theta1
      M0s = 1.0 #amplitude of the RandomConstrainedVectorFieldIC
      component  = 0
    [../]
  [../]
  [./mag_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi1
      theta = polar_theta1
      M0s = 1.0 #amplitude of the RandomConstrainedVectorFieldIC
      component  = 1
    [../]
  [../]
  [./mag_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomConstrainedVectorFieldIC
      phi = azimuth_phi1
      theta = polar_theta1
      M0s = 1.0 #amplitude of the RandomConstrainedVectorFieldIC
      component  = 2
    [../]
  [../]

  [./potential_H_int]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -5.0
      max = 5.0
    [../]
  [../]

[]

[AuxVariables]
  #--------------------------------------------#
  #                                            #
  #  note that LAGRANGE/FIRST is not default   #
  #  for AuxVariables                          #
  #                                            #
  #--------------------------------------------#
  [./c1]
    family = LAGRANGE
    order = FIRST
  [../]

  #--------------------------------------------#
  #                                            #
  #  field to seed IC that obeys constraint    #
  #                                            #
  #--------------------------------------------#
  [./azimuth_phi1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 1.2
      seed = 1
    [../]
  [../]
  [./polar_theta1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 1.2
      seed = 2
    [../]
  [../]

[]

[AuxKernels]
  # This Auxkernel is not used for anything just for visualization purposes
  [./c1_define]
    type = ParsedAux # https://mooseframework.org/source/auxkernels/ParsedAux.html
    variable = c1
    function = '1 - c2 - c3'
    args = 'c2 c3'
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  block = 0
[]

[ICs]
  [./c_CrIC]
    type = RandomIC # https://www.mooseframework.org/source/ics/RandomIC.html
    variable = c2
    # Average Composition 40%
    min = 0.400
    max = 0.401
  [../]

  [./c_CoIC]
    type = RandomIC
    variable = c3
    # Average Composition 40%
    min = 0.400
    max = 0.401
  [../]
[]

[Kernels]
# Kernels for Reverse Split Cahn-Hilliard Equations: https://mooseframework.inl.gov/modules/phase_field/Phase_Field_Equations.html

# Implementing Off-diagonal Onsager Matrix with Koyama 2004, Equation 8:
# (Example shown in: https://mooseframework.org/source/kernels/SplitCHWRes.html)
  [./c2_res]
    # Kernel for the chemical potential w2
    type = SplitCHParsed # https://www.mooseframework.org/source/kernels/SplitCHParsed.html
    variable = c2
    f_name = F_sum # Total Free Energy in the system
    kappa_name = kappa_c
    w = w2
    args = 'c3'
  [../]

  [./w22_res]
    # 2nd term of top Eq. 8 in Koyama 2004 paper, del*L_22*del(dG_sys/dc2)
    type = SplitCHWRes # https://mooseframework.org/source/kernels/SplitCHWRes.html
    variable = w2
    mob_name = M_22
  [../]

  [./w23_res]
    # 3rd term of top Eq. 8 in Koyama 2004 paper, del*L_23*del(dG_sys/dc3)
    type = SplitCHWRes
    variable = w2
    w = w3
    mob_name = M_23
  [../]

  [./c3_res]
    # Kernel for the chemical potential w3
    type = SplitCHParsed
    variable = c3
    f_name = F_sum
    kappa_name = kappa_c
    w = w3
    args = 'c2'
  [../]

  [./w33_res]
    # 2nd term of bottom Eq. 8 in Koyama 2004 paper, del*L_32*del(dG_sys/dc2)
    type = SplitCHWRes
    variable = w3
    mob_name = M_33
  [../]

  [./w32_res]
    # 3rd term of bottom Eq. 8 in Koyama 2004 paper, del*L_33*del(dG_sys/dc3)
    type = SplitCHWRes
    variable = w3
    w = w2
    mob_name = M_23
  [../]

  [./time_c2]
    # 1st term of top Eq. 8 in Koyama 2004 paper, dc2/dt
    type = CoupledTimeDerivative # https://mooseframework.org/source/kernels/CoupledTimeDerivative.html
    variable = w2
    v = c2
  [../]

  [./time_c3]
    # 1st term of bottom Eq. 8 in Koyama 2004 paper, dc3/dt
    type = CoupledTimeDerivative
    variable = w3
    v = c3
  [../]

  [./TensorMechanics]
    displacements = 'disp_x disp_y'
  [../]

 #---------------------------------------#
  #                                       #
  #          Time dependence              #
  #                                       #
  #---------------------------------------#

  [./mag1_x_time]
    type = TimeDerivative
    variable = mag_x
  [../]
  [./mag1_y_time]
    type = TimeDerivative
    variable = mag_y
  [../]
  [./mag1_z_time]
    type = TimeDerivative
    variable = mag_z
  [../]


  #---------------------------------------#
  #                                       #
  #    Magnetostatic Poisson equation     #
  #                                       #
  #---------------------------------------#

  [./int_pot_lap]
    type = Electrostatics
    variable = potential_H_int
    block = '0'
    permittivity = 1.0 #a dummy variable at the moment since we use the "electrostatics" kernel
  [../]
  [./int_bc_pot_lap]
    type = MagHStrongCartFeCrCoAlloy
    variable = potential_H_int
    block = '0'
    mag_x = mag_x
    mag_y = mag_y
    mag_z = mag_z
    c1 = c1
    c2 = c2
    c3 = c3
  [../]
[]

[Materials]
  [./C_Fe] # Defining c_1 as a material property so it can be input into the free energy equations
    type = DerivativeParsedMaterial # https://www.mooseframework.org/source/materials/DerivativeParsedMaterial.html
    f_name = c_1
    function = '1 - c2 - c3'
    args = 'c2 c3'
    derivative_order = 2
  [../]

  [./constants] # Constants used in other material properties
    type = GenericConstantMaterial
    prop_names = 'T    bohrM     mu0'
    prop_values = '913 5.788e-5 1.0'
  [../]

  [./atomic_magnetic_moment] # Koyama 2004, Pg.2, 5th equation after Eq. 2: beta^alpha
    type = ParsedMaterial
    f_name = beta
    function = '(2.22*(1-c2-c3) - 0.01*c2 + 1.35*c3 - 0.85*(1-c2-c3)*c2 + (2.4127 + 0.2418*(c3-(1-c2-c3)))*(1-c2-c3)*c3)'
    args = 'c2 c3'
  [../]

  [./Intensity_of_magnetization] # Koyama 2004, 2nd Page (Between Eq 6 & 7)
    type = ParsedMaterial
    f_name = I
    function='if(tau>0.9, 2^(-(2+10*(tau-1))), 1-1/7*(5*tau^4+2*tau^20))'
    args = 'c2 c3'
    material_property_names = 'tau(c2,c3)'
  [../]

  [./magnetic_moment] # Koyama 2004, Equation 7
    type = ParsedMaterial
    f_name = Ms
    function = 'bohrM*beta*I'
    args = 'c2 c3'
    material_property_names = 'bohrM beta(c2,c3) I(c1,c2,c3)'
  [../]

  [./gradient_coef] # Koyama 2004, Table 1: kappa_c = 1.0e-14 (J*m^2/mol)
    type = GenericFunctionMaterial # Used to make a constant property name made up from multiplying several terms together
    prop_names = 'kappa_c'
    # Multiplying kappa_c with eV_J, nm^2, to turn units to (eV*nm^2/mol) and scaling factor d = 1e-27 for better convergence
    prop_values ='1.0e-14*6.24150934e+18*1e+09^2*1e-27'
  [../]

  [./RT_clnc] # Koyama 2004, Second term of G^alpha_c (Equation 2 contribution of the free energy)
    type = DerivativeParsedMaterial
    f_name = RT_clnc
    constant_names = '     T    eV_J            R'
    constant_expressions ='913  6.24150934e+18  8.31446261815324'
    material_property_names = 'c_1(c2,c3)'
    function = 'eV_J*(R*T*(c_1*log(c_1)+c2*log(c2)+c3*log(c3)))'
    args = 'c2 c3'
    derivative_order = 2
  [../]

  [./heat_of_mixing] # Koyama 2004, Third term of G^alpha_c (Equation 2 contribution of the free energy)
    type = DerivativeParsedMaterial
    f_name = E_G
    constant_names = '     eV_J            T'
    constant_expressions ='6.24150934e+18  913'
    material_property_names = 'c_1(c2,c3)'
    function = 'eV_J*((20500 - 9.68*T)*c_1*c2 + (-23669 + 103.9627*T - 12.7886*T*log(T))*c_1*c3 + ((24357 - 19.797*T) - 2010*(c3 - c2))*c2*c3)'
    args = 'c2 c3'
    derivative_order = 2
  [../]

  [./magnetic_contribution_to_Gibbs] # Koyama 2004, Equation 2: ^mgG^alpha
    type = DerivativeParsedMaterial
    f_name = mg_G
    material_property_names = 'f_tau(c2,c3) c_1(c2,c3)'
    constant_names = '     eV_J            T    R'
    constant_expressions ='6.24150934e+18  913  8.31446261815324'
    function = 'eV_J*(R*T*log((2.22*c_1 - 0.01*c2 + 1.35*c3 - 0.85*c_1*c2 + (2.4127 + 0.2418*(c3-c_1))*c_1*c3)+1)*f_tau)' #f_tau Defined below
    args = 'c2 c3'
    derivative_order = 2
  [../]

  [./tau] # Koyama 2004, Equation 2: tau = T/Curie_Temperature
    type = DerivativeParsedMaterial
    f_name = tau
    material_property_names = 'c_1(c2,c3)'
    function = '913/(1043*c_1 - 311.5*c2 + 1450*c3 + (1650 + 550*(c2-c_1))*c_1*c2 + 590*c_1*c3)'
    args = 'c2 c3'
    derivative_order = 2
  [../]

  [./tau_function] # Koyama 2004, Equation 2 ^mgG^alpha term: f(tau)
  # Koyama does not explicitly define f(tau)
  # Expression for f(tau) was obtained by Xiong 2012 (http://dx.doi.org/10.1016/j.calphad.2012.07.002), Equation 5
    type = DerivativeParsedMaterial
    f_name = f_tau
    constant_names = '     p    A'
    constant_expressions ='0.4  1.55828482003'
    function='if(tau<1.00, 1-1/A*(79*tau^-1/(140*p)+474/497*(1/p-1)*(tau^3/6+tau^9/135+tau^15/600)), -1/A*(1/10*tau^-5+1/315*tau^-15+1/1500*tau^-25))'
    args = 'c2 c3'
    material_property_names = 'tau(c2,c3) c_1(c2,c3)'
    derivative_order = 2
  [../]

  # Onsager coefficients: "Mobility_Lij":
  [./Mobility_M22] # Koyama 2004, After Equation 8: L_22
    type = DerivativeParsedMaterial
    f_name = M_22
    constant_names = '     nm_m   eV_J            d      T    R'
    constant_expressions ='1e+09  6.24150934e+18  1e-27  913  8.31446261815324'
    material_property_names = 'c_1(c2,c3)'
    # D*_1, D*_2, and D*_3 where explicitly writting out in this mobility equation
    function = 'nm_m^2/eV_J/d*(c_1*c2*(1.0e-4*exp(-294000/(8.31446261815324*913))) + (1-c2)^2*(2.0e-5*exp(-308000/(8.31446261815324*913))) + c2*c3*(1.0e-4*exp(-294000/(8.31446261815324*913))))*c2/(R*T)'
    args = 'c2 c3'
    derivative_order = 1
  [../]

  [./Mobility_M23] # Koyama 2004, After Equation 8: L_23
    type = DerivativeParsedMaterial # L_23 = L_32
    f_name = M_23
    constant_names = '     nm_m   eV_J            d      T    R'
    constant_expressions ='1e+09  6.24150934e+18  1e-27  913  8.31446261815324'
    material_property_names = 'c_1(c2,c3)'
    # D*_1, D*_2, and D*_3 where explicitly writting out in this mobility equation
    # Many terms in this equation come from Koyama 2004, Table 1
    function = 'nm_m^2/eV_J/d*(c_1*(1.0e-4*exp(-294000/(8.31446261815324*913))) - (1-c2)*(2.0e-5*exp(-308000/(8.31446261815324*913))) - (1-c3)*(1.0e-4*exp(-294000/(8.31446261815324*913))))*c2*c3/(R*T)'
    args = 'c2 c3'
    derivative_order = 1
  [../]

  [./Mobility_M33] # Koyama 2004, After Equation 8: L33
    type = DerivativeParsedMaterial
    f_name = M_33
    #material_property_names = 'D_1  D_2  D_3'
    constant_names = '     nm_m   eV_J            d      T    R'
    constant_expressions ='1e+09  6.24150934e+18  1e-27  913  8.31446261815324'
    material_property_names = 'c_1(c2,c3)'
    # D*_1, D*_2, and D*_3 where explicitly writting out in this mobility equation
    # Many terms in this equation come from Koyama 2004, Table 1
    function = 'nm_m^2/eV_J/d*(c_1*c3*(1.0e-4*exp(-294000/(8.31446261815324*913))) + c2*c3*(2.0e-5*exp(-308000/(8.31446261815324*913))) + (1-c3)^2*(1.0e-4*exp(-294000/(8.31446261815324*913))))*c3/(R*T)'
    args = 'c2 c3'
    derivative_order = 1
  [../]
################################## TensorMechanics Kernels Start Here ##################################

  [./C_Cr] # Putting c2 as a material property for TensorMechanics
    type = DerivativeParsedMaterial
    f_name = c_2
    function = 'c2'
    args = 'c2'
    derivative_order = 2
  [../]
  [./C_Co] # Putting c3 as a material property for TensorMechanics
    type = DerivativeParsedMaterial
    f_name = c_3
    function = 'c3'
    args = 'c3'
    derivative_order = 2
  [../]

  [./Stiffness_Matrix_Fe] # Koyama 2004, Table 1 Elastic Coefficients C^Fe_ij
    type = ComputeElasticityTensor # https://mooseframework.org/source/materials/ComputeElasticityTensor.html
    # The elasticity tensors below were converted from GPa to eV/nm^3.
    C_ijkl = '1454.89 845.35 845.35 1454.89 845.35 1454.89 735.43 735.43 735.43'
    base_name = Fe_mat
    fill_method = symmetric9
  [../]
  [./Stiffness_matrix_Cr] # Koyama 2004, Table 1 Elastic Coefficients C^Cr_ij
    type = ComputeElasticityTensor
    base_name = Cr_mat
    C_ijkl = '2184.53 423.17 423.17 2184.53 423.17 2184.53 629.14 629.14 629.14'
    fill_method = symmetric9
  [../]
  [./Stiffness_Matrix_Co] # Koyama 2004, Table 1 Elastic Coefficients C^Co_ij = C^Fe_ij
    type = ComputeElasticityTensor
    C_ijkl = '1454.89 845.35 845.35 1454.89 845.35 1454.89 735.43 735.43 735.43'
    base_name = Co_mat
    fill_method = symmetric9
  [../]
  [./Composite_Tensor] # This kernel takes all 3 elasticity tensors above and combines them into one
    type = CompositeElasticityTensor # https://www.mooseframework.org/source/materials/CompositeElasticityTensor.html
    args = 'c2 c3'
    tensors = 'Fe_mat Cr_mat Co_mat' # Name of 3 Stiffness Matrix above
    weights = 'c_1    c_2    c_3' # These are supposed to be material properties, reason for c2 and c3 -> MP
  [../]
  [./stress]
    type = ComputeLinearElasticStress # https://mooseframework.org/source/materials/ComputeLinearElasticStress.html
  [../]
  [./strain]
    type = ComputeSmallStrain # https://www.mooseframework.org/source/materials/ComputeSmallStrain.html
    eigenstrain_names = 'eigenstrain'
  [../]
  [./StressFreeStrain_Prefactor] # Koyama 2004, Equation 4
    type = DerivativeParsedMaterial
    f_name = stress_free_strain
    constant_names = '     epsilon2  epsilon3  c2_0  c3_0' #c2_0 & c3_0 are ICs for c2 and c3 = 40% and 40% respectivly
    constant_expressions ='6.1e-3    -7.1e-3   0.40  0.40' #epsilon2 and epsilon3 from Koyama 2004, Table 1 Lattice mismatch
    function = 'epsilon2*(c2 - c2_0) + epsilon3*(c3 - c3_0)'# Equation 4
    args = 'c2 c3'
    derivative_order = 2
  [../]
  [./eigen_strain] # Koyama 2004, Equation 4 d_ij (Kronecker's delta)
    type = ComputeVariableEigenstrain
    eigen_base = '1 1 0 0 0 0'
    prefactor = stress_free_strain
    args = 'c2 c3'
    eigenstrain_name = 'eigenstrain'
  [../]
  [./Elastic_Free_Energy] # Koyama 2004, Top Equation 3
    type = ElasticEnergyMaterial # https://mooseframework.org/source/materials/ElasticEnergyMaterial.html
    f_name = E_str
    args = 'c2 c3'
    derivative_order = 2
  [../]
  ################################## TensorMechanics Kernels End Here ##################################

  [./F_system] # Total free energy of the system. Koyama, 2004 Eq. 1
    type = DerivativeSumMaterial
    f_name = F_sum
    prefactor =     '1e-27  1e-27    1e-27  1e-27' # mutiplying d=1e-27 to each free energy term
    sum_materials = 'mg_G   RT_clnc  E_G    E_str'
    args = 'c2 c3'
    derivative_order = 2
  [../]
[]

[BCs]
  [./Periodic] # Periodic Boundary Conditions explained: https://en.wikipedia.org/wiki/Periodic_boundary_conditions
    [./c_bcs]
      auto_direction = 'x y potential_H_int'
    [../]
  [../]
[]

[Postprocessors]
  [./<M_s>]
    type = ElementIntegralMaterialProperty
    mat_prop = Ms
  [../]
  [./Felast]
    type = ElementIntegralMaterialProperty
    mat_prop = E_str
  [../]
  [./Fmag]
    type = ElementIntegralMaterialProperty
    mat_prop = mg_G
  [../]
  [./Fmix]
    type = ElementIntegralMaterialProperty
    mat_prop = E_G
  [../]
  [./F_RTclnc]
    type = ElementIntegralMaterialProperty
    mat_prop = RT_clnc
  [../]
  [./Ftot]
    type = ElementIntegralMaterialProperty
    mat_prop = F_sum
  [../]
[]

[Preconditioning]
  [./coupled]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason'
    petsc_options_iname = '-pc_type -ksp_gmres_restart '
    petsc_options_value = 'bjacobi       31          '
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'implicit-euler'
  verbose = true

  solve_type = PJFNK

  l_max_its = 30
  l_tol = 1e-6
  nl_max_its = 15
  nl_abs_tol = 1e-9

  end_time = 180000 #50hrs
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    cutback_factor = 0.8
    cutback_factor_at_failure = 0.5
    growth_factor = 1.5
    optimal_iterations = 7
  [../]
  [./Adaptivity]
    coarsen_fraction = 0.1
    refine_fraction = 0.7
    max_h_level = 2
  [../]
[]

[Outputs]
  print_linear_residuals = false
  perf_graph = true
  [./out]
    type = Exodus
    file_base = out_FeCrCo_elastic_uncoupled
    elemental_as_nodal = true
    interval = 1
  [../]
  [./outCSV]
    type = CSV
    file_base = out_FeCrCo_elastic_uncoupled
  [../]
[]
