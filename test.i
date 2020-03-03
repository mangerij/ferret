[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 40
    ny = 40
    xmin = -40.0
    xmax = 40.0
    ymin = -40.0
    ymax = 40.0
  []
  [./cnode]
    input = gen
    type = ExtraNodesetGenerator
    coord = '0.0 0.0'
    new_boundary = 100
  [../]
[]


[GlobalParams]

  permittivity = 1.0           # dummy scalar variable for Electrostatics [Magnetostatics] Kernel

  c1 = c1
  c2 = c2
  c3 = c3

  T = 913
  bohrM = 5.7883818012e-5
  displacements = 'u_x u_y'
[]


[ICs] 
  #------------------------------------------------------#
  #                                                      #
  # Initial Conditions for the concentration variables   #
  # Koyama 2004, Fig 3 & 4 = Fe-40at.%Cr-40at.%Co        #
  #                                                      #
  #------------------------------------------------------#

  [./c_CrIC]
    type = RandomIC
    variable = c2
    min = 0.39
    max = 0.41
  [../]

  [./c_CoIC]
    type = RandomIC
    variable = c3
    min = 0.39
    max = 0.41
  [../]
[]



[Variables] 
  #---------------------------------------------------#
  # c1 = Fe (aux), c2 = Cr, c3 = Co                   #
  #---------------------------------------------------#

  #---------------------------------------------------#  
  #                                                   #
  # Mole fraction of Cr (unitless)                    #
  #                                                   #
  #---------------------------------------------------#  

  [./c2]
    order = FIRST
    family = LAGRANGE
  [../]

  [./w2]  # Chemical potential of Cr (eV/mol)
    order = FIRST
    family = LAGRANGE
  [../]

  [./c3]  # Mole fraction of Co (unitless)
    order = FIRST
    family = LAGRANGE
  [../]

  [./w3]  # Chemical potential of Co (eV/mol)
    order = FIRST
    family = LAGRANGE
  [../]

  #---------------------------------------------------#  
  #                                                   #
  # 2D elastic model                                  #
  #                                                   #
  #---------------------------------------------------#  

  [./u_x]
    order = FIRST
    family = LAGRANGE
  [../]

  [./u_y]
    order = FIRST
    family = LAGRANGE
  [../]

  [./global_strain]
    order = THIRD
    family = SCALAR
  [../]

  #---------------------------------------------------#  
  #                                                   #
  # Magnetostatic potential                           #
  #                                                   #
  #---------------------------------------------------#  

  [./potential_H_int]
    order = FIRST
    family = LAGRANGE
  [../]
[]


[AuxVariables]
  [./c1]
    order = FIRST
    family = LAGRANGE
  [../]

  [./disp_x]
  [../]
  [./disp_y]
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

  #---------------------------------------------------#  
  #                                                   #
  # Magnetic field calculated from concentrations     #
  #                                                   #
  #---------------------------------------------------#  

  [./magnetic_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./magnetic_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./magnetic_z]
    order = FIRST
    family = LAGRANGE
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



  [./calculate_c1]
    type = ParsedAux
    variable = c1
    function = '1.00 - (c2 + c3)' # Function for computing c1 (composition of Fe)
    args = 'c2 c3'
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  [../]

  #-------------------------------------------------------#  
  #                                                       #
  # The following AuxKernels calculate the magnetization  #
  # as a function of the concentration using Eq. (7)      #
  # from Koyama. This is then fed to the Poisson solver   #
  # at every time step to calculate the scalar field      #
  # arising from the magnetic dipoles.                    #
  #                                                       #
  #-------------------------------------------------------#

  [./calculate_mx_c]
    type = CalcMagFeCrCo
    variable = magnetic_x
    component = 0
  [../]
  [./calculate_my_c]
    type = CalcMagFeCrCo
    variable = magnetic_y
    component = 1
  [../]
  [./calculate_mz_c]
    type = CalcMagFeCrCo
    variable = magnetic_z
    component = 2
  [../]
[]

[Kernels]
# Implementing Off-diagonal Onsager Matrix with Koyama 2004, Equation 8:
# (Example shown in: https://mooseframework.org/source/kernels/SplitCHWRes.html)
  [./c2_res]
    type = ADSplitCHParsed
    variable = c2
    f_name = F_total # G_sys
    kappa_name = kappa_c
    w = w2
    args = 'c3 c1'
  [../]

  [./w22_res]
    type = ADSplitCHWRes
    variable = w2
    mob_name = L_22
  [../]

  [./w23_res]
    type = ADSplitCHWRes
    variable = w2
    w = w3
    mob_name = L_23
  [../]

  [./c3_res]
    type = ADSplitCHParsed
    variable = c3
    f_name = F_total
    kappa_name = kappa_c
    w = w3
    args = 'c2 c1'
  [../]

  [./w33_res]
    type = ADSplitCHWRes
    variable = w3
    mob_name = L_33
  [../]

  [./w32_res]
    type = ADSplitCHWRes
    variable = w3
    w = w2
    mob_name = L_32
  [../]

  [./time_c2]
    type = ADCoupledTimeDerivative
    variable = w2
    v = c2
  [../]

  [./time_c3]
    type = ADCoupledTimeDerivative
    variable = w3
    v = c3
  [../]

  [./TensorMechanics]
    displacements = 'u_x u_y'
  [../]

  [./int_pot_lap]
    type = Electrostatics                  # Note that this is actually a kernel for the magnetostatic Poisson equation.
    variable = potential_H_int
  [../]

  [./int_bc_pos]
    type = MagHStrongCartFeCrCoAlloy       # this is the RHS of the Poisson equation for the specific magnetization 
    variable = potential_H_int             # from the Koyama paper. Some comments: the automatic jacobian debugger will
    mag_x = magnetic_x                     # indicate that this is incorrect but that is not true.
    mag_y = magnetic_y                     # any other functional dependence of c1,c2,c3 on magnetization is hard-coded
    mag_z = magnetic_z
  [../]
[]

[Materials]
  #-------------------------------------------------------------------#
  #                                                                   #
  # We need to change the length scale to units of nanometers.        #
  # and change the energy scale to units of electron volts.           #
  # The conversion from meters to nanometers is 1e+09                 #
  # The conversion from joules to electron volts is 6.24150934e+18    #
  #                                                                   #
  # Koyama 2004, Table 1: kappa_c = 1.0e-14 (J*m^2/mol)               #
  #                                                                   #
  #-------------------------------------------------------------------#

  [./gradient_coef_kappa_c]
    type = GenericFunctionMaterial
    prop_names = 'kappa_c'
    prop_values = '1.0e-14*6.24150934e+18*1e+09^2*1e-27' # eV*nm^2/mol
  [../]

  [./constants] # Constants used in other material properties
    type = GenericConstantMaterial
    prop_names = ' p    T    nm_m   eV_J            d      R                 Q_1     Q_2     Q_3     D01     D02     D03     bohr'
    prop_values = '0.4  913  1e+09  6.24150934e+18  1e-27  8.31446261815324  294000  308000  294000  1.0e-4  2.0e-5  1.0e-4  5.7883818012e-5'
  [../] # bohr is Bohr magneton = eV/T

  [./RT_clnc] # Koyama 2004, Second term of G^alpha_c (Equation 2 contribution of the free energy)
    type = DerivativeParsedMaterial
    f_name = RT_clnc
    material_property_names = 'T eV_J d R'
    args = 'c2 c3 c1'
    function = '1*eV_J*d*(R*T*(c1*log(c1)+c2*log(c2)+c3*log(c3)))'
    derivative_order = 2
  [../]

  [./heat_of_mixing] # Koyama 2004, Third term of G^alpha_c (Equation 2 contribution of the free energy)
    type = DerivativeParsedMaterial
    f_name = E_G
    material_property_names = 'T eV_J d Lalpha_12(T) Lalpha_13(T) Lalpha_23(c2,c3,T)'
    function = 'eV_J*d*(Lalpha_12*c1*c2 + Lalpha_13*c1*c3 + Lalpha_23*c2*c3)'
    args = 'c2 c3 c1'
    derivative_order = 2
  [../]

  [./magnetic_contribution_to_Gibbs] # Koyama 2004, Equation 2: ^mgG^alpha
    type = DerivativeParsedMaterial
    f_name = mg_G
    material_property_names = 'g(c1,c2,c3) beta(c1,c2,c3) T eV_J d R'
    function = '1*eV_J*d*(R*T*log(beta+1)*g)' #g is a tau function. Defined below
    args = 'c2 c3 c1'
    derivative_order = 2
  [../]

  [./tau] # Koyama 2004, Equation 2: tau = T/(T^alpha_C)
    type = ParsedMaterial
    f_name = tau
    function = 'T/Tc'
    material_property_names = 'Tc(c1,c2,c3) T'
    args = 'c1 c2 c3'
    outputs = exodus
  [../]

  [./tau_function] # Koyama 2004, Equation 2: f(tau)
  # Koyama does not explicitly define f(tau)
  # Expression for f(tau) was obtained by Xiong 2012, Equation 5
    type = ParsedMaterial
    f_name = g
    function='if(tau<1, 1-1/A*(79*tau^-1/(140*p)+474/497*(1/p-1)*(tau^3/6+tau^9/135+tau^15/600)), -1/A*(1/10*tau^-5+1/315*tau^-15+1/1500*tau^-25))'
    args = 'c1 c2 c3'
    material_property_names = 'p A tau(c1,c2,c3)'
  [../]

  [./A] #Xiong 2012, Equation 6
  # A is a variable used to calculate tau function
    type = ParsedMaterial
    f_name = A
    function = '518/1125 + 11692/15975*(1/p-1)'
    material_property_names = 'p'
  [../]

  [./Interaction_parameter_1_2] # Koyama 2004, Pg.2, 1st equation after Eq. 2: L^alpha_1,2
    type = DerivativeParsedMaterial
    f_name = Lalpha_12
    material_property_names = 'T'
    function = '20500 - 9.68*T'
    derivative_order = 2
  [../]

  [./Interaction_parameter_1_3] # Koyama 2004, Pg.2, 2nd equation after Eq. 2: L^alpha_1,3
    type = DerivativeParsedMaterial
    f_name = Lalpha_13
    material_property_names = 'T'
    function = '-23669 + 103.9627*T - 12.7886*T*log(T)'
    derivative_order = 2
  [../]

  [./Interaction_parameter_2_3] # Koyama 2004, Pg.2, 3rd equation after Eq. 2: L^alpha_2,3
    type = DerivativeParsedMaterial
    f_name = Lalpha_23
    material_property_names = 'T'
    function = '(24357 - 19.797*T) - 2010*(c3 - c2)'
    args = 'c2 c3'
    derivative_order = 2
  [../]

  [./Curie_Temperature] # Koyama 2004, Pg.2, 4th equation after Eq. 2: T^alpha_C
    type = ParsedMaterial
    f_name = Tc
    function = '1043*c1 - 311.5*c2 + 1450*c3 + (1650 + 550*(c2-c1))*c1*c2 + 590*c1*c3'
    args = 'c1 c2 c3'
    outputs = exodus
  [../]
  [self_diffusion_coef_1] # Koyama 2004, Table 1
    type = DerivativeParsedMaterial
    f_name = D_1
    material_property_names = 'D01 Q_1 R T'
    function = 'D01*exp(-Q_1/(R*T))'
    derivative_order = 1
  []

  [self_diffusion_coef_2] # Koyama 2004, Table 1
    type = DerivativeParsedMaterial
    f_name = D_2
    material_property_names = 'D02 Q_2 R T'
    function = 'D02*exp(-Q_2/(R*T))'
    derivative_order = 1
  []

  [self_diffusion_coef_3] # Koyama 2004, Table 1
    type = DerivativeParsedMaterial
    f_name = D_3
    material_property_names = 'D03 Q_3 R T'
    function = 'D03*exp(-Q_3/(R*T))'
    derivative_order = 1
  []

# Onsager coefficients: "Mobility_Lij":
  [Mobility_L22] # Koyama 2004, After Equation 8: L22
    type = DerivativeParsedMaterial
    f_name = L_22
    material_property_names = 'nm_m  eV_J  d  R  T  D_1(D01,Q_1,R,T) D_2(D02,Q_2,R,T) D_3(D03,Q_3,R,T)'
    function = 'nm_m^2/eV_J/d*(c1*c2*D_1 + (1-c2)^2*D_2 + c2*c3*D_3)*c2/(R*T)'
    args = 'c1 c2 c3'
    derivative_order = 1
  []

  [Mobility_L23] # Koyama 2004, After Equation 8: L23
    type = DerivativeParsedMaterial
    f_name = L_23
    material_property_names = 'nm_m  eV_J  d  R  T  D_1(D01,Q_1,R,T) D_2(D02,Q_2,R,T) D_3(D03,Q_3,R,T)'
    function = 'nm_m^2/eV_J/d*(c1*D_1 - (1-c2)*D_2 - (1-c3)*D_3)*c2*c3/(R*T)'
    args = 'c1 c2 c3'
    derivative_order = 1
  []

  [./Mobility_L32] # Koyama 2004, After Equation 8: L32 = L23
    type = DerivativeParsedMaterial
    f_name = L_32
    material_property_names = 'nm_m  eV_J  d  R  T  D_1(D01,Q_1,R,T) D_2(D02,Q_2,R,T) D_3(D03,Q_3,R,T)'
    function = 'nm_m^2/eV_J/d*(c1*D_1 - (1-c2)*D_2 - (1-c3)*D_3)*c2*c3/(R*T)'
    args = 'c1 c2 c3'
    derivative_order = 1
  [../]

  [./Mobility_L33] # Koyama 2004, After Equation 8: L33
    type = DerivativeParsedMaterial
    f_name = L_33
    material_property_names = 'nm_m  eV_J  d  R  T  D_1(D01,Q_1,R,T) D_2(D02,Q_2,R,T) D_3(D03,Q_3,R,T)'
    function = 'nm_m^2/eV_J/d*(c1*c3*D_1 + c2*c3*D_2 + (1-c3)^2*D_3)*c3/(R*T)'
    args = 'c1 c2 c3'
    derivative_order = 1
  [../]

##################### START OF ELASTIC STRAIN ENERGY CONTRIBUTION #####################
# Elastic coefficients converted from GPa to eV/nm^2

  [./elasticity_tensor_Fe] # Koyama 2004, Table 1
    type = ComputeElasticityTensor
    C_ijkl = '1454.89 845.35 845.35 1454.89 845.35 1454.89 735.43 735.43 735.43'
    fill_method = symmetric9
    base_name = iron_el
  [../]

  [./elasticity_tensor_Cr]  # Koyama 2004, Table 1
    type = ComputeElasticityTensor
    C_ijkl = '2184.53 423.17 423.17 2184.53 423.17 2184.53 629.14 629.14 629.14'
    fill_method = symmetric9
    #base_name = chromium_el
  [../]

  [./elasticity_tensor_Co]  # Koyama 2004, Table 1
    type = ComputeElasticityTensor
    base_name = cobalt_el
    C_ijkl = '1454.89 845.35 845.35 1454.89 845.35 1454.89 735.43 735.43 735.43'
    fill_method = symmetric9
  [../]

  [./prefactorCr] # Koyama 2004, Equation 4 = [epsilon2{c2-c20}...]
    type = DerivativeParsedMaterial
    f_name = prefactorCr
    constant_names       = 'epsilon2   c20'
    constant_expressions = '6.1e-3     .4'
    args = 'c2'
    function = 'epsilon2*(c2 - c20)'
  [../]

  [./prefactorCo] # Koyama 2004, Equation 4 = [... + epsilon3{c3-c30}]*del_ij
    type = DerivativeParsedMaterial
    f_name = prefactorCo
    constant_names       = 'epsilon3    c30'
    constant_expressions = '-7.1e-3     .4'
    args = 'c3'
    function = 'epsilon3*(c3 - c30)'
  [../]

  [./StressFreeStrainFe]
    type = ComputeEigenstrain
    eigen_base = 0
    eigenstrain_name = 'eigenstrainFe'
    base_name = iron_el
  [../]

  [./StressFreeStrainCr]
    type = ComputeVariableEigenstrain
    prefactor = prefactorCr
    eigen_base = 1
    args = 'c2'
    eigenstrain_name = 'eigenstrainCr'
    #base_name = chromium_el
  [../]

  [./StressFreeStrainCo]
    type = ComputeVariableEigenstrain
    prefactor = prefactorCo
    eigen_base = 1
    args = 'c3'
    eigenstrain_name = 'eigenstrainCo'
    base_name = cobalt_el
  [../]

  [./stressCr]
    type = ComputeLinearElasticStress
#    base_name = chromium_el
  [../]

  [./stressFe]
    type = ComputeLinearElasticStress
    base_name = iron_el
  [../]

  [./stressCo]
    type = ComputeLinearElasticStress
    base_name = cobalt_el
  [../]

  [./strainFe]
    type = ComputeSmallStrain
    displacements = 'u_x u_y'
    eigenstrain_names = 'eigenstrainFe'
    base_name = iron_el
  [../]

  [./strainCr]
    type = ComputeSmallStrain
    displacements = 'u_x u_y'
    eigenstrain_names = 'eigenstrainCr'
    #base_name = chromium_el
  [../]

  [./strainCo]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    eigenstrain_names = 'eigenstrainCo'
    base_name = cobalt_el
  [../]

  [./global_strain]
    type = ComputeGlobalStrain
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
  [../]

  [./ElasticEnergyFe]
    type = ElasticEnergyMaterial
    f_name = ElasticFe
    base_name = iron_el
    args = 'c1'
    derivative_order = 2
  [../]

  [./ElasticEnergyCr]
    type = ElasticEnergyMaterial
    f_name = ElasticCr
    args = 'c2'
    derivative_order = 2
  [../]

  [./ElasticEnergyCo]
    type = ElasticEnergyMaterial
    f_name = ElasticCo
    base_name = cobalt_el
    args = 'c3'
    derivative_order = 2
  [../]

  ##################### END OF ELASTIC STRAIN ENERGY CONTRIBUTION #####################

  [./G_system] # Total free energy of the system
    type = DerivativeSumMaterial
    f_name = F_total
    sum_materials = 'E_G RT_clnc mg_G ElasticCr ElasticFe ElasticCo'
    args = 'c1 c2 c3'
    outputs = exodus
    derivative_order = 2
  [../]
[]


[BCs] # Boundary Conditions
  [Periodic]
    [c_bcs]
      variable = 'c2 c3 u_x u_y'
      auto_direction = 'x y'
    []
  []

  #------------------------------------------------------#
  #                                                      #
  # BC for magnetostatic potential                       #
  # Essentially "grounded" values at the boundaries.     #
  #                                                      #
  #------------------------------------------------------#

  [./bc_int_pot_R]
    type = DirichletBC
    variable = potential_H_int
    value = 0.0
    boundary = 'top bottom left right' #+y,-y,-x,+x

  [../]

  #------------------------------------------------------#
  #                                                      #
  # what about global strains???                         #
  #                                                      #
  #                                                      #
  #------------------------------------------------------#

  # fix center point location
  [./centerfix_x]
    type = DirichletBC
    boundary = 100
    variable = u_x
    value = 0.0
  [../]
  [./centerfix_y]
    type = DirichletBC
    boundary = 100
    variable = u_y
    value = 0.0
  [../]
[]

[ScalarKernels]
  [./global_strain]
    type = GlobalStrain
    variable = global_strain
    global_strain_uo = global_strain_uo
  [../]
[]

[UserObjects]
  [./global_strain_uo]
    type = GlobalStrainUserObject
    execute_on = 'Initial Linear Nonlinear'
  [../]
[]

[Postprocessors]
  [./Ftot]
    type = ElementIntegralVariablePostprocessor
    variable = F_total
    execute_on = 'initial timestep_end'
  [../]
[]


[Preconditioning]
  [./coupled]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew -snes_converged_reason'
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    31               1e-10      1e-8      1e-6       bjacobi'
  [../]
[]


[Debug]
  show_var_residual_norms = false
[]

[Executioner]
  type = Transient
  steady_state_detection = true

  scheme = bdf2
  solve_type = 'PJFNK'

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -ksp_gmres_restart  -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm         31     lu      3'


  #automatic_scaling = true
  l_max_its = 30
  l_tol = 1e-6
  nl_max_its = 50
  nl_abs_tol = 1e-8
  end_time = 300   # 1h

  verbose = true

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    cutback_factor = 0.8
    growth_factor = 1.2
    optimal_iterations = 9
  []

  #[Adaptivity]
  #  coarsen_fraction = 0.1
  #  refine_fraction = 0.7
  #  max_h_level = 2
  #[]
[]

[Outputs]
  print_linear_residuals = false
  exodus = true
[]
