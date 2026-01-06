

  ##############################
  ##
  ## UNITS:
  ##
  ##   gamma = (2.2101*10^5 / 2 pi) m/C
  ##
  ##   gamma/(mu*0Ms) = 48291.9 nm*mus/pg
  ##
  ##   coefficients given in (pg/nm*mus)
  ##   which is equivalent to an energy/vol
  ##
  ##   Effective fields are 1/(mu0*Ms)*coeff
  ##   which gives units of aC/(nm*mus)
  ##
  ##   Energies are natively printed in units
  ##    of 0.160218 pg nm^2 / mus^2
  ##    or 1.60218*10^{-22} J
  ##    or 0.001 eV
  ##
  ##############################


Dedef = 3.7551

D0def = 0.003

K1def = -5.0068

K1cdef = -0.0550748

Ktdef = -0.00365997

alphadef = 0.8

#freq = 0.2

corr = 0.0005
dwid = 0.5

hamp = 0.35355339059327373
fDW = 7.88
sDW = 23.53

[Mesh]
  [fileload]
    type = FileMeshGenerator
    file = BFO_dwP1A1_100_ref.e
    use_for_exodus_restart = true
  []
[]


[GlobalParams]
  mag1_x = mag1_x
  mag1_y = mag1_y
  mag1_z = mag1_z

  mag2_x = mag2_x
  mag2_y = mag2_y
  mag2_z = mag2_z

  antiphase_A_x = antiphase_A_x
  antiphase_A_y = antiphase_A_y
  antiphase_A_z = antiphase_A_z

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z


[]

[Materials]
  ############################################################################
  ##==
  ##       material constants used.
  ##
  ##
  ############################################################################

  [./constants]
    type = GenericConstantMaterial
    prop_names = ' alpha           De       D0             g0mu0Ms      g0          K1         K1c      Kt     permittivity Ae     Ms '
    prop_values = '${alphadef}   ${Dedef} ${D0def}      48291.9      48291.9     ${K1def}  ${K1cdef} ${Ktdef}     1.0      0.75    1.0 '
  [../]

  [./a_long]
    type = GenericFunctionMaterial
    prop_names = 'alpha_long'
    prop_values = 'bc_func_1'
  [../]

[]

[Functions]

  ##############################
  ##
  ## Define the ramping function
  ## expression to be used
  ##==
  ##############################

  [./bc_func_1]
    type = ParsedFunction
    expression = 'st'
    symbol_names = 'st'
    symbol_values = '1e8'
  [../]


  [./stripem1x]
    type = ParsedFunction
    expression = 'if(x < 15.707963267948966,(${hamp}+${corr})*tanh(${dwid}*(x-${fDW}))+${hamp},-(${hamp}+${corr})*tanh(${dwid}*(x-${sDW}))+${hamp})'
  [../]

  [./stripem1y]
    type = ParsedFunction
    expression = 'if(x < 15.707963267948966,-(2*${hamp}+${corr})*tanh(${dwid}*(x-${fDW})),(2*${hamp}+${corr})*tanh(${dwid}*(x-${sDW})))'
  [../]
  [./stripem1z]
    type = ParsedFunction
    expression = 'if(x < 15.707963267948966,-(${hamp}+${corr})*tanh(${dwid}*(x-${fDW}))+${hamp},(${hamp}+${corr})*tanh(${dwid}*(x-${sDW}))+${hamp})'
  [../]

  [./stripem2x]
    type = ParsedFunction
    expression = 'if(x < 15.707963267948966,-(${hamp}-${corr})*tanh(${dwid}*(x-${fDW}))-${hamp},(${hamp}-${corr})*tanh(${dwid}*(x-${sDW}))-${hamp})'
  [../]
  [./stripem2y]
    type = ParsedFunction
    expression = 'if(x < 15.707963267948966,(2*${hamp}-${corr})*tanh(${dwid}*(x-${fDW})),-(2*${hamp}-${corr})*tanh(${dwid}*(x-${sDW})))'
  [../]
  [./stripem2z]
    type = ParsedFunction
    expression = 'if(x < 15.707963267948966,(${hamp}-${corr})*tanh(${dwid}*(x-${fDW}))-${hamp},-(${hamp}-${corr})*tanh(${dwid}*(x-${sDW}))-${hamp})'
  [../]
[]
#==

[Variables]
  [./mag1_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = stripem1x
    [../]
  [../]
  [./mag1_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = stripem1y
    [../]
  [../]
  [./mag1_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = stripem1z
    [../]
  [../]

  [./mag2_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = stripem2x
    [../]
  [../]
  [./mag2_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = stripem2y
    [../]
  [../]
  [./mag2_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = stripem2z
    [../]
  [../]
[]

[AuxVariables]

  [./mag1_s]
    order = FIRST
    family = LAGRANGE
  [../]
  [./mag2_s]
    order = FIRST
    family = LAGRANGE
  [../]

  [./Neel_L_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Neel_L_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Neel_L_z]
    order = FIRST
    family = LAGRANGE
  [../]


  [./SSMag_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./SSMag_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./SSMag_z]
    order = FIRST
    family = LAGRANGE
  [../]


  [./antiphase_A_x]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = antiphase_A_x
    initial_from_file_timestep = 'LATEST'
  [../]
  [./antiphase_A_y]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = antiphase_A_y
    initial_from_file_timestep = 'LATEST'
  [../]
  [./antiphase_A_z]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = antiphase_A_z
    initial_from_file_timestep = 'LATEST'
  [../]


  [./polar_x]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_x
    initial_from_file_timestep = 'LATEST'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_y
    initial_from_file_timestep = 'LATEST'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_z
    initial_from_file_timestep = 'LATEST'
  [../]

  [./dSSMag_dt_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dSSMag_dt_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dSSMag_dt_z]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./dL_dt_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dL_dt_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dL_dt_z]
    order = CONSTANT
    family = MONOMIAL
  [../]

  #############################################################
  ##
  ##  other angles
  ##
  ###########################################################

  [./ph]
    order = FIRST
    family = LAGRANGE
  [../]

  [./th1]
    order = FIRST
    family = LAGRANGE
  [../]
  [./th2]
    order = FIRST
    family = LAGRANGE
  [../]

  ######################
  ##                  ##
  ##  Energy density  ##
  ##                  ##
  ######################

  [./Edmi]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Esupexch]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Enlexch]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Eepa1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Eepa2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Eca1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Eca2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Etot]
    order = CONSTANT
    family = MONOMIAL
  [../]

[]

[Kernels]
  #---------------------------------------#
  #                                       #
  #          Time dependence              #
  #                                       #
  #---------------------------------------#

  [./mag1_x_time]
    type = TimeDerivative
    variable = mag1_x
  [../]
  [./mag1_y_time]
    type = TimeDerivative
    variable = mag1_y
  [../]
  [./mag1_z_time]
    type = TimeDerivative
    variable = mag1_z
  [../]

  [./mag2_x_time]
    type = TimeDerivative
    variable = mag2_x
  [../]
  [./mag2_y_time]
    type = TimeDerivative
    variable = mag2_y
  [../]
  [./mag2_z_time]
    type = TimeDerivative
    variable = mag2_z
  [../]


  #---------------------------------------#
  #                                       #
  #     AFM sublattice exchange           #
  #                                       #
  #---------------------------------------#

  [./afmex1_x]
    type = AFMSublatticeSuperexchange
    variable = mag1_x
    mag_sub = 0
    component = 0
  [../]
  [./afmex1_y]
    type = AFMSublatticeSuperexchange
    variable = mag1_y
    mag_sub = 0
    component = 1
  [../]
  [./afmex1_z]
    type = AFMSublatticeSuperexchange
    variable = mag1_z
    mag_sub = 0
    component = 2
  [../]

  [./afmex2_x]
    type = AFMSublatticeSuperexchange
    variable = mag2_x
    mag_sub = 1
    component = 0
  [../]
  [./afmex2_y]
    type = AFMSublatticeSuperexchange
    variable = mag2_y
    mag_sub = 1
    component = 1
  [../]
  [./afmex2_z]
    type = AFMSublatticeSuperexchange
    variable = mag2_z
    mag_sub = 1
    component = 2
  [../]

  #---------------------------------------#
  #                                       #
  #     AFM sublattice DMI                #
  #                                       #
  #---------------------------------------#

  [./afmdmi1_x]
    type = AFMSublatticeDMInteraction
    variable = mag1_x
    mag_sub = 0
    component = 0
  [../]
  [./afmdmi1_y]
    type = AFMSublatticeDMInteraction
    variable = mag1_y
    mag_sub = 0
    component = 1
  [../]
  [./afmdmi1_z]
    type = AFMSublatticeDMInteraction
    variable = mag1_z
    mag_sub = 0
    component = 2
  [../]

  [./afmdmi2_x]
    type = AFMSublatticeDMInteraction
    variable = mag2_x
    mag_sub = 1
    component = 0
  [../]
  [./afmdmi2_y]
    type = AFMSublatticeDMInteraction
    variable = mag2_y
    mag_sub = 1
    component = 1
  [../]
  [./afmdmi2_z]
    type = AFMSublatticeDMInteraction
    variable = mag2_z
    mag_sub = 1
    component = 2
  [../]

  #---------------------------------------#
  #                                       #
  #   Magnetocrystalline anisotropy for   #
  #   the AFM sublattice in easy-plane    #
  #                                       #
  #---------------------------------------#

  [./afma1_x]
    type = AFMEasyPlaneAnisotropy
    variable = mag1_x
    mag_sub = 0
    component = 0
  [../]
  [./afma1_y]
    type = AFMEasyPlaneAnisotropy
    variable = mag1_y
    mag_sub = 0
    component = 1
  [../]
  [./afma1_z]
    type = AFMEasyPlaneAnisotropy
    variable = mag1_z
    mag_sub = 0
    component = 2
  [../]


  [./afma2_x]
    type = AFMEasyPlaneAnisotropy
    variable = mag2_x
    mag_sub = 1
    component = 0
  [../]
  [./afma2_y]
    type = AFMEasyPlaneAnisotropy
    variable = mag2_y
    mag_sub = 1
    component = 1
  [../]
  [./afma2_z]
    type = AFMEasyPlaneAnisotropy
    variable = mag2_z
    mag_sub = 1
    component = 2
  [../]


  #---------------------------------------#
  #                                       #
  #   Single-ion anisotropy environment   #
  #   for the AFM sublattice in the       #
  #   degenerate easy-plane               #
  #                                       #
  #---------------------------------------#

  [./afmsia1_x]
    type = AFMSingleIonCubicSixthAnisotropy
    variable = mag1_x
    mag_sub = 0
    component = 0
  [../]
  [./afmsia1_y]
    type = AFMSingleIonCubicSixthAnisotropy
    variable = mag1_y
    mag_sub = 0
    component = 1
  [../]
  [./afmsia1_z]
    type = AFMSingleIonCubicSixthAnisotropy
    variable = mag1_z
    mag_sub = 0
    component = 2
  [../]

  [./afmsia2_x]
    type = AFMSingleIonCubicSixthAnisotropy
    variable = mag2_x
    mag_sub = 1
    component = 0
  [../]
  [./afmsia2_y]
    type = AFMSingleIonCubicSixthAnisotropy
    variable = mag2_y
    mag_sub = 1
    component = 1
  [../]
  [./afmsia2_z]
    type = AFMSingleIonCubicSixthAnisotropy
    variable = mag2_z
    mag_sub = 1
    component = 2
  [../]

  #---------------------------------------#
  #                                       #
  #          LLB constraint terms         #
  #                                       #
  #---------------------------------------#

  [./llb1_x]
    type = LongitudinalLLB
    variable = mag1_x
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    component = 0
  [../]
  [./llb1_y]
    type = LongitudinalLLB
    variable = mag1_y
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    component = 1
  [../]

  [./llb1_z]
    type = LongitudinalLLB
    variable = mag1_z
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    component = 2
  [../]

  [./llb2_x]
    type = LongitudinalLLB
    variable = mag2_x
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    component = 0
  [../]
  [./llb2_y]
    type = LongitudinalLLB
    variable = mag2_y
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    component = 1
  [../]

  [./llb2_z]
    type = LongitudinalLLB
    variable = mag2_z
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    component = 2
  [../]

  #---------------------------------------#
  #                                       #
  #    Local magnetic exchange            #
  #                                       #
  #---------------------------------------#

  [./dllg1_x_exch]
    type = ExchangeCartLL
    variable = mag1_x
    component = 0
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
  [../]
  [./dllg1_y_exch]
    type = ExchangeCartLL
    variable = mag1_y
    component = 1
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
  [../]
  [./dllg1_z_exch]
    type = ExchangeCartLL
    variable = mag1_z
    component = 2
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
  [../]

  [./dllg2_x_exch]
    type = ExchangeCartLL
    variable = mag2_x
    component = 0
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
  [../]
  [./dllg2_y_exch]
    type = ExchangeCartLL
    variable = mag2_y
    component = 1
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
  [../]
  [./dllg2_z_exch]
    type = ExchangeCartLL
    variable = mag2_z
    component = 2
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
  [../]

  #---------------------------------------#
  #                                       #=
  #    Local AFM magnetic exchange        #
  #           (Local-term)                #
  #---------------------------------------#

  [./dllg1afm_x_exch]
    type = AFMLocalSublatticeExchangeCartLL
    variable = mag1_x
    mag_sub = 0
    component = 0
  [../]
  [./dllg1afm_y_exch]
    type = AFMLocalSublatticeExchangeCartLL
    variable = mag1_y
    mag_sub = 0
    component = 1
  [../]
  [./dllg1afm_z_exch]
    type = AFMLocalSublatticeExchangeCartLL
    variable = mag1_z
    mag_sub = 0
    component = 2
  [../]

  [./dllg2afm_x_exch]
    type = AFMLocalSublatticeExchangeCartLL
    variable = mag2_x
    mag_sub = 1
    component = 0
  [../]
  [./dllg2afm_y_exch]
    type = AFMLocalSublatticeExchangeCartLL
    variable = mag2_y
    mag_sub = 1
    component = 1
  [../]
  [./dllg2afm_z_exch]
    type = AFMLocalSublatticeExchangeCartLL
    variable = mag2_z
    mag_sub = 1
    component = 2
  [../]
[]



[AuxKernels]
  [./mag1_mag]
    type = VectorMag
    variable = mag1_s
    vector_x = mag1_x
    vector_y = mag1_y
    vector_z = mag1_z
    execute_on = 'initial timestep_end final'
  [../]

  [./mag2_mag]
    type = VectorMag
    variable = mag2_s
    vector_x = mag2_x
    vector_y = mag2_y
    vector_z = mag2_z
    execute_on = 'initial timestep_end final'
  [../]


  [./Neel_Lx]
    type = VectorDiffOrSum
    variable = Neel_L_x
    var1 = mag1_x
    var2 = mag2_x
    diffOrSum = 0
    execute_on = 'initial timestep_end final'
  [../]
  [./Neel_Ly]
    type = VectorDiffOrSum
    variable = Neel_L_y
    var1 = mag1_y
    var2 = mag2_y
    diffOrSum = 0
    execute_on = 'initial timestep_end final'
  [../]
  [./Neel_Lz]
    type = VectorDiffOrSum
    variable = Neel_L_z
    var1 = mag1_z
    var2 = mag2_z
    diffOrSum = 0
    execute_on = 'initial timestep_end final'
  [../]

  [./smallSignalMag_x]
    type = VectorDiffOrSum
    variable = SSMag_x
    var1 = mag1_x
    var2 = mag2_x
    diffOrSum = 1
    execute_on = 'initial timestep_end final'
  [../]
  [./smallSignalMag_y]
    type = VectorDiffOrSum
    variable = SSMag_y
    var1 = mag1_y
    var2 = mag2_y
    diffOrSum = 1
    execute_on = 'initial timestep_end final'
  [../]
  [./smallSignalMag_z]
    type = VectorDiffOrSum
    variable = SSMag_z
    var1 = mag1_z
    var2 = mag2_z
    diffOrSum = 1
    execute_on = 'initial timestep_end final'
  [../]

  [./phc]
    type = AngleBetweenTwoVectors
    variable = ph
    var1x = mag1_x
    var1y = mag1_y
    var1z = mag1_z
    var2x = mag2_x
    var2y = mag2_y
    var2z = mag2_z

    execute_on = 'initial timestep_end final'
  [../]

  [./th1c]
    type = AngleBetweenTwoVectors
    variable = th1
    var1x = mag1_x
    var1y = mag1_y
    var1z = mag1_z
    var2x = polar_x
    var2y = polar_y
    var2z = polar_z

    execute_on = 'initial timestep_end final'
  [../]

  [./th2c]
    type = AngleBetweenTwoVectors
    variable = th2
    var1x = mag2_x
    var1y = mag2_y
    var1z = mag2_z
    var2x = polar_x
    var2y = polar_y
    var2z = polar_z

    execute_on = 'initial timestep_end final'
  [../]

  #---------------------------------------#
  #                                       #
  #          Energy density               #
  #                                       #
  #---------------------------------------#

  [./cEdmi]
    type = AFMSublatticeDMInteractionEnergyDensity
    variable = Edmi
    execute_on = 'initial timestep_end final'
    energy_scale = -6241.51
  [../]
  [./cEsupexch]
    type = AFMSublatticeSuperexchangeEnergyDensity
    variable = Esupexch
    execute_on = 'initial timestep_end final'
    energy_scale = 6241.51
  [../]
  [./cEnlexch]
    type = AFMExchangeStiffnessEnergyDensity
    variable = Enlexch
    Neel_L_x = Neel_L_x
    Neel_L_y = Neel_L_y
    Neel_L_z = Neel_L_z
    execute_on = 'initial timestep_end final'
    energy_scale = 6241.51
  [../]

  [./cEepa1]
    type = AFMEasyPlaneAnisotropyEnergyDensity
    variable = Eepa1
    execute_on = 'initial timestep_end final'
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    energy_scale = 6241.51
  [../]
  [./cEepa2]
    type = AFMEasyPlaneAnisotropyEnergyDensity
    variable = Eepa2
    execute_on = 'initial timestep_end final'
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    energy_scale = 6241.51
  [../]
  [./cEca1]
    type = AFMSingleIonCubicSixthAnisotropyEnergyDensity
    variable = Eca1
    execute_on = 'initial timestep_end final'
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    energy_scale = 6241.51
  [../]
  [./cEca2]
    type = AFMSingleIonCubicSixthAnisotropyEnergyDensity
    variable = Eca2
    execute_on = 'initial timestep_end final'
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    energy_scale = 6241.51
  [../]

  [./cEdtot]
    type = AFMTotalEnergyDensity
    variable = Etot
    execute_on = 'initial timestep_end final'
    Edmi = Edmi
    Esupexch = Esupexch
    Enlexch = Enlexch
    Eepa1 = Eepa1
    Eepa2 = Eepa2
    Eca1 = Eca1
    Eca2 = Eca2
  [../]

[]


[BCs]
  #---------------------------------------#
  #                                       #
  #  periodic magnetization distribution  #
  #                                       #
  #---------------------------------------#

  [./Periodic]
    [./x]
      auto_direction = 'x'
      variable = 'mag1_x mag1_y mag1_z mag2_x mag2_y mag2_z'
    [../]
  [../]

[]

[Postprocessors]
   [./dt]
     type = TimestepSize
   [../]

  #---------------------------------------#
  #                                       #
  #     Average Mk = |m_k| and along      #
  #     other directions                  #
  #                                       #
  #---------------------------------------#

  [./M1]
    type = ElementAverageValue
    variable = mag1_s
    execute_on = 'initial timestep_end final'
  [../]
  [./M2]
    type = ElementAverageValue
    variable = mag2_s
    execute_on = 'initial timestep_end final'
  [../]

  [./<Lx>]
    type = ElementAverageValue
    variable = Neel_L_x
    execute_on = 'initial timestep_end final'
  [../]
  [./<Ly>]
    type = ElementAverageValue
    variable = Neel_L_y
    execute_on = 'initial timestep_end final'
  [../]
  [./<Lz>]
    type = ElementAverageValue
    variable = Neel_L_z
    execute_on = 'initial timestep_end final'
  [../]

  [./<m1x>]
    type = ElementAverageValue
    variable = mag1_x
    execute_on = 'initial timestep_end final'
  [../]
  [./<m1y>]
    type = ElementAverageValue
    variable = mag1_y
    execute_on = 'initial timestep_end final'
  [../]
  [./<m1z>]
    type = ElementAverageValue
    variable = mag1_z
    execute_on = 'initial timestep_end final'
  [../]

  [./<m2x>]
    type = ElementAverageValue
    variable = mag2_x
    execute_on = 'initial timestep_end final'
  [../]
  [./<m2y>]
    type = ElementAverageValue
    variable = mag2_y
    execute_on = 'initial timestep_end final'
  [../]
  [./<m2z>]
    type = ElementAverageValue
    variable = mag2_z
    execute_on = 'initial timestep_end final'
  [../]

  [./<ph>]
    type = ElementAverageValue
    variable = ph
    execute_on = 'initial timestep_end final'
  [../]
  [./<th1>]
    type = ElementAverageValue
    variable = th1
    execute_on = 'initial timestep_end final'
  [../]
  [./<th2>]
    type = ElementAverageValue
    variable = th2
    execute_on = 'initial timestep_end final'
  [../]

  #---------------------------------------#
  #                                       #
  #   Calculate exchange energy of        #
  #   the magnetic body                   #
  #                                       #
  #---------------------------------------#

  [./FafmSLexch]
    type = AFMSublatticeSuperexchangeEnergy
    execute_on = 'initial timestep_end final'
    mag1_x = mag1_x
    mag1_y = mag1_y
    mag1_z = mag1_z
    mag2_x = mag2_x
    mag2_y = mag2_y
    mag2_z = mag2_z
    energy_scale = 0.001

  [../]
  [./FafmSLdmi]
    type = AFMSublatticeDMInteractionEnergy
    execute_on = 'initial timestep_end final'
    energy_scale = 0.001
  [../]

  [./FexchS]
    type = AFMExchangeStiffnessEnergy
    execute_on = 'initial timestep_end final'
    energy_scale = 0.001
    Neel_L_x = Neel_L_x
    Neel_L_y = Neel_L_y
    Neel_L_z = Neel_L_z
  [../]

  #---------------------------------------#
  #                                       #
  #   Calculate excess energy from missed #
  #   LLB targets                         #
  #                                       #
  #---------------------------------------#

  [./Fllb1]
    type = MagneticExcessLLBEnergy
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    execute_on = 'initial timestep_end final'
  [../]
  [./Fllb2]
    type = MagneticExcessLLBEnergy
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    execute_on = 'initial timestep_end final'
  [../]

  #---------------------------------------#
  #                                       #
  #   Calculate the anisotropy energy     #
  #                                       #
  #---------------------------------------#

  [./Fa1]
    type = AFMEasyPlaneAnisotropyEnergy
    execute_on = 'initial timestep_end final'
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    energy_scale = 0.001
  [../]
  [./Fa2]
    type = AFMEasyPlaneAnisotropyEnergy
    execute_on = 'initial timestep_end final'
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    energy_scale = 0.001
  [../]

  [./Fsia1]
    type = AFMSingleIonCubicSixthAnisotropyEnergy
    execute_on = 'initial timestep_end final'
    mag_x = mag1_x
    mag_y = mag1_y
    mag_z = mag1_z
    energy_scale = 0.001
  [../]
  [./Fsia2]
    type = AFMSingleIonCubicSixthAnisotropyEnergy
    execute_on = 'initial timestep_end final'
    mag_x = mag2_x
    mag_y = mag2_y
    mag_z = mag2_z
    energy_scale = 0.001
  [../]

  #---------------------------------------#
  #                                       #
  #   add all the energy contributions    #
  #   and calculate their percent change  #
  #                                       #
  #---------------------------------------#

  [./Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'FafmSLexch FafmSLdmi Fa1 Fa2 Fsia1 Fsia2 FexchS'
    pp_coefs = ' 1.0 1.0 1.0 1.0 1.0 1.0 1.0'
    execute_on = 'initial timestep_end final'
  [../]

  [./FtotLLB]
    type = LinearCombinationPostprocessor
    pp_names = 'Fllb1 Fllb2'
    pp_coefs = ' 1.0 1.0'
    execute_on = 'initial timestep_end final'
  [../]

  [./perc_change]
    type = EnergyRatePostprocessor
    postprocessor = Ftot
    dt = dt
    execute_on = 'timestep_end final'
  [../]

  [./elapsed]
    type = PerfGraphData
    section_name = "Root"  # for profiling the problem
    data_type = total
  [../]
[]

[UserObjects]
  [./kill]
    type = Terminator
    expression = 'perc_change <= 1.0e-8'
  [../]
[]

[Preconditioning]

  #---------------------------------------#
  #                                       #
  #            Solver options             #
  #                                       #
  #---------------------------------------#

  [./smp]
    type = SMP
    full = true
    petsc_options_iname = ' -ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type '
    petsc_options_value = '    250               1e-8      1e-8      1e-5     bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'

  [./TimeIntegrator]
    type = NewmarkBeta
  [../]

  dtmin = 1e-18
  dtmax = 1e-7

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 25  #usually 10
    linear_iteration_ratio = 100
    dt = 1e-9
    growth_factor = 1.1
  [../]

  num_steps = 18001
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_BFO_P111bA111b-P111bA111b_m1
    elemental_as_nodal = true
    time_step_interval = 20
  [../]
  [./outCSV]
    type = CSV
    file_base = out_BFO_P111bA111b-P111bA111b_m1
  [../]
[]
