

val = -0.2652339 #note that this quantity depends on the amplitude of the applied field
freq = 2.5e10
amplitude = 0.01

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

[GlobalParams]

  len_scale = 1.0

  #########################################
  ##
  ## Gradient and Landau coefficients from
  ## Marton and Hlinka
  ##    Phys. Rev. B. 74, 104014, (2006)
  ##
  #########################################

  alpha1 = -0.02772
  alpha11 = -0.6476
  alpha111 = 8.004
  alpha12 = 0.323
  alpha112 = 4.47
  alpha123 = 4.919

  G110 = 1.0
  G11_G110 = 0.51
  G12_G110 = 0.02
  G44_G110 = 0.02
  G44P_G110 = 0.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_E_int = potential_int
[]


[Variables]

  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-5
      max = 0.1e-5
      legacy_generator = true
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-5
      max = 0.1e-5
      legacy_generator = true
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = ${val}
      max = ${val}
      legacy_generator = true
    [../]
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./Ps] # Just a placer variable to make postprocessing easier
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = RandomIC
      min = ${val}
      max = ${val}
      legacy_generator = true
    [../]
  [../]
  [./Ez]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./cEz]
    type = ElecFieldAux
    component = 2
    variable = Ez
  [../]
[]

[Kernels]

  #########################################
  ##
  ## Landau's problem (no elastic coupling:
  ##
  #########################################

  [./bed_x]
    type = BulkEnergyDerivativeSixth
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = BulkEnergyDerivativeSixth
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = BulkEnergyDerivativeSixth
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

  #########################################
  ##
  ## Poisson's equation and P*E interaction
  ##
  #########################################

  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     permittivity = 0.0637501464
     #NOTE: This is a static permittivity contribution from core-electrons.
     #      This effectively screens the electrostatic interactions.
     #      For BTO, this value is about 7*e_0, where e_0 is the permitivitty of the vacuum.
     #      See the brief discussion before sec. IV on pp. 4 of Phys. Rev. B. 74, 104014, (2006)
  [../]
  [./polar_electric_px]
     type = PolarElectricPStrong
     variable = polar_x
     component = 0
  [../]
  [./polar_electric_py]
     type = PolarElectricPStrong
     variable = polar_y
     component = 1
  [../]
  [./polar_electric_pz]
     type = PolarElectricPStrong
     variable = polar_z
     component = 2
  [../]

  #########################################
  ##
  ## Time dependence
  ##
  #########################################

  [./polar_x_time]
     type = TimeDerivativeScaled
     variable = polar_x
     # Time scale estimate for BTO, from Hlinka (2007)
     # We use seconds here
     time_scale = 1e-12
  [../]
  [./polar_y_time]
     type = TimeDerivativeScaled
     variable = polar_y
     time_scale = 1e-12
  [../]
  [./polar_z_time]
     type = TimeDerivativeScaled
     variable = polar_z
     time_scale = 1e-12
  [../]
[]

[Functions]

  ##############################
  ##
  ## Define the electric field
  ## expression to be used below
  ##
  ##############################

  [./bc_func_1]
    type = ParsedFunction
    value = 'amplitude*sin(freq*t)'
    vars = 'freq amplitude'
    vals = '${freq}  ${amplitude}'
  [../]
[]

[BCs]

  ##############################
  ##
  ## Boundary Condition System
  ##
  ## This corresponds to a small
  ## ac field along the direction
  ## of the spontaneous polarization
  ## that was ostensibly put there
  ## by a strong dc bias.
  ##
  ##############################

  [./front_pot]
    type = FunctionDirichletBC
    variable = potential_int
    boundary = 'front'
    function = bc_func_1
  [../]
  [./back_pot]
    type = DirichletBC
    variable = potential_int
    boundary = 'back'
    value = 0.0
  [../]
[]

[Postprocessors]
  [./avePz]
    type = ElementAverageValue
    variable = polar_z
    execute_on = 'initial timestep_end'
  [../]
  [./cPs]
    type = ElementAverageValue
    variable = Ps
    execute_on = 'initial timestep_end'
  [../]
  [./Ea]
    type = ElementAverageValue
    variable = Ez
    execute_on = 'initial timestep_end'
  [../]
  [./inducedP]
    type = LinearCombinationPostprocessor
    pp_names = 'avePz cPs'
    pp_coefs = ' 1 -1'
    execute_on = 'initial timestep_end'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart  -snes_atol -ksp_rtol -pc_type'
    petsc_options_value = '    121                1e-10     1e-8     bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  scheme = 'implicit-euler'
  dtmin = 1e-20
  dt = 1e-12
  dtmax = 1.0
  num_steps = 15
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = outBTO_0_adef
    elemental_as_nodal = true
    #execute_on = 'initial final'
  [../]
[]
