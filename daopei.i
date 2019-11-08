

val1 = -0.2652339 #note that this quantity depends on the amplitude of the applied field
val2 = -0.2652338
freq = 200000000.0
amplitude = 1e-05

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 4
  ny = 4
  nz = 4
  xmin = -1
  xmax = 1
  ymin = -1
  ymax = 1
  zmin = -1
  zmax = 1
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

  alpha111 = 3.96

  alpha112 = 4.47
  alpha123 = 4.919
  alpha1111 = 0.0
  alpha1112 = 0.0
  alpha1122 = 0.0
  alpha1123 = 0.0

  G110 = 1.0
  G11_G110 = 0.51
  G12_G110 = -0.02
  G44_G110 = 0.02
  G44P_G110 = 0.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  alpha1 = alpha1
  alpha11 = alpha11
  alpha12 = alpha12
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
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-5
      max = 0.1e-5
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = ${val1}
      max = ${val2}
    [../]
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./alpha1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.0005
      max = 0.0022
      seed = 1
    [../]
  [../]
  [./alpha11]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 8.3
      max = 8.8
      seed = 1
    [../]
  [../]
  [./alpha12]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = 0.30
      max = 0.36
      seed = 1
    [../]
  [../]
  [./Ps] # Just a placer variable to make postprocessing easier
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = RandomIC
      min = ${val1}
      max = ${val2}
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
    type = SDBulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = SDBulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = SDBulkEnergyDerivativeEighth
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
  dt = 1.57079632679e-10
  dtmax = 1.0
  num_steps = 20000
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = frequency=2.0e08
    elemental_as_nodal = true
    execute_on = 'initial final'
  [../]
  [./outCSV]
    type = CSV
    new_row_tolerance = 1e-16
    file_base = frequency=2.0e08
  [../]
[]
