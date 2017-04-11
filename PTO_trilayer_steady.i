
[Mesh]
  file = exodus_trilayer_12_12_8.e
[]

[GlobalParams]
  #NOTE: We use a nanometer, nanonewton, attocoulomb unit system

  len_scale = 1.0
  alpha1 = -0.1722883 # T = 298 K
  alpha11 = -0.07253
  alpha111 = 0.26
  alpha12 = 0.75
  alpha112 = 0.61
  alpha123 = -3.67

  G110 = 0.173
  G11/G110 = 0.6
  G12/G110 = 0
  G44/G110 = 0.3
  G44P/G110 = 0.3

  #miscellaneous GlobalParams options
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int = potential_int

  epsilon = 0.008 #negative = tension, positive = compression
  T = 298
[]

[Variables]
  [./polar_x]
    block = '1'
    order = FIRST
    family = LAGRANGE
  [../]
  [./polar_y]
    block = '1'
    order = FIRST
    family = LAGRANGE
  [../]
  [./polar_z]
    block = '1'
    order = FIRST
    family = LAGRANGE
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
  [../]
[]


[AuxVariables]
  #semiconducting charge carriers (store their values)
  [./nm]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pp]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./NAm]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rho]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]


[AuxKernels]
  
  #calculate the semiconducting charge carriers
  [./nm_calc]
    type = SemiconductingChargeCarriersAux
    charge_type = 0
    variable = nm
    kT = 0.0041124
    q = -0.16
    NA = 1e-8 #1e20 / m^3
    NC = 0.0000666
    NV = 0.0015
    EA = -0.00801
    EC = -0.66483
    EV = -0.71289 #-4.45 eV
    EF = -0.71289 #~-4.45 eV
    block  = '2'
    execute_on = 'timestep_end'
  [../]
  [./pp_calc]
    type = SemiconductingChargeCarriersAux
    charge_type = 1
    variable = pp
    kT = 0.0041124
    q = -0.16
    NA = 1e-8 #1e20 / m^3
    NC = 0.0000666
    NV = 0.0015
    EA = -0.00801
    EC = -0.66483
    EV = -0.71289 #-4.45 eV
    EF = -0.71289 #~-4.45 eV
    block  = '2'
    execute_on = 'timestep_end'
  [../]
  [./NAm_calc]
    type = SemiconductingChargeCarriersAux
    charge_type = 2
    variable = NAm
    kT = 0.0041124
    q = -0.16
    NA = 1e-8 #1e20 / m^3
    NC = 0.0000666
    NV = 0.0015
    EA = -0.00801
    EC = -0.66483
    EV = -0.71289 #-4.45 eV
    EF = -0.71289 #~-4.45 eV
    block  = '2'
    execute_on = 'timestep_end'
  [../]
  [./rho_calc]
    type = SemiconductingChargeCarriersAux
    charge_type = 3
    variable = rho
    kT = 0.0041124
    q = -0.16
    NA = 1e-8 #1e20 / m^3
    NC = 0.0000666
    NV = 0.0015
    EA = -0.00801
    EC = -0.66483
    EV = -0.71289 #-4.45 eV
    EF = -0.71289 #~-4.45 eV
    block  = '2'
    execute_on = 'timestep_end'
  [../]
[]


[Kernels]
  #Bulk energy density
  [./bed_x]
    type = RenormalizedFreeEnergy
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = RenormalizedFreeEnergy
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = RenormalizedFreeEnergy
    variable = polar_z
    component = 2
  [../]
  ##Wall energy penalty
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

  ##Electrostatics
  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
     #permittivity = 0.0575522155 #er = 6.5
     block = '1'
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     block = '1'
     permittivity = 0.0575522155 #er = 6.5
  [../]

  [./semi_E_int]
     type = Electrostatics
     variable = potential_int
     block  = '2'
     permittivity = 0.44270935 #er = 50
  [../]

  [./semiconducting_charge_carriers]
     type = SemiconductorChargeCarriers
     variable = potential_int
     kT = 0.0041124
     q = -0.16
     NA = 1e-9 #1e20 / m^3 1e-7 fact of 10..
     NC = 0.00000666
     NV = 0.00015
     EA = -0.00801
     EC = -0.66483
     EV = -0.71289 #-4.45 eV
     EF = -0.71289 #~-4.45 eV
     block  = '2'
     #note: can get elast constants for Sb2Te3 from
     #      Comp. Mater. Sci. 96, 342–347, (2015)
  [../]
  #
  [./FE_charge_carriers]
    type = SemiconductorChargeCarriers
    variable = potential_int
    kT = 0.0041124
     q = 0.16
     NA = 0.0 #1e20 / m^3
     NC = 0.00025
     NV = 0.00025
     EA = 0.0
     EC = -0.5607
     EV = -1.10538 #-4.45 eV
     EF = -0.92916 #~-5.8 eV
     block  = '1'
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

[]


[BCs]

  [./bot_potential_int]
    variable = potential_int
    type = DirichletBC
    value = 0.0
    boundary = '1'
  [../]

  [./top_potential_int]
    variable = potential_int
    type = DirichletBC
    value = 0.0
    boundary = '2'
  [../]
[]


[Postprocessors]
   [./Fbulk]
      type = BulkEnergy
      block = '1'
      execute_on = 'timestep_end'
    [../]
    [./Fwall]
      type = WallEnergy
      block = '1'
      execute_on = 'timestep_end'
    [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason -ksp_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart  -snes_atol -snes_rtol -pc_type'
    petsc_options_value = '       200            1e-12     1e-10         lu'
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = outPTO_trilayer_steady
    elemental_as_nodal = true
    interval = 1
  [../]
[]
