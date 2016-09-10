[Mesh]
  file = mie_test.e
  #type = GeneratedMesh
  #dim = 3
  #nx = 50
  #ny = 50
  #nz = 50
  #xmin = -2
  #xmax = 2
  #ymin = -2
  #ymax = 2
  #zmin = -2
  #zmax = 2
  #elem_type = HEX8
[]

[GlobalParams]
  a = 1.0
  nh =  1.458
  omega = 1.0 #These aren't used at the moment.
  c = 1.0
  epsilonI = 1.0 
  sigmaI = 0.0
  epsilonII = 1.0
  sigmaII = 0.0
  L = 1.5
  order = 20 #IIRC convergence should be around ~20.
  scale = 1.0 #possibly not needed.
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./ReE_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./ReE_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./ReE_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./ImagE_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./ImagE_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./ImagE_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./intensity]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./mie_field_ReE_x]
    type = MieElecFieldReals
    variable = ReE_x
    execute_on = 'timestep_end'
    component = 0
    block = '2'
  [../]
  [./mie_field_ReE_y]
    type = MieElecFieldReals
    variable = ReE_y
    execute_on = 'timestep_end'
    component = 1
    block = '2'
  [../]
  [./mie_field_ReE_z]
    type = MieElecFieldReals
    variable = ReE_z
    execute_on = 'timestep_end'
    component = 2
    block = '2'
  [../]

  [./mie_field_ImagE_x]
    type = MieElecFieldImag
    variable = ImagE_x
    execute_on = 'timestep_end'
    component = 0
    block = '2'
  [../]
  [./mie_field_ImagE_y]
    type = MieElecFieldImag
    variable = ImagE_y
    execute_on = 'timestep_end'
    component = 1
    block = '2'
  [../]
  [./mie_field_ImagE_z]
    type = MieElecFieldImag
    variable = ImagE_z
    execute_on = 'timestep_end'
    component = 2
    block = '2'
  [../]

  [./mie_field_intensity]
    type = Intensity
    variable = intensity
    ReE_x = ReE_x
    ReE_y = ReE_y
    ReE_z = ReE_z
    ImagE_x = ImagE_x
    ImagE_y = ImagE_y
    ImagE_z = ImagE_z
    execute_on = 'timestep_end'
    block = '2'
  [../]
[]

[Kernels]
  [./u]
    type = Diffusion
    variable = u
  [../]
[]

[BCs]
  [./top]
    type = DirichletBC
    variable = u
    boundary = '1'
    value = 0.0
  [../]
  [./bot]
    type = DirichletBC
    variable = u
    boundary = '2'
    value = 1
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-snes_rtol -ksp_rtol'
    petsc_options_value = '  1e-1      1e-1 '
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
[]

[Outputs]
  print_perf_log = true
  [./out]
    type = Exodus
    elemental_as_nodal = true
    execute_on = 'timestep_end'
  [../]
[]

