[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  nx = 10
  ny = 10
  elem_type = QUAD4
[]


[GlobalParams]
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  antiferrodis_A_x = antiferrodis_A_x
  antiferrodis_A_y = antiferrodis_A_y
  antiferrodis_A_z = antiferrodis_A_z

  e_xx = e_xx
  e_yy = e_yy
  e_zz = e_zz
  e_xy = e_xy
  e_yz = e_yz
  e_zx = e_zx

  Ap = -1.09028
  Bp = 6.68020e-1
  Cp = -4.67264e-1

  Ar = -5.33957e-3
  Br = 1.35373e-5
  Cr = -7.73973e-6

  Bpr = 6.99618e-3
  Cpr = -2.14508e-3
  Cprp = -1.40134e-2

  C11 = 114.409
  C12 = 45.5681
  C44 = 28.709


  ##############################################################
  ##
  ## NOTE: Sign convention in **this implementation**
  ##       for the electrostrictive coeff. is multiplied by
  ##       an overall factor of (-1). Note that other elastic
  ##       coupling Kernels/Materials in Ferret DO NOT have the 
  ##       (-1) prefactor. Please be careful here.
  ##
  ###############################################################
  
  q11 = -5.89457
  q12 = -0.971881
  q44 = -2.01751

  r11 = -7.35454e-3
  r12 = 7.23046e-4
  r44 = 7.21126e-3

  G = 1.0
[]


[Variables]
  [./diffused]
    order = FIRST
    family = LAGRANGE
  [../]


  # ODE variables
  [./polar_x]
    family = SCALAR
    order = FIRST
    initial_condition = 0.4
  [../]
  [./polar_y]
    family = SCALAR
    order = FIRST
    initial_condition = 0.4
  [../]
  [./polar_z]
    family = SCALAR
    order = FIRST
    initial_condition = 0.4
  [../]
  [./antiferrodis_A_x]
    family = SCALAR
    order = FIRST
    initial_condition = 7.0
  [../]
  [./antiferrodis_A_y]
    family = SCALAR
    order = FIRST
    initial_condition = 7.0
  [../]
  [./antiferrodis_A_z]
    family = SCALAR
    order = FIRST
    initial_condition = 7.0
  [../]

  [./e_xx]
    family = SCALAR
    order = FIRST
    initial_condition = 0.001
  [../]
  [./e_yy]
    family = SCALAR
    order = FIRST
    initial_condition = 0.001
  [../]
  [./e_zz]
    family = SCALAR
    order = FIRST
    initial_condition = 0.001
  [../]

  [./e_xy]
    family = SCALAR
    order = FIRST
    initial_condition = 0.0001
  [../]
  [./e_yz]
    family = SCALAR
    order = FIRST
    initial_condition = 0.0001
  [../]
  [./e_zx]
    family = SCALAR
    order = FIRST
    initial_condition = 0.0001
  [../]


[]

[Kernels]
  [./td]
    type = TimeDerivative
    variable = diffused
  [../]
  [./diff]
    type = Diffusion
    variable = diffused
  [../]
[]

[ScalarKernels]
  [./td_px]
    type = ODETimeDerivative
    variable = polar_x
  [../]
  [./td_py]
    type = ODETimeDerivative
    variable = polar_y
  [../]
  [./td_pz]
    type = ODETimeDerivative
    variable = polar_z
  [../]

  [./td_Ax]
    type = ODETimeDerivative
    variable = antiferrodis_A_x
  [../]
  [./td_Ay]
    type = ODETimeDerivative
    variable = antiferrodis_A_y
  [../]
  [./td_Az]
    type = ODETimeDerivative
    variable = antiferrodis_A_z
  [../]

  [./td_exx]
    type = ODETimeDerivative
    variable = e_xx
  [../]
  [./td_eyy]
    type = ODETimeDerivative
    variable = e_yy
  [../]
  [./td_ezz]
    type = ODETimeDerivative
    variable = e_zz
  [../]

  [./td_exy]
    type = ODETimeDerivative
    variable = e_xy
  [../]
  [./td_eyz]
    type = ODETimeDerivative
    variable = e_yz
  [../]
  [./td_ezx]
    type = ODETimeDerivative
    variable = e_zx
  [../]


  [./bulk_px]
    type = ScalarBulkEnergyP
    variable = polar_x
    component = 0
  [../]
  [./bulk_py]
    type = ScalarBulkEnergyP
    variable = polar_y
    component = 1
  [../]
  [./bulk_pz]
    type = ScalarBulkEnergyP
    variable = polar_z
    component = 2
  [../]

  [./bulk_Ax]
    type = ScalarBulkEnergyA
    variable = antiferrodis_A_x
    component = 0
  [../]
  [./bulk_Ay]
    type = ScalarBulkEnergyA
    variable = antiferrodis_A_y
    component = 1
  [../]
  [./bulk_Az]
    type = ScalarBulkEnergyA
    variable = antiferrodis_A_z
    component = 2
  [../]


  [./roto_Px]
    type = ScalarRotopolarEnergy
    variable = polar_x
    component = 0
  [../]
  [./roto_Py]
    type = ScalarRotopolarEnergy
    variable = polar_y
    component = 1
  [../]
  [./roto_Pz]
    type = ScalarRotopolarEnergy
    variable = polar_z
    component = 2
  [../]
  [./roto_Ax]
    type = ScalarRotopolarEnergy
    variable = antiferrodis_A_x
    component = 3
  [../]
  [./roto_Ay]
    type = ScalarRotopolarEnergy
    variable = antiferrodis_A_y
    component = 4
  [../]
  [./roto_Az]
    type = ScalarRotopolarEnergy
    variable = antiferrodis_A_z
    component = 5
  [../]

  [./estric_Px]
    type = ScalarElectrostrictiveEnergy
    variable = polar_x
    component = 0
  [../]
  [./estric_Py]
    type = ScalarElectrostrictiveEnergy
    variable = polar_y
    component = 1
  [../]
  [./estric_Pz]
    type = ScalarElectrostrictiveEnergy
    variable = polar_z
    component = 2
  [../]

  [./estric_exx]
    type = ScalarElectrostrictiveEnergy
    variable = e_xx
    component = 3
  [../]
  [./estric_eyy]
    type = ScalarElectrostrictiveEnergy
    variable = e_yy
    component = 4
  [../]
  [./estric_ezz]
    type = ScalarElectrostrictiveEnergy
    variable = e_zz
    component = 5
  [../]
  [./estric_exy]
    type = ScalarElectrostrictiveEnergy
    variable = e_xy
    component = 6
  [../]
  [./estric_eyz]
    type = ScalarElectrostrictiveEnergy
    variable = e_yz
    component = 7
  [../]
  [./estric_ezx]
    type = ScalarElectrostrictiveEnergy
    variable = e_zx
    component = 8
  [../]


  [./rstric_Ax]
    type = ScalarRotostrictiveEnergy
    variable = antiferrodis_A_x
    component = 0
  [../]
  [./rstric_Ay]
    type = ScalarRotostrictiveEnergy
    variable = antiferrodis_A_y
    component = 1
  [../]
  [./rstric_Az]
    type = ScalarRotostrictiveEnergy
    variable = antiferrodis_A_z
    component = 2
  [../]

  [./rstric_exx]
    type = ScalarRotostrictiveEnergy
    variable = e_xx
    component = 3
  [../]
  [./rstric_eyy]
    type = ScalarRotostrictiveEnergy
    variable = e_yy
    component = 4
  [../]
  [./rstric_ezz]
    type = ScalarRotostrictiveEnergy
    variable = e_zz
    component = 5
  [../]
  [./rstric_exy]
    type = ScalarRotostrictiveEnergy
    variable = e_xy
    component = 6
  [../]
  [./rstric_eyz]
    type = ScalarRotostrictiveEnergy
    variable = e_yz
    component = 7
  [../]
  [./rstric_ezx]
    type = ScalarRotostrictiveEnergy
    variable = e_zx
    component = 8
  [../]

  [./elast_exx]
    type = ScalarElasticEnergy
    variable = e_xx
    component = 0
  [../]
  [./elast_eyy]
    type = ScalarElasticEnergy
    variable = e_yy
    component = 1
  [../]
  [./elast_ezz]
    type = ScalarElasticEnergy
    variable = e_zz
    component = 2
  [../]
  [./elast_exy]
    type = ScalarElasticEnergy
    variable = e_xy
    component = 3
  [../]
  [./elast_eyz]
    type = ScalarElasticEnergy
    variable = e_yz
    component = 4
  [../]
  [./elast_exz]
    type = ScalarElasticEnergy
    variable = e_zx
    component = 5
  [../]

[]


[BCs]
  [./right]
    type = DirichletBC
    variable = diffused
    boundary = 1
    value = 0.0
  [../]

  [./left]
    type = DirichletBC
    variable = diffused
    boundary = 0
    value = 1.0
  [../]
[]

[Postprocessors]
  [./Px]
    type = ScalarVariable
    variable = polar_x
    execute_on = timestep_end
  [../]
  [./Py]
    type = ScalarVariable
    variable = polar_y
    execute_on = timestep_end
  [../]
  [./Pz]
    type = ScalarVariable
    variable = polar_z
    execute_on = timestep_end
  [../]

  [./Ax]
    type = ScalarVariable
    variable = antiferrodis_A_x
    execute_on = timestep_end
  [../]
  [./Ay]
    type = ScalarVariable
    variable = antiferrodis_A_y
    execute_on = timestep_end
  [../]
  [./Az]
    type = ScalarVariable
    variable = antiferrodis_A_z
    execute_on = timestep_end
  [../]


  [./exx]
    type = ScalarVariable
    variable = e_xx
    execute_on = timestep_end
  [../]
  [./eyy]
    type = ScalarVariable
    variable = e_yy
    execute_on = timestep_end
  [../]
  [./ezz]
    type = ScalarVariable
    variable = e_zz
    execute_on = timestep_end
  [../]

  [./exy]
    type = ScalarVariable
    variable = e_xy
    execute_on = timestep_end
  [../]
  [./eyz]
    type = ScalarVariable
    variable = e_yz
    execute_on = timestep_end
  [../]
  [./ezx]
    type = ScalarVariable
    variable = e_zx
    execute_on = timestep_end
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121            1e-8          1e-8       1e-6     bjacobi'
  [../]
[]


[Executioner]
  type = Transient
  start_time = 0
  dt = 0.1

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'
  num_steps = 50
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_full_LK
  [../]
[]
