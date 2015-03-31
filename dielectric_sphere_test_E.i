
[Mesh]
  file = sphere_150_25_coarse.e
  #uniform_refine=1
[]
[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block='2'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block='2'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block='2'
  [../]
  [./potential_int]
    order=FIRST
    family = LAGRANGE
  [../]
  [./potential_ext]
    order=FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./SurfCharge]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Ex]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ey]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ez]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./depol_Ex]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./depol_Ey]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./depol_Ez]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
 [./polar_electric_E]
     type=PolarElectricEStrong
     variable=potential_int
     block='2'
     permittivity = 1  # permittivity = 8.85*e-12
     polar_x = polar_x
     polar_y = polar_y
     polar_z = polar_z
     #implicit=false
  [../]

  [./diffusion_E_int_block2]
    type=Electrostatics
     permittivity = 4 #permittivity = 4*8.85*e-12
     variable=potential_int
     block='2'
  [../]
  [./diffusion_E_int_block1]
    type=Electrostatics
     permittivity = 1 #permittivity = 4*8.85*e-12
     variable=potential_int
     block='1'
  [../]


  [./diffusion_E_Ext_block2]
     type=Electrostatics
     permittivity = 4 #permittivity = 8.85*e-12
     variable=potential_ext
     block='2'
  [../]

  [./diffusion_E_Ext_block1]
     type=Electrostatics
     permittivity = 1 #permittivity = 8.85*e-12
     variable=potential_ext
     block='1'
  [../]
  

  [./polar_electric_px]
     type=PolarElectricPStrong
     variable=polar_x
     permittivity_ext = 1
     permittivity_int = 4
     block = '2'
     potential_ext = potential_ext
     potential_int = potential_int
     component=0
     #implicit=false
  [../]

  [./polar_electric_py]
     type=PolarElectricPStrong
     variable=polar_y
     permittivity_ext = 1
     permittivity_int = 4
     block = '2'
     potential_ext = potential_ext
     potential_int = potential_int
     component=1
     #implicit=false
  [../]
  [./polar_electric_pz]
     type=PolarElectricPStrong
     variable=polar_z
     block = '2'
     permittivity_ext = 1
     permittivity_int = 4
     potential_ext = potential_ext
     potential_int = potential_int
     component=2
     #implicit=false
  [../]
  [./polar_x_time]
     type=TimeDerivative
     variable=polar_x
  [../]
  [./polar_y_time]
     type=TimeDerivative
     variable=polar_y
  [../]
  [./polar_z_time]
     type=TimeDerivative
     variable=polar_z
  [../]
[]

[AuxKernels]
   [./surfacechargeaux]
    type = SurfaceChargeAux
    variable = SurfCharge
    boundary = '7'
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
  [../]
  [./Ex_fieldAux]
    type = Ex_fieldAux
    variable = Ex
    potential_int = potential_int
    potential_ext = potential_ext
  [../]
  [./Ey_fieldAux]
    type = Ey_fieldAux
    variable = Ey
    potential_int = potential_int
    potential_ext = potential_ext
  [../]
  [./Ez_fieldAux]
    type = Ez_fieldAux
    variable = Ez
    potential_int = potential_int
    potential_ext = potential_ext
  [../]

  [./Depol_x_fieldAux]
    type = Depol_x_fieldAux
    variable = depol_Ex
    block = '2'
    potential_int = potential_int
  [../]
  [./Depol_y_fieldAux]
    type = Depol_y_fieldAux
    variable = depol_Ey
    block = '2'
    potential_int = potential_int
  [../]
  [./Depol_z_fieldAux]
    type = Depol_z_fieldAux
    variable = depol_Ez
    block = '2'
    potential_int = potential_int
  [../]
[]

[BCs]
  [./potential_ext_1]
    type = NeumannBC  
    variable = potential_ext
    boundary = '1'
    value = 1.0
  [../]
  [./potential_ext_2]
    type = NeumannBC
    variable = potential_ext
    boundary = '2'
    value = -1.0
  [../]
  [./potential_ext_3]
    type = NeumannBC
    variable = potential_ext
    boundary = '3'
    value = 0.0
  [../]
  [./potential_ext_4]
    type = NeumannBC
    variable = potential_ext
    boundary = '4'
    value = 0.0
  [../]
   [./potential_ext_5]
    type = NeumannBC
    variable = potential_ext
    boundary = '5'
    value = 0.0
  [../]
  [./potential_ext_6]
    type = NeumannBC
    variable = potential_ext
    boundary = '6'
    value = 0.0
  [../]

  [./potential_int_1]
    type = DirichletBC
    variable = potential_int
    boundary = '1'
    value = 0
  [../]
  [./potential_int_2]
    type = DirichletBC
    variable = potential_int
    boundary = '2'
    value = 0
  [../]
  [./potential_int_3]
    type = DirichletBC
    variable = potential_int
    boundary = '3'
    value = 0
  [../]
  [./potential_int_4]
    type = DirichletBC
    variable = potential_int
    boundary = '4'
    value = 0
  [../]
  [./potential_int_5]
    type = DirichletBC
    variable = potential_int
    boundary = '5'
    value = 0
  [../]
  [./potential_int_6]
    type = DirichletBC
    variable = potential_int
    boundary = '6'
    value = 0
  [../]
[]

[ICs]
  active='polar_x_constic polar_y_constic polar_z_constic'
  [./polar_x_constic]
     type=ConstantIC
     variable=polar_x
     block = '2'
     value = -0.3
  [../]
  [./polar_y_constic]
     type=ConstantIC
     variable=polar_y
     block = '2'
     value = -0.3
  [../]
  [./polar_z_constic]
     type=ConstantIC
     variable=polar_z
     block = '2'
     value = -0.3
  [../]
[]

[Preconditioning]
   [./smp]
     type=SMP
     full=true   #to use every off diagonal block
     pc_side=left
   [../]
[]

[Executioner]
  #type = Steady
  type=Transient
  solve_type=newton
  scheme=implicit-euler
  #l_max_its=1500
  dt=1e0
  num_steps=15 #converges about 15 iterations for strong coupling
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname='-gmres_restart -ksp_type  -pc_type -snes_linesearch_type -pc_factor_zeropivot'
  petsc_options_value='1000             gmres     jacobi       basic                1e-50'
[]

#[Debug]
#  show_parser = true
#[]

[Outputs]
  file_base = outlin_die_sph_strong_implic_dt0_n15_er4_E0-1_size25_150
  output_initial = true
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    elemental_as_nodal = true
  [../]
[]
