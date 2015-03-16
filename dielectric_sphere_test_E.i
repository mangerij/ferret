
[Mesh]
  file = cube_sphere.e
 # uniform_refine=1
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
    famiy = LAGRANGE
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
[]

[Kernels]
  [./polar_electric_E]
     type=PolarElectricE
     variable=potential_int
     block='2'
     permittivity = 8.85*e-12
     polar_x = polar_x
     polar_y = polar_y
     polar_z = polar_z
     implicit=false
  [../]
  [./diffusion_E]
     type=Electrostatics
     permittivity = 8.85*e-12
     variable=potential_int
     block='1 2'
  [../]
  [./diffusion_E_Ext]
     type=Electrostatics
     #type=Diffusion
     permittivity = 8.85*e-12
     variable=potential_ext
     block='1 2'
  [../]
  [./polar_electric_px]
     type=PolarElectricP
     variable=polar_x
     potential_ext = potential_ext
     potential_int = potential_int
     component=0
     implicit=false
  [../]
  [./polar_electric_py]
     type=PolarElectricP
     variable=polar_y
     potential_ext = potential_ext
     potential_int = potential_int
     component=1
     implicit=false
  [../]
  [./polar_electric_pz]
     type=PolarElectricP
     variable=polar_z
     potential_ext = potential_ext
     potential_int = potential_int
     component=2
     implicit=false
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
  [./surfacechargeaux
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
[]

[BCs]
  [./potential_ext_1]
    type = NeumannBC
    variable = potential_ext
    boundary = '1'
    value = -1.0e5
  [../]
  [./potential_ext_2]
    type = NeumannBC
    variable = potential_ext
    boundary = '2'
    value = 1.0e5
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
[]

[ICs]
  active='polar_x_constic polar_y_constic polar_z_constic'
  [./polar_x_constic]
     type=ConstantIC
     variable=polar_x
     block = '2'
     value=1.0
  [../]
  [./polar_y_constic]
     type=ConstantIC
     variable=polar_y
     block = '2'
     value=1.0
  [../]
  [./polar_z_constic]
     type=ConstantIC
     variable=polar_z
     block = '2'
     value=1.0
  [../]
[]

[Preconditioning]
   [./smp]
     type=SMP   #or SMP
     #off_diag_row='var_name'
     #off_diag_column='var_name'
     full=true   #to use every off diagonal block
     pc_side=left
     #petsc_options='snes_mf_operator'
     #petsc_options_iname = '-pc_type -mat_fd_coloring_err -mat_fd_type'
     #petsc_options_value = 'lu       1e-6                 ds'
   [../]
[]

[Executioner]
  #type = Steady
  type=Transient
  solve_type=newton
  scheme=explicit-euler     #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dt=1e-1
  nl_max_its=500
  num_steps=200
  #petsc_options="-snes_monitor -snes_converged_reason -ksp_monitor -ksp_converged_reason"
 # petsc_options='-snes_monitor -snes_converged_reason -ksp_monitor -ksp_converged_reason'
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname='-gmres_restart -snes_rtol -ksp_type  -ksp_rtol -pc_type -snes_linesearch_type -pc_factor_zeropivot'
  petsc_options_value='1000            1e-8        gmres       1e-8      jacobi       basic                1e-50'
  #petsc_options_iname='-snes_rtol'
  #petsc_options_value='1e-16'
[]

[Debug]
  show_parser = true
[]

[Outputs]
  file_base = out_cube_die_sphere_explicit
  output_initial = true
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    elemental_as_nodal = true
    output_nodal = true
  [../]
[]
