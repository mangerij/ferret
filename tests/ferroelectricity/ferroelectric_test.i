[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 2
  ny = 2
  nz = 2
  ymax = 2
  ymin = -2
  xmin = -2
  xmax = 2
  zmin = -2
  zmax = 2
[]

[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
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

[GlobalParams]
   len_scale=1e-9
   alpha1=-1.8202e8 # 3.8(T-479)*10^5 C^{-2}m^2 (T=0 K)
   alpha11=-7.3e7
   alpha111=2.6e8
   alpha12=7.5e8
   alpha112=6.1e8
   alpha123=-3.7e9
   G110=0.6e-10
   G11/G110=0.6
   G12/G110=0.0
   G44/G110=0.3
   G44P/G110=0.3
   permittivity=8.85e-12
   polar_x=polar_x
   polar_y=polar_y
   polar_z=polar_z
   potential_int=potential_int
   potential_ext=potential_ext
   time_scale = 1.0e-29
[]

[Kernels]
  [./bed_x]
    type = BulkEnergyDerivative
    variable = polar_x
    component=0
  [../]
  [./bed_y]
    type = BulkEnergyDerivative
    variable = polar_y
    component=1
  [../]
  [./bed_z]
    type = BulkEnergyDerivative
    variable = polar_z
    component = 2
  [../]
  [./walled_x]
     type=WallEnergyDerivative
     variable=polar_x
     component = 0
  [../]
  [./walled_y]
     type=WallEnergyDerivative
     variable=polar_y
     component = 1
  [../]
  [./walled_z]
     type=WallEnergyDerivative
     variable=polar_z
     component = 2
  [../]
  [./polar_electric_E]
     type=PolarElectricEStrong
     variable=potential_int
  [../]
  [./FE_E_int]
     type=Electrostatics
     variable=potential_int
  [../]
  [./FE_E_ext]
     type=Electrostatics
     variable=potential_ext
  [../]
  [./polar_electric_px]
     type=PolarElectricPStrong
     variable=polar_x
     component = 0
  [../]
  [./polar_electric_py]
     type=PolarElectricPStrong
     variable = polar_y
     component = 1
  [../]
  [./polar_electric_pz]
     type=PolarElectricPStrong
     variable=polar_z
     component = 2
  [../]
  [./polar_x_time]
     type=TimeDerivativeScaled
     variable = polar_x
  [../]
  [./polar_y_time]
     type=TimeDerivativeScaled
     variable = polar_y
  [../]
  [./polar_z_time]
     type=TimeDerivativeScaled
     variable = polar_z
  [../]
[]

[ICs]
  [./polar_x_constic_rand]
     type = RandomIC
     variable = polar_x
     min = 0.6
     max = 0.75
  [../]
  [./polar_y_constic_rand]
     type = RandomIC
     variable = polar_y
     min = 0.6
     max = 0.75
  [../]
  [./polar_z_constic_rand]
     type = RandomIC
     variable = polar_z
     min = 0.0
     max = 0.0
  [../]
[]

[BCs]
  [./potential_int_upz]
    type = DirichletBC
    variable = potential_int
    boundary = 'top'
    value = 0.0
  [../]
  [./potential_int_downz]
    type = DirichletBC
    variable = potential_int
    boundary = 'bottom'
    value = 0.0
  [../]
  [./potential_ext_upz]
    type = DirichletBC
    variable = potential_ext
    boundary = 'top'
    value = 0.0
  [../]
  [./potential_ext_downz]
    type = DirichletBC
    variable = potential_ext
    boundary = 'bottom'
    value = 0.0
  [../]
[]

#[Postprocessors]
#   [./_pps_percent]
#     type = PercentChangePostprocessor
#     postprocessor = total_energy
#   [../]
#[]

#[UserObjects]
#  [./kill]
#    type = Terminator
#    expression = '_pps_percent <= 4.0e-7'
#  [../]
#[]

[Executioner]
  type=Transient

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.1e-15
    optimal_iterations = 3
    growth_factor = 1.001
    cutback_factor =  0.999
  [../]
  scheme = 'explicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dtmin=1.0e-29
  dtmax=1.70e-15
  num_steps = 350
  petsc_options='-snes_converged_reason'
  petsc_options_iname='-ksp_type -snes_type  -snes_rtol -ksp_rtol -pc_type -pc_hypre_boomeramg_strong_threshold -snes_linesearch_type -pc_factor_zeropivot'
  petsc_options_value=' gmres     newtonls       1e-8     1e-12      hypre     0.5                                 basic         1e-50  '
[]


[Outputs]
  print_linear_residuals = true
  print_perf_log = true

  [./out]
    type = Exodus
    file_base = out_ferroelectric_test
    output_initial = true
    elemental_as_nodal = false
    interval = 50
  [../]

[]
