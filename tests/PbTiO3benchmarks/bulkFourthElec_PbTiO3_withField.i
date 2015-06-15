# This input file tests the BulkEnergyDerivative Kernel for fourth order and time dependent kernel. The solution will be polarized domains in preferred directions
# for PbTiO3 at 20 K. We won't be able to use ExoDiff because the solution will be different everytime it is ran due to the nonlinearity of the problem
# Electrostatics is turned off.

[Mesh]
  file = cube.e
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
    order = FIRST
    family = LAGRANGE
  [../]
  [./potential_ext]
    order = FIRST
    family = LAGRANGE
  [../]
[]


[GlobalParams]
   len_scale=1e-9
   alpha1=-1.7252e8 # 3.8(T-479)*10^5 C^{-2}m^2
   alpha11=-7.3e7
   alpha111=2.6e8
   alpha12=7.5e8
   alpha112=6.1e8
   alpha123=-3.7e9
#   G110=4e-10
#   #G110=0.0
#   G11/G110=0.6
#   G12/G110=0.0
#   G44/G110=0.3
#   G44P/G110=0.3 #is this a problem?
   polar_x=polar_x
   polar_y=polar_y
   polar_z=polar_z
   permittivity=8.85e-12
   potential_ext = potential_ext
   potential_int = potential_int
[]

[Kernels]
  [./bed_x]
    type = BulkEnergyDerivative_nosixth
    variable = polar_x
    component=0
    implicit=false
  [../]
  [./bed_y]
    type = BulkEnergyDerivative_nosixth
    variable = polar_y
    component=1
    implicit=false
  [../]
  [./bed_z]
    type = BulkEnergyDerivative_nosixth
    variable = polar_z
    component=2
    implicit=false
  [../]
  [./polar_electric_E]
     type=PolarElectricEStrong
     variable=potential_int
     permittivity=8.85e-12
     implicit=false
  [../]
  [./FE_E_int]
     type=Electrostatics
     variable=potential_int
     permittivity=8.85e-12
  [../]
  [./E_int]
     type=Electrostatics
     variable=potential_int
     permittivity=8.85e-12
  [../]
  [./E_ext]
     type=Electrostatics
     variable=potential_ext
     permittivity=8.85e-12
  [../]
  [./FE_E_ext]
     type=Electrostatics
     #type=Diffusion
     variable=potential_ext
     permittivity=8.85e-12
  [../]
  [./polar_electric_px]
     type=PolarElectricPStrong
     variable=polar_x
     component=0
     implicit=false
  [../]
  [./polar_electric_py]
     type=PolarElectricPStrong
     variable=polar_y
     component=1
     implicit=false
  [../]
  [./polar_electric_pz]
     type=PolarElectricPStrong
     variable=polar_z
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

[ICs]
  active='polar_x_constic polar_y_constic polar_z_constic'
  [./polar_x_constic]
     type=ConstantIC
     variable=polar_x
     value = 0.6
  [../]
  [./polar_y_constic]
     type=ConstantIC
     variable=polar_y
     value = 0.6
  [../]
  [./polar_z_constic]
     type=ConstantIC
     variable=polar_z
     value = 0.6
  [../]
[]

[BCs]

  [./potential_int_upz]
    type = DirichletBC
    variable = potential_int
    boundary = '1'
    value = 0.0
  [../]

  [./potential_int_downz]
    type = DirichletBC
    variable = potential_int
    boundary = '2'
    value = 0.0
  [../]


   [./potential_ext_upz]
    type = NeumannBC
    variable = potential_ext
    boundary = '1'
  #  value = 1.0e-13 #1e-12 gives a field of ~1e7 V/m ~ 100 kV/cm for size 1?, what about size 100?
    value = 0.0
  [../]
  [./potential_ext_downz]
    type = NeumannBC
    variable = potential_ext
    boundary = '2'
  #    value = -1.0e-13
    value = 0.0
  [../]
[]


[Executioner]
  type=Transient
  solve_type=newton
  scheme=explicit-euler
  dt=1e17
  num_steps = 500
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname='-ksp_type -pc_type -snes_linesearch_type -pc_factor_zeropivot'
  petsc_options_value='  gmres     jacobi       basic                1e-50'
[]

[Postprocessors]
  [./bulk_energy]
    type=BulkEnergy
  [../]
[]


[Outputs]
  file_base = out_PbTiO3_cube_bulk_len-9_icConst_dt1e17_nosixth
  output_initial = true
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    elemental_as_nodal = true
    interval = 1
  [../]
[]
