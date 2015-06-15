[Mesh]
  file = cube.e
  uniform_refine=0
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

[]

[AuxVariables]
#[./auxv_es_energy_density_e] #es for electrostatic
#     order=CONSTANT
#     family=MONOMIAL
#  [../]
#  [./auxv_es_energy_density]
#     order=CONSTANT
#     family=MONOMIAL
#  [../]
#  [./auxv_bulk_energy_density]
#     order=CONSTANT
#     family=MONOMIAL
#  [../]
[]

[GlobalParams]
   len_scale=1e-9
   #len_scale=1.0
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
   G44P/G110=0.3
#   permittivity=8.85e-12
   #permittivity=1.0
   polar_x=polar_x
   polar_y=polar_y
   polar_z=polar_z
#   potential_int=potential_int
#   potential_ext=potential_ext
[]

[AuxKernels]
#  #active='diff'
#  [./auxk_electrostatic_energy_density_e]
#    type =ElectrostaticEnergyDensityE
#    variable =auxv_es_energy_density_e
#    potential=potential_int
#  [../]
#  [./auxk_electrostatic_energy_density]
#    type =ElectrostaticEnergyDensity
#    variable =auxv_es_energy_density
#  [../]
#  [./auxk_bulk_energy_density]
#    type =BulkEnergyDensity
#    variable =auxv_bulk_energy_density
#  [../]
[]

[Kernels]
  [./bed_x]
    type = BulkEnergyDerivative
    variable = polar_x
    component=0
    implicit=false
  [../]
  [./bed_y]
    type = BulkEnergyDerivative
    variable = polar_y
    component=1
    implicit=false
  [../]
  [./bed_z]
    type = BulkEnergyDerivative
    variable = polar_z
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
     value = 0.1
  [../]
  [./polar_y_constic]
     type=ConstantIC
     variable=polar_y
     value = -1.0
  [../]
  [./polar_z_constic]
     type=ConstantIC
     variable=polar_z
     value = 1.0
  [../]
[]

[BCs]

#  [./potential_int_upz]
#    type = DirichletBC
#    variable = potential_int
#    boundary = '1'
#    value = 0.0
#    #implicit=false
#  [../]
#  [./potential_int_downz]
#    type = DirichletBC
#    variable = potential_int
#    boundary = '2'
#    value = 0.0
#    #implicit=false
#  [../]
#   [./potential_ext_upz]
#    type = DirichletBC
#    variable = potential_ext
#    boundary = '1'
#    value = 0.0
#    #implicit=false
#  [../]
#  [./potential_ext_downz]
#    type = DirichletBC
#    variable = potential_ext
#    boundary = '2'
#    value = 0.0
#    #implicit=false
#  [../]
[]

[Executioner]
  #type = Steady
  type=Transient
  solve_type=newton
  scheme=explicit-euler     #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dt=1e17
  nl_max_its=100
  num_steps=3500
  #petsc_options="-snes_monitor -snes_converged_reason -ksp_monitor -ksp_converged_reason"
 # petsc_options='-snes_monitor -snes_converged_reason -ksp_monitor -ksp_converged_reason'
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname='-gmres_restart -snes_rtol -ksp_type  -ksp_rtol -pc_type -snes_linesearch_type -pc_factor_zeropivot'
  petsc_options_value='1000            1e-6        gmres       1e-10      jacobi       basic                1e-50'
  #petsc_options_iname='-snes_rtol'
  #petsc_options_value='1e-16'
[]
[Postprocessors]
  [./bulk_energy]
    type=BulkEnergy
   [../]
#   [./wall_energy]
#    type=WallEnergy
#   [../]
#   [./electric_energy]
#    type=ElectrostaticEnergy
#    [../]
#    [./total_energy]
#    type=TotalEnergy
#    bulk_energy=bulk_energy
#    wall_energy=wall_energy
#    electric_energy=electric_energy
#    [../]
[]


[Outputs]
  file_base = out_PbTiO3_cube_bulk_len-9_icConst_dt1e17_n3500_run1
  output_initial = true
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    elemental_as_nodal = true
    interval = 1
  [../]
  [./debug]
    type = VariableResidualNormsDebugOutput
  [../]
[]
