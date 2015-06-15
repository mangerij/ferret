[Mesh]
  file = cube.e
  uniform_refine=1
[]

[Variables]
  [./polar_x]
    #scaling=1e-14
    order = FIRST
    family = LAGRANGE
    block='1'
  [../]
  [./polar_y]
    #scaling=1e-14
    order = FIRST
    family = LAGRANGE
    block='1'
  [../]
  [./polar_z]
    #scaling=1e-14
    order = FIRST
    family = LAGRANGE
    block='1'
  [../]
[]


[GlobalParams]
   len_scale=2e-8
   #len_scale=1.0
   alpha1=-1.7252e8 # 3.8(T-479)*10^5 C^{-2}m^2 for T=0
   alpha11=-7.3e7
   alpha111=2.6e8
   alpha12=7.5e8
   alpha112=6.1e8
   alpha123=-3.7e9
   G110=4e-10  #not sure what the 10^{4} term was, since Chen's paper points to this value as being correct
   #G110=0.0
   G11/G110=0.6
   G12/G110=0.0
   G44/G110=0.3
   G44P/G110=0.3
   permittivity=8.85e-12
   #permittivity=1.0
   polar_x=polar_x
   polar_y=polar_y
   polar_z=polar_z
   potential_int=potential_int
   potential_ext=potential_ext
[]


[Kernels]
#  active='diffusion_E diffusion_E_Ext  polar_electric_E polar_electric_px polar_electric_py polar_electric_pz polar_x_time polar_y_time polar_z_time'
  [./bed_x]
    type = BulkEnergyDerivative_scaled
    variable = polar_x
    component=0
    implicit=false
  [../]
  [./bed_y]
    type = BulkEnergyDerivative_scaled
    variable = polar_y
    component=1
    implicit=false
  [../]
  [./bed_z]
    type = BulkEnergyDerivative_scaled
    variable = polar_z
    component=2
    implicit=false
  [../]
  [./walled_x]
     type=WallEnergyDerivative_scaled
     variable=polar_x
     component=0
  [../]
  [./walled_y]
     type=WallEnergyDerivative_scaled
     variable=polar_y
     component=1
  [../]
  [./walled_z]
     type=WallEnergyDerivative_scaled
     variable=polar_z
     component=2
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
     type=RandomIC
     variable=polar_x
     min=-4.0
     max=4.0
  [../]
  [./polar_y_constic]
     type=RandomIC
     variable=polar_y
     min=-4.0
     max=4.0
  [../]
  [./polar_z_constic]
     type=RandomIC
     variable=polar_z
     min=-4.0
     max=4.0
  [../]
[]


[Executioner]
  type=Transient            #Transient, Steady
  solve_type=newton
  scheme=explicit-euler     #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dt=1e9
  nl_max_its=30
  num_steps=1500
  #petsc_options="-snes_monitor -snes_converged_reason -ksp_monitor -ksp_converged_reason"
 # petsc_options='-snes_monitor -snes_converged_reason -ksp_monitor -ksp_converged_reason'
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname='-gmres_restart -snes_rtol -ksp_type  -ksp_rtol -pc_type -snes_linesearch_type -pc_factor_zeropivot'
  petsc_options_value='10          1e-8      gmres       1e-10    jacobi      basic                1e-50'
  #petsc_options_iname='-snes_rtol'
  #petsc_options_value='1e-16'
[]
[Postprocessors]
   [./bulk_energy]
    type=BulkEnergy
   [../]
   [./wall_energy]
    type=WallEnergy
   [../]
[]

[Outputs]
  file_base = out_PbTiO3_cube_len(2e-8)_dimL10W10H1_dt1e9_n1500_E0const_nowall_coarse_scaled
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
