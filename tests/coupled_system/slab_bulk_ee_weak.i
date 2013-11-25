##This example will eventually failed with diverge line search. However, when it happens, we have reached the expected solution.
[Mesh]
  file=slab.e
[]
[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block='interior exterior'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block='interior exterior'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block='interior exterior'
  [../]
 []

[AuxVariables]
  [./auxv_bulk_energy_density]
     order=CONSTANT
     family=MONOMIAL
  [../]
[]

[GlobalParams]
   len_scale=1e-7
   energy_scale=1e12
   alpha1=-1.7252e8 # 3.8(T-479)*10^5 C^{-2}m^2
   alpha11=-7.3e7
   alpha111=2.6e8
   alpha12=7.5e8
   alpha112=6.1e8
   alpha123=-3.7e9
   permittivity=8.85e-12
   polar_x=polar_x
   polar_y=polar_y
   polar_z=polar_z
[]

[AuxKernels]
   [./auxk_bulk_energy_density]
    type =BulkEnergyDensity
    variable =auxv_bulk_energy_density
  [../]
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
  [./polar_x_constic]
     type=ConstantIC
     variable=polar_x
     value=0.6
  [../]
  [./polar_y_constic]
     type=ConstantIC
     variable=polar_y
     value=0.6
  [../]
  [./polar_z_constic]
     type=ConstantIC
     variable=polar_z
     value=0.6
  [../]
[]

[Preconditioning]
   [./smp]
     type=SMP   #or SMP
     full=true
     pc_side=left
   [../]
[]

[Executioner]
  #type = Steady
  type=Transient
  solve_type=newton
  scheme=explicit-euler     #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dt=0.1
  nl_max_its=100
  num_steps=800
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname='-gmres_restart -snes_rtol -ksp_type  -ksp_rtol -pc_type -snes_linesearch_type -pc_factor_zeropivot'
  petsc_options_value='1000            1e-5        gmres       1e-8      jacobi       basic                1e-50'
[]

[Postprocessors]
   [./bulk_energy]
      type=BulkEnergy
   [../]
[]

[Output]
  output_initial=1
  exodus = true
  perf_log = true
[]