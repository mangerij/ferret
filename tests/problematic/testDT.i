[Mesh]
   file = poissonstripe_coarse.e
   uniform_refine=0
[]
# [Mesh]
#   type=GeneratedMesh
#   dim=2
#   nx=5
#   ny=5
#   xmin=0.0
#   xmax=1.0
#   ymin=0.0
#   ymax=1.0
# []
[Variables]
  [./potential_int]
    order=FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./diffusion_E]
     type=Diffusion
     variable=potential_int
  [../]
 [./diffusion_E_time]
    type=TimeDerivative
    variable=potential_int
 [../]
[]

[ICs]
  [./potential_int_ic]
     type=FunctionIC
     variable=potential_int
     function=initial_cond_func
  [../]
[]

[Functions]
  [./initial_cond_func]
     type=ParsedFunction
     value=z*z
  [../]
[]

[BCs]
  [./potential_int_upz]
    type = DirichletBC
    variable = potential_int
    boundary = 'upz'
    value = 1.0
  [../]
  [./potential_int_downz]
    type = DirichletBC
    variable = potential_int
    boundary = 'downz'
    value = 0.0
  [../]
[]
[Preconditioning]
   [./smp]
     type=SMP   #or SMP
     full=true   #to use every off diagonal block
     pc_side=left
   [../]
[]

[Executioner]
  #type = Steady
  type=Transient
  solve_type=newton
  scheme=explicit-euler     #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  nl_max_its=100
  num_steps=3
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname='-gmres_restart -snes_rtol -ksp_type  -ksp_rtol -pc_type -snes_linesearch_type'
  petsc_options_value='1000            1e-5        gmres       1e-8      jacobi       basic'
  [./TimeStepper]
     type=CustomDT
     #type=ConstantDT
     dt=1
     increase_rate=1.00001
     #increase_rate=100.0
     postprocessor=total_energy
  [../]
[]
[Postprocessors]
    [./total_energy]
    type=ElementAverageValue
    variable=potential_int
    [../]
[]

[Outputs]
  output_initial = true
  exodus = true
  [./console]
    type = Console
    perf_log = true
  [../]
[]