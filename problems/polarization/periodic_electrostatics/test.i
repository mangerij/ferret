[Mesh]
  type=GeneratedMesh
  dim=2
  nx=2
  ny=1
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
[]

[Variables]
active='potential_int'
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
[]


[BCs]
  [./potential_int_bottomy]
    type = DirichletBC
    variable = potential_int
    boundary ='bottom'
    value = 0.0
  [../]
  [./potential_int_topy]
    type = DirichletBC
    variable = potential_int
    boundary ='top'
    value = 0.0
  [../]
  [./Periodic]
    active='potential_int_x'
    [./potential_int_x]
       variable = potential_int
       primary ='left'
       secondary ='right'
       translation = '1 0 0'
    [../]
  [../]
[]
[Preconditioning]
   [./smp]
     type=SMP
     pc_side=left
   [../]
[]

[Executioner]
  type = Steady
  solve_type=newton
  # num_steps=400
  petsc_options='-snes_monitor -snes_view -snes_converged_reason  -ksp_monitor_true_residual -snes_linesearch_monitor -options_left -snes_test_display'
  petsc_options_iname='-snes_linesearch_type -snes_type'
  petsc_options_value='basic                       test'
[]

[Output]
  output_initial=1
  exodus = true
  perf_log = true
[]