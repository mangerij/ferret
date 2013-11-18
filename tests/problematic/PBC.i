[Mesh]
  type=GeneratedMesh
  dim=3
  nx=4
  ny=4
  nz=4
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  zmin=0.0
  zmax=1.0
[]

[Variables]
active='u'
  [./u]
    order=FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./diffusion_E]
     type=Diffusion
     variable=u
  [../]
[]


[BCs]
  [./Periodic]
    active='potential_int_x potential_int_y'
    [./potential_int_x]
       variable = u
       primary ='left'
       secondary ='right'
       translation = '1 0 0'
    [../]
    [./potential_int_y]
       variable = u
       primary = 'bottom'
       secondary ='top'
       translation = '0 1 0'
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
  petsc_options='-snes_monitor -snes_view -snes_converged_reason  -ksp_monitor_true_residual -snes_linesearch_monitor -options_left'
  petsc_options_iname='-snes_linesearch_type -snes_type'
  petsc_options_value='basic                test'
[]

[Output]
  output_initial=1
  exodus = true
  perf_log = true
[]