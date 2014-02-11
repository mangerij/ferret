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
  [./u]
    order=FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./diffusion]
     type=Diffusion
     variable=u
  [../]
[]


[BCs]
  [./u_bottomy]
    type = DirichletBC
    variable = u
    boundary ='bottom'
    value = 0.0
  [../]
  [./u_topy]
    type = DirichletBC
    variable = u
    boundary ='top'
    value = 0.0
  [../]
  [./Periodic]
    [./u_perx]
       variable = u
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
  petsc_options='-snes_test_display'
  petsc_options_iname='-snes_type'
  petsc_options_value='      test'
[]

[Output]
  output_initial=1
  exodus = true
  perf_log = true
[]