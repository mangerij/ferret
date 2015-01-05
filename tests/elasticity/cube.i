[Mesh]
  #file = mug.e
  type=GeneratedMesh
  dim=3
  nx=5
  ny=5
  nz=5
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  zmin=0.0
  zmax=1.0
[]

[Variables]

  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]

  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]

  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]

[] # Variables

[TensorMechanics]
  [./solid]
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
  [../]
[]

[BCs]
  [./anchor_upz_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'front'
    value = 0
  [../]
  [./anchor_upz_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'front'
    value = 0
  [../]
  [./anchor_upz_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'front'
    value = 0.1
  [../]
  [./anchor_downz_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'back'
    value = 0
  [../]
  [./anchor_downz_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'back'
    value = 0
  [../]
  [./anchor_downz_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = -0.1
  [../]
[]

[Materials]
  [./whole]
    type = LinearElasticMaterial
    block = '0'
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z

    all_21 = false

    #Steel: Young's Modulus 200GPa, Poisson ratio 0.29
    C_ijkl='262 107 107 262 107 262 77 77 77'

    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]
[] # Materials

[Preconditioning]
   [./smp]
     type=SMP
     full=true
     pc_side=left
   [../]
[]

[Executioner]
  type = Steady
  solve_type=newton
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
[]

[Outputs]
  exodus = true
  output_on = 'initial timestep_end'
  [./console]
    type = Console
    perf_log = true
    output_on = 'timestep_end failed nonlinear'
  [../]
[]