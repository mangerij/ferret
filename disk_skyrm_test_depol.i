
[Mesh]
  file = exodus_disk_r8_2d.e
[]

[Problem]
  type = FEProblem
  coord_type = RZ
  rz_coord_axis = X
[]

[MeshModifiers]
  [./centernodeset_1]
    type = AddExtraNodeset
    new_boundary = 'center_node_1'
    coord = '0.0 0.0'
  [../]
[]

[GlobalParams]
  len_scale = 1.0
  P = P
  theta = theta
[]

[Variables]
  [./P]
    order = THIRD
    family = HERMITE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = -0.5e-3
      max = 0.5e-3
    [../]
  [../]
  [./theta]
    order = THIRD
    family = HERMITE
    block = '1'
    [./InitialCondition]
      type = RandomIC
      min = 0.0
      max = 5.0
    [../]
  [../]
[]

[Kernels]
    [./theta_lap]
       type = EulerSkyrmionThetaTerm
       variable = theta
       xi0 = 1.0
    [../]
    [./P_lap]
       type = EulerSkyrmionPTerm
       variable = P
       xi0 = 1.0
    [../]
    [./thetaK]
      type = EulerSkyrmionThetaKappaTerm
      variable = theta
      xi0 = 1.0
      kappa = 0.4
    [../]
    [./thetadepol]
      type = EulerSkyrmionThetaDepolTerm
      variable = theta
      edep = 0.02
    [../]
    [./Pdepol]
      type = EulerSkyrmionPDepolTerm
      variable = P
      edep = 0.02
    [../]
    [./Pc_term]
      type = EulerSkyrmionPCubeTerm
      variable = P
      P0 = 1.0
    [../]
    [./Ptemp_term]
      type = EulerSkyrmionPTempTerm
      variable = P
      kappa = 0.4
      t = -0.65
      xi0 = 1.0
    [../]
[]


[BCs]
   [./middle_P_condition]
     type = PresetBC
     value = 0.75
     boundary = 'center_node_1'
     variable = 'P'
   [../]
   [./middle_th_condition]
     type = PresetBC
     value = 0.0
     boundary = 'center_node_1'
     variable = 'theta'
   [../]
   [./radial_P_condition]
     type = NeumannBC
     value = 0.0
     boundary = '1'
     variable = 'P'
   [../]
   [./radial_th_condition]
     type = NeumannBC
     value = 0.0
     boundary = '1'
     variable = 'theta'
   [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_view -snes_linesearch_monitor -snes_converged_reason -ksp_converged_reason'
    petsc_options_iname = '-ksp_gmres_restart  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '        121                1e-8      1e-14    lu     '
  [../]
[]

[Executioner]

  type = Steady
  solve_type = 'NEWTON'       #"PJFNK, JFNK, NEWTON"
  #line_search = none
[]

[Outputs]
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    file_base = out_nanodisk_scyl_skyrm
    elemental_as_nodal = true
    interval = 1
  [../]
[]
