[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 10
    nz = 20
    xmin = -10.0
    xmax = 10.0
    ymin = -10.0
    ymax = 10.0
    zmin = -10.0
    zmax = 10.0
    elem_type = HEX8
  []
  [./cnode]
    input = gen
    type = ExtraNodesetGenerator
    coord = '-10.0 -10.0 -10.0'
    new_boundary = 100
  [../]

  [subdomains]
    type = SubdomainBoundingBoxGenerator
    input = cnode
    bottom_left = '-10.0 -10.0 -10.0'
    block_id = 1
    top_right = '10.0 10.0 0.0'
    location = INSIDE
  []
  [film_interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = subdomains
    primary_block = 0
    paired_block = 1
    new_boundary = 52
  []
[]

[GlobalParams]
  displacements = 'u_x u_y u_z'
[]

[Variables]
  [./global_strain_film]
    order = SIXTH
    family = SCALAR
    block = '0'
  [../]

  [./global_strain_sub]
    order = SIXTH
    family = SCALAR
    block = '1'
  [../]

  [./u_x]
    order = FIRST
    family = LAGRANGE
    block = '0 1'
  [../]
  [./u_y]
    order = FIRST
    family = LAGRANGE
    block = '1 0'
  [../]
  [./u_z]
    order = FIRST
    family = LAGRANGE
    block = '1 0'
  [../]
[]


[Materials]

  [./eigen_strain]
    type = ComputeEigenstrain
    #    eigen_base = 'x y z yz xz xy'
    eigen_base = '0.005 0.005 0.003725 0 0 0'
    eigenstrain_name = epitaxy
    block = '0'
  [../]

  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '179.073 66.71 66.71 179.073 66.71 179.073 82.6446 82.6446 82.6446'
  [../]

  [./stress_1]
    type = ComputeLinearElasticStress
  [../]
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [./Ferroelectric_SM]
        strain = SMALL
        add_variables = true
        incremental = false
        global_strain = global_strain
        eigenstrain_names = 'epitaxy'
        generate_output = 'strain_xx strain_xy strain_xz strain_yy strain_yz strain_zz vonmises_stress'
	       block = '0'
      [../]

      [./Substrate_SM]
        strain = SMALL
        add_variables = true
        global_strain = global_strain
  #      eigenstrain_names = 'epitaxy'
        incremental = false
        generate_output = 'strain_xx strain_xy strain_xz strain_yy strain_yz strain_zz vonmises_stress'
        block = '1'
      [../]
    []

    [./GlobalStrain]
      [./global_strain_film]
        scalar_global_strain = global_strain_film
        displacements = 'u_x u_y u_z'
        auxiliary_displacements = 'disp_x disp_y disp_z'
        global_displacements = 'ug_x ug_y ug_z'
        block = '0'
      [../]
      [./global_strain_sub]
        scalar_global_strain = global_strain_sub
        displacements = 'u_x u_y u_z'
        auxiliary_displacements = 'disp_x disp_y disp_z'
        global_displacements = 'ug_x ug_y ug_z'
        block = '1'
      [../]
    [../]
  []
[]

[BCs]
  [./Periodic]
    [./xy]
      auto_direction = 'x y'
      variable = 'u_x u_y u_z'
    [../]
  [../]
  [./centerfix_x]
    type = DirichletBC
    boundary = 100
    variable = u_x
    value = 0
  [../]
  [./centerfix_y]
    type = DirichletBC
    boundary = 100
    variable = u_y
    value = 0
  [../]
  [./centerfix_z]
    type = DirichletBC
    boundary = 100
    variable = u_z
    value = 0
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    160             1e-8        1e-6      1e-5        bjacobi'
  [../]
[]

[Debug]
  show_material_props = true
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  scheme = 'implicit-euler'
  dtmin = 1e-13
  dtmax = 2.0
  end_time = 10
  l_max_its = 200
  automatic_scaling = true
  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 8
    growth_factor = 1.1
    cutback_factor = 0.8
    linear_iteration_ratio = 1000
    dt = 0.5
  [../]
  verbose = true
  nl_max_its = 16
[]

[Outputs]
  print_linear_residuals = false
  perf_graph = false
  [./out]
    type = Exodus
    file_base = out_GS_Test
    interval = 2
  [../]
[]
