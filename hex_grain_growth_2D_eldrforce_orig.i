[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 42
  nz = 0
  xmax = 10.0
  ymax = 8.66
  zmax = 0
  elem_type = QUAD4
  uniform_refine = 2
[]

[GlobalParams]
  op_num = 7
  var_name_base = gr
[]

[Variables]
  [./PolycrystalVariables]
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    coloring_algorithm = bt
    grain_num = 12
    rand_seed = 10
    output_adjacency_matrix = true
  [../]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = grn_36_test2_2D.tex
  [../]
  [./grain_tracker]
    type = GrainTrackerElasticity
    threshold = 0.2
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_begin'
    flood_entity_type = ELEMENTAL
    fill_method = symmetric9
    C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5'
    euler_angle_provider = euler_angle_file
    reserve_op = 1
    reserve_op_threshold = 0.95
  [../]


  [./inserter]
    type = DiscreteNucleationInserter
    hold_time = 5e-10
    probability = 1
    seed = 12346
    radius = 3.27
  [../]
  [./map]
    type = DiscreteNucleationMap
    int_width = 1
    periodic = gr6
    inserter = inserter
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./elastic_strain11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./C1111]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler_angle]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./nuc_map]
  order = CONSTANT
  family = MONOMIAL
[../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
  [./PolycrystalElasticDrivingForce]
  [../]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
  [../]
  [./c_force]
    type = DiscreteNucleationForce
    variable = gr6
    map = map
    no_nucleus_value = 0
    nucleus_value = 1
  [../]
  [./c_react]
    type = Reaction
    variable = gr6
  [../]
[]

[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  [../]
  [./elastic_strain11]
    type = RankTwoAux
    variable = elastic_strain11
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  [../]
  [./elastic_strain22]
    type = RankTwoAux
    variable = elastic_strain22
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./elastic_strain12]
    type = RankTwoAux
    variable = elastic_strain12
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    execute_on = timestep_end
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  [../]
  [./C1111]
    type = RankFourAux
    variable = C1111
    rank_four_tensor = elasticity_tensor
    index_l = 0
    index_j = 0
    index_k = 0
    index_i = 0
    execute_on = timestep_end
  [../]
  [./vonmises_stress]
    type = RankTwoScalarAux
    variable = vonmises_stress
    rank_two_tensor = stress
    scalar_type = VonMisesStress
  [../]
  [./euler_angle]
    type = OutputEulerAngles
    variable = euler_angle
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
  [../]

  [./nuc_mapo]
    type = DiscreteNucleationAux
    variable = nuc_map
    map = map
  [../]

[]

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x'
      variable = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 disp_x disp_y'
    [../]
  [../]
  [./top_displacement]
    type = DirichletBC
    variable = disp_y
    boundary = top
    value = 0.5 #50
  [../]
  [./x_anchor]
    type = DirichletBC
    variable = disp_x
    boundary = 'left right'
    value = 0.0
  [../]
  [./y_anchor]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
[]

[Materials]
  [./Copper]
    type = GBEvolution
    block = 0
    T = 900 # K
    wGB = 0.15 #15 # nm
    GBmob0 = 2.5e-6 # m^4/(Js) from Schoenfelder 1997
    Q = 0.23 # Migration energy in eV
    GBenergy = 0.708 # GB energy in J/m^2
    length_scale = 1e-7
    time_scale = 1
  [../]
  [./ElasticityTensor]
    type = ComputePolycrystalElasticityTensor
    block = 0
    grain_tracker = grain_tracker
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = 0
    displacements = 'disp_x disp_y'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = 0
  [../]
  [./probability]
    # This is a made up toy nucleation rate it should be replaced by
    # classical nucleation theory in a real simulation.
    type = ParsedMaterial
    property_name = P
    coupled_variables = bnds
    expression = 'if(bnds<0.21,bnds*1e-8,0)'
    outputs = exodus
  [../]
[]

[Postprocessors]
  [./dofs]
    type = NumDOFs
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./bnd_length]
    type = GrainBoundaryArea
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    off_diag_row = 'disp_x disp_y'
    off_diag_column = 'disp_y disp_x'
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'
  l_tol = 1.0e-4
  l_max_its = 60
  nl_max_its = 15
  nl_rel_tol = 1.0e-7
  start_time = 0.0
  num_steps = 1000
  automatic_scaling = true
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1E-10
    growth_factor = 1.05
    cutback_factor = 0.85
    optimal_iterations = 7
  [../]
  [./Adaptivity]
    initial_adaptivity = 2
    refine_fraction = 0.5
    coarsen_fraction = 0.2
    max_h_level = 2
  [../]
[]

[Outputs]
  file_base = 1Nano_Vox_Half_Mesh
  exodus = true
  print_linear_residuals = false
  interval = 1
[]
