###############################################################
##
##  Ferroelectric 180-degree domain wall in BaTiO3.
##
##  This input is the NON-FLEXOELECTRIC, UNROTATED reduction of the
##  flexoelectric domain-wall model.  Two systematic changes were made
##  relative to that model:
##
##   1. RotatedCubicParentNVFlexoelectricFullPDerivative ->
##      CubicParentElasticPDerivative.  The "Full" kernel carried BOTH the
##      electrostrictive and the flexoelectric parts of dF/dP (which is why
##      there was no separate elastic P-derivative Kernel).  Every gradient-
##      test term in it carries a flexocoupling factor mu, so as mu -> 0 it
##      reduces exactly to the elastic P-derivative.  Verified numerically:
##      running the flexo input with mu11 = mu12 = mu44 = 0 reproduces this
##      input's elastic energy to 4e-8 relative and total energy to 1.2e-8.
##
##   2. All Rotated* objects -> their unrotated counterparts, and the twist
##      angle / euler angles removed.  Nothing here is rotated, so the crystal
##      axes and the simulation axes coincide.
##
##  The flexoelectric eigenstrain, the flexocoupling constants mu_ij and the
##  elastic compliances s_ij are gone with them.  NOTE that the flexo
##  eigenstrain entered the CompositeEigenstrain with weight -1, so dropping
##  it is not a no-op.
##
##  GRADIENT COEFFICIENT CONVENTION -- read before editing the G block.
##  The unrotated WallEnergyDerivative/WallEnergy take G110 together with the
##  RATIOS G11/G110 etc.  The Rotated variants took bare G11/G12/G44.  Matching
##  the two residuals term-by-term at twist = 0 gives
##      G110 = G11,  G11_G110 = 1,  G12_G110 = G12/G11,  G44_G110 = G44/G11
##  and, crucially, G44P_G110 = 0: the Rotated kernel carries NO antisymmetric
##  G44' term.  Other unrotated Ferret inputs ship a nonzero G44P_G110
##  (surface_charge.i uses 0.3), so copying a G block from one of those would
##  silently add physics this model never had.
##
##  Units in Ferret are nm, kg, seconds and attocoulombs.
##
##  CAVEAT: the Terminator below stops on the ENERGY RATE, which is blind to
##  small polarization components -- it can fire while a component orders of
##  magnitude below the dominant one is still evolving.  If you care about such
##  a component, disable it and use a fixed end_time instead.
##
###############################################################

##  Gradient (wall) energy coefficients.
##  The unrotated WallEnergyDerivative/WallEnergy take G110 together with the
##  ratios G11/G110 etc., NOT the bare G11/G12/G44 that the Rotated variants use.
##  Matching the two residuals term-by-term at twist = 0 gives
##      G110 = G11,  G11_G110 = 1,  G12_G110 = G12/G11,  G44_G110 = G44/G11
##  and, crucially, G44P_G110 = 0: the Rotated kernel carries NO antisymmetric
##  G44' term, so copying a G-block from another unrotated input (surface_charge.i
##  ships G44P_G110 = 0.3) would silently add physics this model never had.
##  Here G11 = 0.5, G12 = -0.02, G44 = 0.02.
G110      = 0.5
G11_G110  = 1.0
G12_G110  = -0.04
G44_G110  = 0.04
G44P_G110 = 0.0

L = 0.25
n_xy = 1
Z = 110.0
n_z = 440.0

period = 0.028559933214452663


[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${n_z}
    ny = ${n_xy}
    nz = ${n_xy}
    xmin = -${Z}
    xmax = ${Z}
    ymin = -${L}
    ymax = ${L}
    zmin = -${L}
    zmax = ${L}
    elem_type = HEX8
  []
  [./cnode]
    input = gen
    type = ExtraNodesetGenerator
    coord = '-${Z} -${L} -${L}'
    new_boundary = 100
  [../]

[]

[GlobalParams]
  #len_scale = 1.0

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z

  displacements = 'u_x u_y u_z'
  
  potential_E_int = potential_E_int
[]

[Functions]
  [./DW_func]
    type = ParsedFunction
    symbol_names = 'A p'
    symbol_values = '0.2 ${period}'
    expression = 'A*cos(p*x)'
  [../]
[]

[Variables]
  [./u_x]
  [../]
  [./u_y]
  [../]
  [./u_z]
  [../]
  [./global_strain]
    order = SIXTH
    family = SCALAR
  [../]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -1e-5
      max = 1e-5
      seed = 1
    [../]
    block = '0'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -1e-5
      max = 1e-5
      seed = 2
    [../]
    block = '0'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = DW_func
    [../]
    block = '0'
  [../]
  [./potential_E_int]
    order = FIRST
    family = LAGRANGE
    block = '0'
  [../]
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [./all_FE]
        strain = SMALL
        add_variables = true
        incremental = false
        eigenstrain_names = 'total_eigenstrain'
        global_strain = global_strain
        generate_output = 'strain_xx strain_yy strain_zz strain_yz strain_xz strain_xy'
	block = '0'
      [../]
    []

    # GlobalStrain action for generating the objects associated with the global
    # strain calculation and associated displacement visualization

    [./GlobalStrain]
      [./global_strain]
        scalar_global_strain = global_strain
        displacements = 'u_x u_y u_z'
        auxiliary_displacements = 'disp_x disp_y disp_z'
        global_displacements = 'ug_x ug_y ug_z'
      [../]
    [../]
  []
[]


[Kernels]

  ### Operators for the polar field: ###
  [./bed_x]
    type = BulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
    block = '0'
  [../]
  [./bed_y]
    type = BulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
    block = '0'
  [../]
  [./bed_z]
    type = BulkEnergyDerivativeEighth
    variable = polar_z
    component = 2
    block = '0'
  [../]

  [./walled_x]
    type = WallEnergyDerivative
    variable = polar_x
    component = 0
    block = '0'
  [../]
  [./walled_y]
    type = WallEnergyDerivative
    variable = polar_y
    component = 1
    block = '0'
  [../]
  [./walled_z]
    type = WallEnergyDerivative
    variable = polar_z
    component = 2
    block = '0'
  [../]
 
  [./elastic_polar_coupled_x]
    type = CubicParentElasticPDerivative
    variable = polar_x
    component = 0
    block = '0'
  [../]
  [./elastic_polar_coupled_y]
    type = CubicParentElasticPDerivative
    variable = polar_y
    component = 1
    block = '0'
  [../]
  [./elastic_polar_coupled_z]
    type = CubicParentElasticPDerivative
    variable = polar_z
    component = 2
    block = '0'
  [../]
  
  [./polar_x_electric_E]
    type = PolarElectricEStrong
    variable = potential_E_int
    block = '0'
  [../]
  [./FE_E_int]
    type = Electrostatics
    variable = potential_E_int
    block = '0'
  [../]

  [./polar_electric_px]
    type = PolarElectricPStrong
    variable = polar_x
    component = 0
    block = '0'
  [../]
  [./polar_electric_py]
    type = PolarElectricPStrong
    variable = polar_y
    component = 1
    block = '0'
  [../]
  [./polar_electric_pz]
    type = PolarElectricPStrong
    variable = polar_z
    component = 2
    block = '0'
  [../]
  
  [./polar_x_time]
    type = TimeDerivativeScaled
    variable=polar_x
    time_scale = 1.0
    block = '0'
  [../]
  [./polar_y_time]
    type = TimeDerivativeScaled
    variable=polar_y
    time_scale = 1.0
    block = '0'
  [../]
  [./polar_z_time]
    type = TimeDerivativeScaled
    variable = polar_z
    time_scale = 1.0
    block = '0'
  [../]
[]


[Materials]

  [./eigen_strain]
    type = ComputeEigenstrain
    #xx yy zz (xy yz xz)
    eigen_base = '0.0 0.0 0.0 0.0 0.0 0.0'
    eigenstrain_name = 'epitaxy'
    block = '0'
  [../]

  [./Landau_P]
    type = GenericConstantMaterial
    prop_names = 'alpha1 alpha11 alpha12 alpha111 alpha112 alpha123 alpha1111 alpha1112 alpha1122 alpha1123'
    prop_values = '-0.03635 -0.2097 0.7974 1.294 -1.95 -2.5 38.63 25.29 16.37 13.67'
    block = '0'
  [../]
  
  [./Landau0_G_FE]
    type = GenericConstantMaterial
    prop_names = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '${G110} ${G11_G110} ${G12_G110} ${G44_G110} ${G44P_G110}'
    block = '0'
  [../]

  [./mat_C]
    type = GenericConstantMaterial
    prop_names = 'C11 C12 C44'
    prop_values = '178.0 96.4 122.0'
    block = '0'
  [../]
 
  [./mat_Q] #this is slightly different from Phys Rev B 89, 174111 (2014) ?
    type = GenericConstantMaterial
    prop_names = 'Q11 Q12 Q44'
    prop_values = '0.10 -0.045 0.029' #0.059/4?
    block = '0'
  [../]
  
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '178.0 96.4 96.4 178.0 96.4 178.0 122.0 122.0 122.0'
    block = '0'
  [../]
  
  [./film_eigenstrain]
    type = CompositeEigenstrain
    ##  NOTE: 'flexo' carried weight2 = -1, so dropping it is not a no-op.
    tensors = 'ferro epitaxy'
    weights = 'weight1 weight3'
    eigenstrain_name = total_eigenstrain
    coupled_variables = 'polar_x polar_y polar_z'
    block = '0'
  [../]
  [./weight1]
    type = DerivativeParsedMaterial
    block = '0'
    expression = '1'
    property_name = weight1
    coupled_variables = 'polar_x polar_y polar_z'
  [../]
  [./weight3]
    type = DerivativeParsedMaterial
    block = '0'
    expression = '1'
    property_name = weight3
  [../]
  
  [./stress_1]
    type = ComputeLinearElasticStress
    block = '0'
  [../]

  [./electrostrictive_eigenstrain]
    type = ComputeCubicParentElectrostrictiveStrain
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
    eigenstrain_name = 'ferro'
    block = '0'
  [../]
  [./permitivitty]

    ###############################################
    ##
    ##  BTO background permittivity 45 (used in flexo paper)
    ##
    ###############################################

    type = GenericConstantMaterial
    prop_names = 'permittivity'
    prop_values = '0.39'
    block = '0'
  [../]
[]

[Postprocessors]
  [./dt]
     type = TimestepSize
  [../]
  [./FbP]
    type = BulkEnergyEighth
    execute_on = 'timestep_end'
    block = '0'
  [../]
  [./FgP]
    type = WallEnergy
    execute_on = 'timestep_end'
    block = '0'
  [../]
  
  [./Fele]
    type = ElectrostaticEnergy
    execute_on = 'timestep_end'
    block = '0'
  [../]
  
  [./Fela]
    type = CubicParentElasticEnergy
    execute_on = 'timestep_end'
    polar_x = polar_x
    polar_y = polar_y
    polar_z = polar_z
    block = '0'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor
    pp_names = 'FbP FgP Fela Fele'
    pp_coefs = ' 1 1 1 1'
    execute_on = 'timestep_end'

    ##########################################
    #
    # NOTE: Ferret output is in attojoules
    #
    ##########################################
  [../]
  
  [./maxPy]
     type = NodalExtremeValue
     variable = polar_y
  [../]
  [./maxPx]
     type = NodalExtremeValue
     variable = polar_x
  [../]
  
  [./perc_change]
    type = EnergyRatePostprocessor
    postprocessor = Ftot
    execute_on = 'timestep_end'
    dt = dt
  [../]
[]


[BCs]
  [./Periodic]
    [./x]
      auto_direction = 'x y z'
      variable = 'u_x u_y u_z polar_x polar_y polar_z potential_E_int'
    [../]
  [../]
  # fix center point location
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

[UserObjects]
  [./kill]
   type = Terminator
   expression = 'perc_change <= 5.0e-7'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -snes_atol  -snes_rtol -ksp_rtol -pc_type'
    petsc_options_value = '    121            1e-8          1e-6      1e-5    bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  scheme = 'bdf2'
  dtmin = 1e-13
  dtmax = 10.0

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 25  #usually 10
    linear_iteration_ratio = 100
    dt = 0.1
    growth_factor = 1.1
  [../]
[]

[Outputs]
  print_linear_residuals = false
  perf_graph_live = false
  [./out]
    type = Exodus
    execute_on = 'INITIAL FINAL'
    file_base = out_ferroelectric_domain_wall
    elemental_as_nodal = true
  [../]
  [./csv]
    type = CSV
    file_base = out_ferroelectric_domain_wall
  [../]
[]

