
[Mesh]
  type = GeneratedMesh
  dim = 2

  ############################################
  ##
  ##  Grid definition. Note that it should be
  ##  nJ = 2*(Jmax-Jmin) for J = x, y
  ##
  ############################################

  nx = 60
  ny = 60

  ############################################
  ##
  ## Actual spatial coordinates of mesh. 
  ## Jmax - Jmin = nJ/2 for J = x, y
  ##  
  ############################################

  xmin = -15
  xmax = 15
  ymin = -15
  ymax = 15

[]

[GlobalParams]

  len_scale = 1.0

  ############################################
  ##
  ## BST Landau coefficients from
  ## Y.H. Huang and L.-Q. Chen et al
  ##    Appl. Phys. Lett. 112, 102901 (2018)
  ##
  ############################################

                        #----------------#
                        #     Units      #
                        #----------------#
  T = 210               # Kelvin


  #concentration dependent coefficients

  alpha01 = 0.08        # nm^2 nN / aC^2
  alpha011 = -0.1154    # nm^6 nN / aC^4
  alpha012 = 0.653      # nm^6 nN / aC^4
  alpha0111 = -2.106    # nm^10 nN / aC^6
  alpha0112 = 4.091     # nm^10 nN / aC^6
  alpha0123 = -6.688    # nm^10 nN / aC^6

  #fixed coefficients 

  alpha1111 = 75.9      # nm^14 nN / aC^8
  alpha1112 = -21.93    # nm^14 nN / aC^8
  alpha1122 = -22.21    # nm^14 nN / aC^8
  alpha1123 = 24.16     # nm^14 nN / aC^8

  # coupling constants 
  b1 = -1.75
  b2 = 1.05
  b3 = 0.483

  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z



  ############################################
  ##
  ## Fourier Noise implementation from 
  ## D. Schwen (INL) - see paper
  ##
  ############################################


  x = conc
 
  lambda = 5.0          # correlation length for random field (nm)
  range = 0.01          # variation of the Sr concentration (plus or minus)
  mid = 0.4             # nominal value of the Sr concentration field

  potential_E_int = potential_int

  seed = 1  #NOTE THAT THE SAME MESH GRID MUST BE USED WITH THE SAME PROCESSOR GROUPING ALWAYS TO COMPARE DOMAIN STRUCTURE EVOLUTION
[]


[Variables]

  [./polar_x]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-5
      max = 0.1e-5
    [../]
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-5
      max = 0.1e-5
    [../]
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = RandomIC
      min = -0.1e-5
      max = 0.1e-5
    [../]
  [../]
  [./potential_int]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Functions]
  [./fn]  
    type = SDFourierNoise

    # Note that FourierNoise automatically makes concentration periodic which in turn makes the Landau coefficients periodic.
    # This then makes the polarization periodic.

  [../]

  ##############################
  ##
  ## Define the electric field
  ## expression to be used below
  ##
  ##############################

  #[./bc_func_1]
  #  type = ParsedFunction
  #  value = 'amplitude*sin(freq*t)'
  #  vars = 'freq amplitude'
  #  vals = '${freq}  ${amplitude}'
  #[../]
[]

[AuxVariables]
  [./conc]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = FunctionIC
      function = fn
    [../]
  [../]

  #[./Ey]
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]

  [./fb]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
 # [./cEy]
 #   type = ElecFieldAux
 #   component = 2
 #   variable = Ey
 # [../]
  [./fbk]
    type = SDBulkEnergyDensity
    variable = fb
  [../]
[]

[Materials]
  ############################################
  ##
  ## Gradient coefficients (assumed BTO) from
  ## Marton and Hlinka
  ##    Phys. Rev. B. 74, 104014, (2006)
  ##
  ############################################

  [./Landau_G]
    type = GenericConstantMaterial
    prop_names = 'G110 G11_G110 G12_G110 G44_G110 G44P_G110'
    prop_values = '1.0 0.51 -0.02 0.02 0.0'
  [../]
[]


[Kernels]

  #########################################
  ##
  ## Landau's problem (no elastic coupling:
  ##
  #########################################

  [./bed_x]
    type = SDBulkEnergyDerivativeEighth
    variable = polar_x
    component = 0
  [../]
  [./bed_y]
    type = SDBulkEnergyDerivativeEighth
    variable = polar_y
    component = 1
  [../]
  [./bed_z]
    type = SDBulkEnergyDerivativeEighth
    variable = polar_z
    component = 2
  [../]
  [./walled_x]
    type = WallEnergyDerivative
    variable = polar_x
    component = 0
 [../]
 [./walled_y]
    type = WallEnergyDerivative
    variable = polar_y
    component = 1
  [../]
  [./walled_z]
     type = WallEnergyDerivative
     variable = polar_z
     component = 2
  [../]

  #########################################
  ##
  ## Poisson's equation and P*E interaction
  ##
  #########################################

  [./polar_x_electric_E]
     type = PolarElectricEStrong
     variable = potential_int
  [../]
  [./FE_E_int]
     type = Electrostatics
     variable = potential_int
     permittivity = 0.0637501464
     #NOTE: This is a static permittivity contribution from core-electrons.
     #      This effectively screens the electrostatic interactions.
     #      For BTO, this value is about 7*e_0, where e_0 is the permitivitty of the vacuum.
     #      See the brief discussion before sec. IV on pp. 4 of Phys. Rev. B. 74, 104014, (2006)
  [../]
  [./polar_electric_px]
     type = PolarElectricPStrong
     variable = polar_x
     component = 0
  [../]
  [./polar_electric_py]
     type = PolarElectricPStrong
     variable = polar_y
     component = 1
  [../]
  [./polar_electric_pz]
     type = PolarElectricPStrong
     variable = polar_z
     component = 2
  [../]

  #########################################
  ##
  ## Time dependence
  ##
  #########################################

  [./polar_x_time]
     type = TimeDerivativeScaled
     variable = polar_x
     # Time scale estimate for BTO, from Hlinka (2007)
     # We use seconds here
     time_scale = 1.0
  [../]
  [./polar_y_time]
     type = TimeDerivativeScaled
     variable = polar_y
     time_scale = 1.0
  [../]
  [./polar_z_time]
     type = TimeDerivativeScaled
     variable = polar_z
     time_scale = 1.0
  [../]
[]


[BCs]

  ##############################
  ##
  ## Boundary Condition System
  ##
  ## This corresponds to a small
  ## ac field along the direction
  ## of the spontaneous polarization
  ## that was ostensibly put there
  ## by a strong dc bias.
  ##
  ##############################

#  [./front_pot]
#    type = FunctionDirichletBC
#    variable = potential_int
#    boundary = 'top'
#    function = bc_func_1
#  [../]
#  [./back_pot]
#    type = DirichletBC
#    variable = potential_int
#    boundary = 'bottom'
#    value = 0.0
#  [../]

  [./top_pot]
    type = DirichletBC
    variable = potential_int
    boundary = 'top'
    value = 0.0
  [../]

  [./back_pot]
    type = DirichletBC
    variable = potential_int
    boundary = 'bottom'
    value = 0.0
  [../]


  [./Periodic]
    [./xy]
      auto_direction = 'x y'
      variable = 'polar_x polar_y polar_z'
    [../]
  [../]

[]

[Postprocessors]
  #[./avePy]
  #  type = ElementAverageValue
  #  variable = polar_y
  #  execute_on = 'initial timestep_end'
  #[../]
  #[./cPs]
  #  type = ElementAverageValue
  #  variable = Ps
  #  execute_on = 'initial timestep_end'
  #[../]
  #[./Ea]
  #  type = ElementAverageValue
  #  variable = Ey
  #  execute_on = 'initial timestep_end'
  #[../]
  #[./inducedP]
  #  type = LinearCombinationPostprocessor
  #  pp_names = 'avePy cPs'
  #  pp_coefs = ' 1 -1'
  #  execute_on = 'initial timestep_end'
  #[../]

  [./Fgrad]
    type = WallEnergy
    execute_on = 'initial timestep_end'
  [../]

  [./Fsdbulk]
    type = SDBulkEnergyEighth
    execute_on = 'initial timestep_end'
    block = '0'
  [../]
  [./Ftot]
    type = LinearCombinationPostprocessor 
    pp_names = 'Fgrad Fsdbulk'
    pp_coefs = ' 1 1 ' 
    execute_on = 'initial timestep_end'
  [../]
  [./perc_change]
    type = PercentChangePostprocessor
    postprocessor = Ftot
    execute_on = 'initial timestep_end'
  [../]
[]

[UserObjects]
  [./kill]
   type = Terminator
   expression = 'perc_change <= 5.0e-5'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_gmres_restart  -snes_atol -ksp_rtol -pc_type'
    petsc_options_value = '    121                1e-10     1e-8     bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  scheme = 'implicit-euler'
  dtmin = 1e-10
  dt = 3.0
  dtmax = 10.0
  num_steps = 2
[]

[Outputs]
  print_linear_residuals = false
  [./out]
    type = Exodus
    file_base = out_test
    elemental_as_nodal = true
  [../]
  [./outCSV]
    type = CSV
    new_row_tolerance = 1e-16
    file_base = out_test
  [../]
[]
