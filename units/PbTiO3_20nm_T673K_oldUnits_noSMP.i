[Mesh]
  file = slab_exodus_coarse_20.e
  uniform_refine = 0
[]

[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block='2'
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block='2'
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block='2'
  [../]
  [./potential_int]
    order=FIRST
    family = LAGRANGE
  [../]
  [./potential_ext]
    order=FIRST
    family = LAGRANGE
  [../]
[]


[GlobalParams]
  len_scale = 1e-9
  #len_scale=1.0
  alpha1 = -3.46849e7 # 3.8(T-479)*10^5 C^{-2}m^2 (T=0 K)
  alpha11 = -7.3e7
  alpha111 = 2.6e8
  alpha12 = 7.5e8
  alpha112 = 6.1e8
  alpha123 = -3.7e9
  G110 = 2.5e-10
#   #G110=0.0
  G11/G110 = 0.6
  G12/G110 = 0.0
  G44/G110 = 0.3
  G44P/G110 = 0.3
  polar_x = polar_x
  polar_y = polar_y
  polar_z = polar_z
  potential_int=potential_int
  potential_ext=potential_ext
  time_scale = 1.0e-31

# this ostensibly defines dt = 1e-14. Note if time_scale=1 then dt=1e17 to see domain relaxation
# Olle suggests that the dynamics should occur on the 10 ps time scale. This parameter allows the time scale to be
# user-input and most likely depends on material and problem in question, maybe even temperature.
[]

[Kernels]
  [./bed_x]
    type = BulkEnergyDerivative
    variable = polar_x
    component = 0
    implicit = false
  [../]
  [./bed_y]
    type = BulkEnergyDerivative
    variable = polar_y
    component = 1
    implicit = false
  [../]
  [./bed_z]
    type = BulkEnergyDerivative
    variable = polar_z
    component = 2
    implicit = false
  [../]
  [./walled_x]
     type = WallEnergyDerivative
     variable = polar_x
     component = 0
  [../]
  [./walled_y]
     type=WallEnergyDerivative
     variable = polar_y
     component = 1
  [../]
  [./walled_z]
     type = WallEnergyDerivative
     variable = polar_z
     component = 2
  [../]
  [./polar_electric_E]
     type=PolarElectricEStrong
     variable = potential_int
     block = '2'
     implicit = false
  [../]
  [./FE_E_int]
     type=Electrostatics
     variable = potential_int
     block = '2'
  [../]
  [./E_int]
     type=Electrostatics
     variable=potential_int
     block='1'
  [../]
  [./E_ext]
     type=Electrostatics
     variable = potential_ext
     block= '1'
  [../]
  [./FE_E_ext]
     type = Electrostatics
     variable = potential_ext
     block = '2'
  [../]
  [./polar_electric_px]
     type=PolarElectricPStrong
     variable=polar_x
     component = 0
     implicit = false
  [../]
  [./polar_electric_py]
     type =PolarElectricPStrong
     variable=polar_y
     component = 1
     implicit = false
  [../]
  [./polar_electric_pz]
     type = PolarElectricPStrong
     variable=polar_z
     component = 2
     implicit = false
  [../]
  [./polar_x_time]
     type = TimeDerivativeScaled
     variable = polar_x
  [../]
  [./polar_y_time]
     type = TimeDerivativeScaled
     variable = polar_y
  [../]
  [./polar_z_time]
     type = TimeDerivativeScaled
     variable = polar_z
  [../]
[]

[ICs]
  #IC cannot be larger than Pmax or BulkEnergy diverges after sufficient num_steps.
  #in fact, it seems that if the initial guess is close to Pmax, then the solution
  #does not waste time steps
  [./polar_x_constic_rand]
     type=RandomIC
     variable=polar_x
     min = 0.6
     max = 0.75
     seed = 1
  [../]
  [./polar_y_constic_rand]
     type=RandomIC
     variable=polar_y
     min = 0.6
     max = 0.75
     seed = 1
  [../]
  [./polar_z_constic_rand]
     type=ConstantIC
     variable=polar_z
     value = 0.0
  [../]
[]

[BCs]
  [./potential_int_upz]
    type = PresetBC
    variable = potential_int
    boundary = '1'
    value = 0.0
  [../]
  [./potential_int_downz]
    type = PresetBC
    variable = potential_int
    boundary = '2'
    value = 0.0
  [../]
  # Applied field: for zero field use NeumannBC on the external potential = 0. A
  # Note that \nabla^2 \Phi_{ext} = 0 is satisfied if \Phi_{ext} = 0, ie Dirichlet and Neumann BC classes are equivalent
  [./potential_ext_upz]
    type = PresetBC
    variable = potential_ext
    boundary = '1'
    #value = 5.0e-11 #1e-12 gives a field of ~1e8 V/m ~ 1000 kV/cm for size 1, what about size 100?
                    # seems that this doesn't matter here? dt needs to change for the size however...
                    # lowering this value, to about 1e-16 is the limit here. 5e-11 gives 5.6e9 V/m for 100 size
    value = 0.0
    #implicit=false
  [../]
  [./potential_ext_downz]
    type = PresetBC
    variable = potential_ext
    boundary = '2'
    #value = -5.0e-11
    value = 0.0
  [../]
[]

[Postprocessors]
   [./bulk_energy]
     type=BulkEnergy
     execute_on = 'timestep_end'
   [../]
   [./wall_energy]
     type=WallEnergy
     execute_on = 'timestep_end'
   [../]
   [./electrostatic_energy]
     type=ElectrostaticEnergy
     execute_on = 'timestep_end'
   [../]
   [./total_energy]
     type=TotalEnergy
     bulk_energy = bulk_energy
     wall_energy = wall_energy
     electrostatic_energy = electrostatic_energy
     elastic_energy = 0.0
     coupled_energy = 0.0
     execute_on = 'timestep_end'
   [../]
  # [./_dt]
  #   type = TimestepSize
  # [../]
  # [./_pps_percent]
  #   type = PercentChangePostprocessor
  #   postprocessor = TotalEnergy
  # [../]
  # [./gradx]
  #   type=TotalEnergyGradient
  #   component = 0
  # [../]
  # [./grady]
  #   type=TotalEnergyGradient
  #   component = 1
  # [../]
  # [./gradz]
  #   type=TotalEnergyGradient
  #   component = 2
  # [../]
  # [./GradE-L2]
  #   type=TotalEnergyGradientL2
  #   gradx = gradx
  #   grady = grady
  #   gradz = gradz
  # [../]
  # [./R(i)]
  #   type=Residual
  # [../]
[]

[Executioner]
  type=Transient

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0e-15
    #timestep_limiting_function = total_energy
    optimal_iterations = 3
    growth_factor = 1.001
    cutback_factor =  0.999
  [../]
  # have to cap dt at 2e-15 since energy will diverge at larger dt at any time step
  # wall energy acts as a damper which forces use to have large num steps


  #[./TimeStepper] #mildly faster than ConstantDT but can cause divergence
  #  type = SolutionTimeAdaptiveDT
  #  dt = 1e-16
  #  dtmin=1e-25
  #[../]
  scheme = 'explicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  #dt=1e-15 #adjustable within a few orders of magnitude. It seems that we need small dt steps but large num_steps.
           #also, seems that 1e-15 is largest time step to be able to be done
  dtmin = 1.0e-25
  dtmax = 1.2e-15
  #NOTE: First time step calculates the depolarization field due to the unphysical initial condition. Energy may increase, which is allowed.
  #presetBC included to reduce this
  num_steps = 15000
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname='-ksp_type -snes_type  -snes_rtol -ksp_rtol -pc_type  -snes_linesearch_type -pc_factor_zeropivot'
  petsc_options_value=' gmres     newtonls         1e-8     1e-8     bjacobi     basic         1e-50  '
[]

#seems that Boomer AMG is good for electrostatics
#-snes_rtol -ksp_rtol
#1e-5     1e-8
#[./Adaptivity]
#    marker = 'marker_x marker_y'
##    cycles_per_step = 2
##    max_h_level = 2
#    [./Indicators]
#      [./indicator_x]
#        type = GradientJumpIndicator
#        variable = polar_x
#      [../]
#      [./indicator_y]
#        type = GradientJumpIndicator
#        variable = polar_y
#      [../]
#    [../]
#    [./Markers]
#      [./marker_x]
#        type = ErrorFractionMarker
#        indicator = indicator_x
#        coarsen = 0.01
#        refine = 0.01
#      [../]
#      [./marker_y]
#        type = ErrorFractionMarker
#        indicator = indicator_y
#        coarsen = 0.01
#        refine = 0.01
#      [../]
#    [../]
#[../]

[Outputs]

print_linear_residuals = true
print_perf_log = true

  [./out]
    type = Exodus
    file_base = out_PbTiO3_T673K_20nm_oldUnits_noSMP
    output_initial = true
    elemental_as_nodal = false
    interval = 1
  [../]

  [./debug]
    type = VariableResidualNormsDebugOutput
  [../]
[]
