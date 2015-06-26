[Mesh]
  file = out_PbTiO3_T_0_E_0_G110_6e-1_size_160.e
  uniform_refine=0
[]

[Variables]
  [./polar_x]
    order = FIRST
    family = LAGRANGE
    block='2'
    initial_from_file_var = polar_x
    initial_from_file_timestep = 45
  [../]
  [./polar_y]
    order = FIRST
    family = LAGRANGE
    block='2'
    initial_from_file_var = polar_y
    initial_from_file_timestep = 45
  [../]
  [./polar_z]
    order = FIRST
    family = LAGRANGE
    block='2'
    initial_from_file_var = polar_z
    initial_from_file_timestep = 45
  [../]
  [./potential_int]
    order=FIRST
    family = LAGRANGE
    initial_from_file_var = potential_int
    initial_from_file_timestep = 45
  [../]
  [./potential_ext]
    order=FIRST
    family = LAGRANGE
    initial_from_file_var = potential_ext
    initial_from_file_timestep = 45
  [../]
[]

[AuxVariables]
  [./surface_charge]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ex]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ey]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ez]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[GlobalParams]
   len_scale=1e-9
   #len_scale=1.0
   alpha1=-1.8202e8 # 3.8(T-479)*10^5 C^{-2}m^2 (T=0 K)
   alpha11=-7.3e7
   alpha111=2.6e8
   alpha12=7.5e8
   alpha112=6.1e8
   alpha123=-3.7e9
   G110=0.6e-10
#   #G110=0.0
   G11/G110=0.6
   G12/G110=0.0
   G44/G110=0.3
   G44P/G110=0.3
   permittivity=8.85e-12
   polar_x=polar_x
   polar_y=polar_y
   polar_z=polar_z
   potential_int=potential_int
   potential_ext=potential_ext
   time_scale = 1.0e-31

# this ostensibly defines dt = 1e-14. Note if time_scale=1 then dt=1e17 to see domain relaxation
# Olle suggests that the dynamics should occur on the 10 ps time scale. This parameter allows the time scale to be
# user-input and most likely depends on material and problem in question, maybe even temperature.
[]

[AuxKernels]
   [./surfcharge_3]
     type = SurfaceChargeAux
     variable = surface_charge
     boundary = '3'
   [../]
   [./surfcharge_4]
     type = SurfaceChargeAux
     variable = surface_charge
     boundary = '4'
   [../]
   [./surfcharge_5]
     type = SurfaceChargeAux
     variable = surface_charge
     boundary = '5'
   [../]
   [./surfcharge_6]
     type = SurfaceChargeAux
     variable = surface_charge
     boundary = '6'
   [../]
   [./surfcharge_7]
     type = SurfaceChargeAux
     variable = surface_charge
     boundary = '7'
   [../]
   [./surfcharge_8]
     type = SurfaceChargeAux
     variable = surface_charge
     boundary = '8'
   [../]
   [./field_x]
     type = ExFieldAux
     variable = Ex
   [../]
   [./field_y]
     type = EyFieldAux
     variable = Ey
   [../]
   [./field_z]
     type = EzFieldAux
     variable = Ez
   [../]
[]

[Kernels]
  [./bed_x]
    type = BulkEnergyDerivative
    variable = polar_x
    component=0
  [../]
  [./bed_y]
    type = BulkEnergyDerivative
    variable = polar_y
    component=1
  [../]
  [./bed_z]
    type = BulkEnergyDerivative
    variable = polar_z
    component=2
  [../]
  [./walled_x]
     type=WallEnergyDerivative
     variable=polar_x
     component=0
  [../]
  [./walled_y]
     type=WallEnergyDerivative
     variable=polar_y
     component=1
  [../]
  [./walled_z]
     type=WallEnergyDerivative
     variable=polar_z
     component=2
  [../]
  [./polar_electric_E]
     type=PolarElectricEStrong
     variable=potential_int
     block='2'

  [../]
  [./FE_E_int]
     type=Electrostatics
     variable=potential_int
     block='2'
  [../]
  [./E_int]
     type=Electrostatics
     variable=potential_int
     block='1'
  [../]
  [./E_ext]
     type=Electrostatics
     variable=potential_ext
     block='1'
  [../]
  [./FE_E_ext]
     type=Electrostatics
     variable=potential_ext
     block='2'
  [../]
  [./polar_electric_px]
     type=PolarElectricPStrong
     variable=polar_x
     component=0
  [../]
  [./polar_electric_py]
     type=PolarElectricPStrong
     variable=polar_y
     component=1
  [../]
  [./polar_electric_pz]
     type=PolarElectricPStrong
     variable=polar_z
     component=2
  [../]
  [./polar_x_time]
     type=TimeDerivativeScaled
     variable=polar_x
  [../]
  [./polar_y_time]
     type=TimeDerivativeScaled
     variable=polar_y
  [../]
  [./polar_z_time]
     type=TimeDerivativeScaled
     variable=polar_z
  [../]
[]

#[ICs]
#  #IC cannot be larger than Pmax or BulkEnergy diverges after sufficient num_steps.
#  #in fact, it seems that if the initial guess is close to Pmax, then the solution
#  #does not waste time steps
#  [./polar_x_constic_rand]
#     type=RandomIC
#     variable=polar_x
#     min = 0.6
#     max = 0.75
#  [../]
#  [./polar_y_constic_rand]
#     type=RandomIC
#     variable=polar_y
#     min = 0.6
#     max = 0.75
#  [../]
#  [./polar_z_constic_rand]
#     type=RandomIC
#     variable=polar_z
#     min = 0.0
#     max = 0.0
#  [../]
#[]

[BCs]
  [./potential_int_upz]
    type = DirichletBC
    variable = potential_int
    boundary = '1'
    value = 0.0
  [../]
  [./potential_int_downz]
    type = DirichletBC
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
    value = 1 #note that this value depends on distance between plates... 2 V / 160 nm = 125 kV/cm
    #value = 0.0
    #implicit=false
  [../]
  [./potential_ext_downz]
    type = PresetBC
    variable = potential_ext
    boundary = '2'
    value = -1
    #value = 0.0
  [../]
[]

[Postprocessors]
   [./bulk_energy]
    type=BulkEnergy
   [../]
   [./wall_energy]
    type=WallEnergy
   [../]
   [./electrostatic_energy]
    type=ElectrostaticEnergy
   [../]
   [./total_energy]
    type=TotalEnergy
    bulk_energy=bulk_energy
    wall_energy=wall_energy
    electrostatic_energy=electrostatic_energy
    elastic_energy=elastic_energy
   [../]
   [./_pps_percent]
     type = PercentChangePostprocessor
     postprocessor = total_energy
   [../]

[]

[UserObjects]
  [./kill]
    type = Terminator
    expression = '_pps_percent <= 4.0e-7'
  [../]
[]

[Executioner]
  type=Transient

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.1e-15
    optimal_iterations = 3
    growth_factor = 1.001
    cutback_factor =  0.999
  [../]

  scheme = 'explicit-euler'   #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"

  dtmin=1.0e-29
  dtmax=1.70e-15

  num_steps=35000
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname='-ksp_type -snes_type  -snes_rtol -ksp_rtol -pc_type -pc_hypre_boomeramg_strong_threshold -snes_linesearch_type -pc_factor_zeropivot'
  petsc_options_value=' gmres     newtonls       1e-8     1e-6      hypre     0.5                                 basic         1e-50  '
[]


[Outputs]

  print_linear_residuals = true
  print_perf_log = true

  [./out]
    type = Exodus
    file_base = out_PbTiO3_T_0_E_100_G110_6e-1_size_160
    output_initial = true
    elemental_as_nodal = false
    interval = 1
  [../]

[]
