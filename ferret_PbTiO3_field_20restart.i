[Mesh]
  file = out_PbTiO3_bulk_len-9_icConst_n3500_0field_1wall_restart.e
  uniform_refine=0
  distribution=serial
[]

[Variables]
#active='polar_x polar_y polar_z'
#active='potential_int potential_ext'
  [./polar_x]
    #scaling=1e14
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_x
    initial_from_file_timestep = 120
    block='2'
  [../]
  [./polar_y]
    #scaling=1e14
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_y
    initial_from_file_timestep = 120
    block='2'
  [../]
  [./polar_z]
    #scaling=1e14
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = polar_z
    initial_from_file_timestep = 120
    block='2'
  [../]
  [./potential_int]
    #scaling=1e7
    order=FIRST
    family = LAGRANGE
    initial_from_file_var = potential_int 
    initial_from_file_timestep = 120
  [../]
  [./potential_ext]
    #scaling=1e7
    order=FIRST
    family = LAGRANGE
    initial_from_file_var = potential_ext 
    initial_from_file_timestep = 120
  [../]
[]

[AuxVariables]
[./auxv_es_energy_density_e] #es for electrostatic
     order=CONSTANT
     family=MONOMIAL
     initial_from_file_var = auxv_es_energy_density_e
     initial_from_file_timestep = 120
  [../]
  [./auxv_es_energy_density]
     order=CONSTANT
     family=MONOMIAL
     initial_from_file_var = auxv_es_energy_density
     initial_from_file_timestep = 120
  [../]
  [./auxv_bulk_energy_density]
     order=CONSTANT
     family=MONOMIAL
     initial_from_file_var = auxv_bulk_energy_density
     initial_from_file_timestep = 120
  [../]
  [./fieldEz]
     order=CONSTANT
     family=MONOMIAL
     initial_from_file_var = fieldEz
     initial_from_file_timestep = 120
  [../]
[]

[GlobalParams]
   len_scale=1e-9
   #len_scale=1.0
   alpha1=-1.7252e8 # 3.8(T-479)*10^5 C^{-2}m^2
   alpha11=-7.3e7
   alpha111=2.6e8
   alpha12=7.5e8
   alpha112=6.1e8
   alpha123=-3.7e9
   G110=1e-10
#   #G110=0.0
   G11/G110=0.6
   G12/G110=0.0
   G44/G110=0.3
   G44P/G110=0.3
   permittivity=8.85e-12
   #permittivity=1.0
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
  #active='diff'
  [./auxk_electrostatic_energy_density_e]
    type =ElectrostaticEnergyDensityE
    variable =auxv_es_energy_density_e
    potential=potential_int
  [../]
  [./auxk_electrostatic_energy_density]
    type =ElectrostaticEnergyDensity
    variable =auxv_es_energy_density
  [../]
  [./auxk_bulk_energy_density]
    type =BulkEnergyDensity
    variable =auxv_bulk_energy_density
  [../]
  [./Ez_fieldaux]
    type = Ez_fieldAux
    variable = fieldEz
    potential_ext = potential_ext
    potential_int = potential_int
  [../]
#  [./Ex_fieldaux]
#    type = Ex_fieldAux
#    variable = fieldEx
#    potential_ext = potential_ext
#    potential_int = potential_int
#  [../]
#  [./Ey_fieldaux]
#    type = Ey_fieldAux
#    variable = fieldEy
#    potential_ext = potential_ext
#    potential_int = potential_int
#  [../]
[]

[Kernels]
  [./bed_x]
    type = BulkEnergyDerivative
    variable = polar_x
    component=0
    implicit=false
  [../]
  [./bed_y]
    type = BulkEnergyDerivative
    variable = polar_y
    component=1
    implicit=false
  [../]
  [./bed_z]
    type = BulkEnergyDerivative
    variable = polar_z
    component=2
    implicit=false
  [../]
  [./walled_x]
     type=WallEnergyDerivative_scaled
     variable=polar_x
     component=0
  [../]
  [./walled_y]
     type=WallEnergyDerivative_scaled
     variable=polar_y
     component=1
  [../]
  [./walled_z]
     type=WallEnergyDerivative_scaled
     variable=polar_z
     component=2
  [../]
  [./polar_electric_E]
     type=PolarElectricEStrong
     variable=potential_int
     block='2'
     implicit=false
  [../]
  [./diffusion_E]
     type=Electrostatics
     #type=Diffusion
     variable=potential_int
     block='1 2'
  [../]
  [./diffusion_E_Ext]
     type=Electrostatics
     #type=Diffusion
     variable=potential_ext
     block='1 2'
  [../]
  [./polar_electric_px]
     type=PolarElectricPStrong
     variable=polar_x
     permittivity_int = 8.85e-12
     permittivity_ext = 8.85e-12
     component=0
     implicit=false
  [../]
  [./polar_electric_py]
     type=PolarElectricPStrong
     variable=polar_y
     permittivity_int = 8.85e-12
     permittivity_ext = 8.85e-12
     component=1
     implicit=false
  [../]
  [./polar_electric_pz]
     type=PolarElectricPStrong
     variable=polar_z
     permittivity_int = 8.85e-12
     permittivity_ext = 8.85e-12
     component=2
     implicit=false
  [../]
  [./polar_x_time]
     type=TimeDerivative_scaled
     variable=polar_x
  [../]
  [./polar_y_time]
     type=TimeDerivative_scaled
     variable=polar_y
  [../]
  [./polar_z_time]
     type=TimeDerivative_scaled
     variable=polar_z
  [../]
[]

#[ICs]
#  [./polar_x_constic]
#     type=ConstantIC
#     variable=polar_x
#     value = 0.6
#  [../]
#  [./polar_y_constic]
#     type=ConstantIC
#     variable=polar_y
#     value = 0.6
#  [../]
#  [./polar_z_constic]
#     type=ConstantIC
#     variable=polar_z
#     value = 0.6
#  [../]
#[]

[BCs]

  [./potential_int_upz]
    type = DirichletBC 
    variable = potential_int
    boundary = '1'
    value = 0.0
    #implicit=false
  [../]
  [./potential_int_downz]
    type = DirichletBC 
    variable = potential_int
    boundary = '2'
    value = 0.0
    #implicit=false
  [../]
# Applied field: for zero field use NeumannBC on the external potential = 0. A
# Note that \nabla^2 \Phi_{ext} = 0 is satisfied if \Phi_{ext} = 0, ie Dirichlet and Neumann BC classes are equivalent
   [./potential_ext_upz]
    type = NeumannBC  
    variable = potential_ext
    boundary = '1'
    #value = 5.0e-11 #1e-12 gives a field of ~1e8 V/m ~ 1000 kV/cm for size 1, what about size 100?
                    # seems that this doesn't matter here? dt needs to change for the size however...
                    # lowering this value, to about 1e-16 is the limit here. 5e-11 gives 5.6e9 V/m for 100 size
    value = 0.0
    #implicit=false
  [../]
  [./potential_ext_downz]
    type = NeumannBC
    variable = potential_ext
    boundary = '2'
    #value = -5.0e-11
    value = 0.0
    #implicit=false
  [../]
[]

[Executioner]
  #type = Steady
  type=Transient
  solve_type=newton
  scheme = explicit-euler     #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dt=3e-15   #adjustable within a few orders of magnitude. It seems that we need small dt steps but large num_steps.
  dtmin=1e-19
  #l_max_its=1200
#NOTE: First time step calculates the depolarization field due to the unphysical initial condition. Energy may increase, which is allowed.
  #nl_max_its=5
  num_steps=2000
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname='-snes_rtol -ksp_type  -ksp_rtol -pc_type -snes_linesearch_type -pc_factor_zeropivot -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value='1e-6        gmres       1e-8      hypre       basic                1e-50                     0.25' #may need thresold to 0.5 for 3D problems?

#  [./TimeStepper]
#    type = TransientHalf
#    ratio = 0.5
#    min_dt = 1e-19
#    dtmin = 1e-19
#    dt = 1e-13
#  [../]
[]

[Postprocessors]
  [./bulk_energy]
    type=BulkEnergy
    initial_from_file_var = bulk_energy
    initial_from_file_timestep = 120
   [../]
   [./wall_energy]
    type=WallEnergy
    initial_from_file_var = wall_energy
    initial_from_file_timestep = 120
   [../]
   [./electric_energy]
    type=ElectrostaticEnergy
    initial_from_file_var = electric_energy
    initial_from_file_timestep = 120
    [../]
    [./total_energy]
    type=TotalEnergy
    bulk_energy=bulk_energy
    wall_energy=wall_energy
    electric_energy=electric_energy
    initial_from_file_var = total_energy
    initial_from_file_timestep = 120
    [../]
    [./gradx]
     type=TotalEnergyGradient
     component = 0
     permittivity_int = 8.85e-12
     permittivity_ext = 8.85e-12
     initial_from_file_var = gradx
     initial_from_file_timestep = 120
    [../]
    [./grady]
     type=TotalEnergyGradient
     component = 1
     permittivity_int = 8.85e-12
     permittivity_ext = 8.85e-12
     initial_from_file_var = grady
     initial_from_file_timestep = 120
    [../]
    [./gradz]
     type=TotalEnergyGradient
     component = 2
     permittivity_int = 8.85e-12
     permittivity_ext = 8.85e-12
     initial_from_file_var = gradz
     initial_from_file_timestep = 120
    [../]
    [./GradE-L2]
     type=TotalEnergyGradientL2
     gradx = gradx
     grady = grady
     gradz = gradz
     initial_from_file_var = GradE-L2
     initial_from_file_timestep = 120
    [../]
    [./R(i)]
     type=Residual
     block = '1 2'
     initial_from_file_var = R(i)
     initial_from_file_timestep = 120
    [../]
[]


[Outputs]

print_linear_residuals = true
print_perf_log = true

  [./out]
    type = Exodus
    file_base = out_PbTiO3_bulk_len-9_icConst_n3500_0field_1wall_adaptTime
    output_initial = true
    elemental_as_nodal = false
    interval = 1
  [../]

#  [./exodus_post]
#    type = Exodus
#    file_base = postprocess
#    interval = 1
#    output_postprocessors_on = 'timestep_end'
#  [../]

  [./debug]
    type = VariableResidualNormsDebugOutput
  [../]
[]
