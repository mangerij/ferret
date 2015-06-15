[Mesh]
  file = slab_exodus_coarse_100.e
  uniform_refine=0
[]

[Variables]
#active='polar_x polar_y polar_z'
#active='potential_int potential_ext'
  [./polar_x]
    #scaling=1e14
    order = FIRST
    family = LAGRANGE
    block='2'
  [../]
  [./polar_y]
    #scaling=1e14
    order = FIRST
    family = LAGRANGE
    block='2'
  [../]
  [./polar_z]
    #scaling=1e14
    order = FIRST
    family = LAGRANGE
    block='2'
  [../]
  [./potential_int]
    #scaling=1e7
    order=FIRST
    family = LAGRANGE
  [../]
  [./potential_ext]
    #scaling=1e7
    order=FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
[./auxv_es_energy_density_e] #es for electrostatic
     order=CONSTANT
     family=MONOMIAL
  [../]
  [./auxv_es_energy_density]
     order=CONSTANT
     family=MONOMIAL
  [../]
  [./auxv_bulk_energy_density]
     order=CONSTANT
     family=MONOMIAL
  [../]
  [./fieldE]
     order=CONSTANT
     family=MONOMIAL
  [../]
[]

[GlobalParams]

# alpha_{ijk} lifted from Chen appendix
   len_scale=1e-9
   #len_scale=1.0
   alpha1=-1.7252e8 # 3.8(T-479)*10^5 C^{-2}m^2
   alpha11=-7.3e7
   alpha111=2.6e8
   alpha12=7.5e8
   alpha112=6.1e8
   alpha123=-3.7e9
#   G110=4e-10
#   #G110=0.0
#   G11/G110=0.6
#   G12/G110=0.0
#   G44/G110=0.3
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
    variable = fieldE
    potential_ext = potential_ext
    potential_int = potential_int
  [../]
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
#  [./walled_x]
#     type=WallEnergyDerivative_scaled
#     variable=polar_x
#     component=0
#  [../]
#  [./walled_y]
#     type=WallEnergyDerivative_scaled
#     variable=polar_y
#     component=1
#  [../]
#  [./walled_z]
#     type=WallEnergyDerivative_scaled
#     variable=polar_z
#     component=2
#  [../]
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

[ICs]
  [./polar_x_constic]
     type=ConstantIC
     variable=polar_x
     value=0.6
  [../]
  [./polar_y_constic]
     type=ConstantIC
     variable=polar_y
     value=-0.6
  [../]
  [./polar_z_constic]
     type=ConstantIC
     variable=polar_z
     value=0.2
  [../]
[]

[BCs]

  [./potential_int_upz]
    type = DirichletBC
    variable = potential_int
    boundary = '1'
    value = 0.0
    #implicit=false
  [../]
  [./potential_int_downz]
    type = NeumannBC
    variable = potential_int
    boundary = '2'
    value = 0.0
    #implicit=false
  [../]
# Applied field: for zero field use DirichletBC on the external potential = 0. 
# Note that \nabla^2 \Phi_{ext} = 0 is satisfied if \Phi_{ext} = 0
   [./potential_ext_upz]
    type = NeumannBC
    variable = potential_ext
    boundary = '1'
    value = 0.0
    #value = 5.0e-11 #1e-12 gives a field of ~1e8 V/m ~ 1000 kV/cm for size 1, what about size 100?
                    # seems that this doesn't matter here? dt needs to change for the size however...
                    # lowering this value, to about 1e-16 is the limit here. 5e-11 gives 5.6e9 V/m for 100 size
    #implicit=false
  [../]
  [./potential_ext_downz]
    type = DirichletBC
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
  scheme=explicit-euler     #"implicit-euler, explicit-euler, crank-nicolson, bdf2, rk-2"
  dt=5e-14                  #adjustable within a few orders of magnitude due to scaling constant \tau
  l_max_its=1200
  nl_max_its=5
  num_steps=3500 
  #petsc_options="-snes_monitor -snes_converged_reason -ksp_monitor -ksp_converged_reason"
 # petsc_options='-snes_monitor -snes_converged_reason -ksp_monitor -ksp_converged_reason'
  petsc_options='-ksp_monitor_true_residual -snes_monitor -snes_view -snes_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname='-gmres_restart -snes_rtol -ksp_type  -ksp_rtol -pc_type -snes_linesearch_type -pc_factor_zeropivot'
  petsc_options_value='1000            1e-6        gmres       1e-8      hypre       basic                1e-50'
  #petsc_options_iname='-snes_rtol'
  #petsc_options_value='1e-16'
[]

[Postprocessors]
  [./bulk_energy]
    type=BulkEnergy
   [../]
#   [./wall_energy]
#    type=WallEnergy
#   [../]
   [./electric_energy]
    type=ElectrostaticEnergy
    [../]
#    [./total_energy]
#    type=TotalEnergy
#    bulk_energy=bulk_energy
#    wall_energy=wall_energy
#    electric_energy=electric_energy
#    [../]
[]


[Outputs]
  file_base = out_PbTiO3_bulk_len-9_icConst_dt1e17_n3500_coarse_Neumannfield-0
  output_initial = true
  print_linear_residuals = true
  print_perf_log = true
  [./out]
    type = Exodus
    elemental_as_nodal = true
    interval = 1
  [../]
  [./debug]
    type = VariableResidualNormsDebugOutput
  [../]
[]
